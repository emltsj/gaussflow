"""DeepSeek LLM integration: natural-language molecule generation.

Provider is fixed to DeepSeek (OpenAI-compatible ``/chat/completions``).
The model returns strict JSON with one or more molecules; every SMILES is
validated with RDKit before it is accepted, so the LLM can never inject an
invalid structure into the calculation queue.

Typical use: "生成 XXX 电解质分子，含 1 个 F" or, for incremental fluorine
scans, "F 原子个数 0-3 各生成一个" -> one molecule per F count.
"""
from __future__ import annotations

import json
import os
import re
import urllib.error
import urllib.request

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .structure import infer_charge_mult

DEEPSEEK_BASE = "https://api.deepseek.com"
# DeepSeek V4 ids.  "flash" is the fast/cheap default; "pro" is the strongest.
MODEL_FLASH = "deepseek-v4-flash"
MODEL_PRO = "deepseek-v4-pro"
MODEL_CHOICES = (MODEL_FLASH, MODEL_PRO)
MODEL_LABELS = {
    MODEL_FLASH: "V4 Flash（快 / 默认）",
    MODEL_PRO: "V4 Pro（最强 / 较慢）",
}
DEEPSEEK_MODEL = MODEL_FLASH

_MODEL_ALIASES = {
    "v4flash": MODEL_FLASH, "v4-flash": MODEL_FLASH,
    "deepseek-v4-flash": MODEL_FLASH, "deepseek-flash": MODEL_FLASH,
    "flash": MODEL_FLASH, "chat": MODEL_FLASH, "deepseek-chat": MODEL_FLASH,
    "v4pro": MODEL_PRO, "v4-pro": MODEL_PRO,
    "deepseek-v4-pro": MODEL_PRO, "deepseek-reasoner": MODEL_PRO,
    "pro": MODEL_PRO, "reasoner": MODEL_PRO,
}


class LlmError(RuntimeError):
    """Raised for configuration / API / parsing problems."""


def resolve_model(name: str | None = None, cfg=None) -> str:
    """Normalise a user-facing model name (v4pro / V4 Pro / …) to its API id.

    Unknown ids are passed through unchanged so private/compatible endpoints
    can still use their own model names.
    """
    raw = str(name or (cfg.get("llm.model") if cfg is not None else "")
              or MODEL_FLASH).strip()
    key = raw.lower().replace("_", "-").replace(" ", "")
    return _MODEL_ALIASES.get(key, raw)


def model_label(name: str | None = None) -> str:
    model = resolve_model(name)
    return MODEL_LABELS.get(model, model)


def model_from_label(label: str | None, cfg=None) -> str:
    """Reverse :func:`model_label` (accepts labels or raw ids)."""
    text = str(label or "").strip()
    for model, shown in MODEL_LABELS.items():
        if text == shown:
            return model
    return resolve_model(text or None, cfg)


# ---------------------------------------------------------------------------
_SYSTEM_PROMPT = (
    "你是计算化学分子建模助手。用户会用自然语言描述目标分子（例如某种电解质、"
    "含 F 的溶剂等）。你的任务是给出可直接用于量子化学计算的分子结构。\n"
    "只输出一个 JSON 对象，不要输出任何解释性文字或 Markdown 代码块，格式必须是：\n"
    '{"molecules":[{"name":"短名称",'
    '"smiles":"标准 SMILES",'
    '"charge":0,"multiplicity":1,'
    '"note":"一句话说明（例如 F 原子个数）"}],'
    '"explanation":"整体说明"}\n'
    "规则：\n"
    "1. SMILES 必须是 RDKit 能解析的合法结构；用标准有机 SMILES，不要写坐标。\n"
    "2. 默认中性闭壳层：charge=0, multiplicity=1；离子/自由基才标非零电荷和多重度。\n"
    "3. 如果用户要求多个分子（例如 F 原子个数 0,1,2,3），必须为每个数量各给一个分子，"
    "保持同一骨架，只改 F 的个数；name 用「基名-F0」「基名-F1」这种形式。\n"
    "4. 不要编造不存在的元素或价态；分子应尽量接近真实电解质/溶剂化学结构。\n"
    "5. molecules 数量不超过用户要求；宁可少而准确。"
)


def _api_key(cfg) -> str:
    key = (cfg.get("llm.api_key") or os.environ.get("DEEPSEEK_API_KEY") or "").strip()
    if not key:
        raise LlmError(
            "未配置 DeepSeek API Key。请在「设置 → DeepSeek 大模型」填写，"
            "或设置环境变量 DEEPSEEK_API_KEY。")
    return key


def chat(cfg, messages: list[dict], *, temperature: float | None = None,
         response_json: bool = True, timeout: float | None = None,
         thinking: str | None = None,
         reasoning_effort: str | None = None) -> dict:
    """Call DeepSeek chat/completions; returns {'content', 'usage', 'model'}."""
    key = _api_key(cfg)
    base = (cfg.get("llm.base_url") or DEEPSEEK_BASE).rstrip("/")
    url = base + "/chat/completions"
    think = str(thinking if thinking is not None
                else (cfg.get("llm.thinking") or "disabled")).lower()
    payload = {
        "model": resolve_model(cfg=cfg),
        "messages": messages,
        "stream": False,
        "thinking": {"type": "enabled" if think == "enabled" else "disabled"},
    }
    if think == "enabled":
        # DeepSeek V4 thinking mode: temperature is ignored; effort controls depth
        payload["reasoning_effort"] = str(
            reasoning_effort or cfg.get("llm.reasoning_effort") or "high")
    else:
        payload["temperature"] = float(
            temperature if temperature is not None
            else cfg.get("llm.temperature") or 0.2)
    if response_json:
        payload["response_format"] = {"type": "json_object"}
    req = urllib.request.Request(
        url, data=json.dumps(payload, ensure_ascii=False).encode("utf-8"),
        headers={"Content-Type": "application/json",
                 "Authorization": f"Bearer {key}"},
        method="POST")
    try:
        with urllib.request.urlopen(
                req, timeout=float(timeout or cfg.get("llm.timeout_seconds") or 180)) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
    except urllib.error.HTTPError as exc:
        detail = ""
        try:
            detail = exc.read().decode("utf-8", errors="replace")[:400]
        except Exception:
            pass
        raise LlmError(f"DeepSeek API HTTP {exc.code}: {detail or exc.reason}") from exc
    except urllib.error.URLError as exc:
        raise LlmError(f"无法连接 DeepSeek API: {exc.reason}") from exc
    except TimeoutError as exc:
        raise LlmError("DeepSeek API 超时") from exc
    try:
        message = data["choices"][0]["message"]
        content = message.get("content") or ""
    except (KeyError, IndexError, TypeError) as exc:
        raise LlmError(f"DeepSeek 返回格式异常: {str(data)[:300]}") from exc
    return {"content": content, "reasoning": message.get("reasoning_content"),
            "usage": data.get("usage"), "model": data.get("model"), "raw": data}


# ---------------------------------------------------------------------------
def parse_f_range(text: str | None) -> list[int]:
    """Parse '0-3', '0,1,2,3', '3' into a sorted unique list of F counts."""
    if text is None:
        return []
    text = str(text).strip()
    if not text:
        return []
    values: set[int] = set()
    for part in re.split(r"[,\s;，、]+", text):
        if not part:
            continue
        m = re.fullmatch(r"(\d+)\s*[-~]\s*(\d+)", part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            if a > b:
                a, b = b, a
            values.update(range(a, b + 1))
        elif part.isdigit():
            values.add(int(part))
    return sorted(v for v in values if 0 <= v <= 20)


def _extract_json(content: str) -> dict:
    text = (content or "").strip()
    text = re.sub(r"^```(?:json)?\s*", "", text)
    text = re.sub(r"\s*```$", "", text)
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass
    start, end = text.find("{"), text.rfind("}")
    if start >= 0 and end > start:
        try:
            return json.loads(text[start:end + 1])
        except json.JSONDecodeError:
            pass
    raise LlmError("大模型返回的内容不是合法 JSON，请重试或换一种描述。")


def validate_molecules(payload: dict, *, max_molecules: int = 12) -> dict:
    """Validate/normalise the LLM JSON payload into RDKit-checked molecules."""
    items = payload.get("molecules") if isinstance(payload, dict) else None
    if not isinstance(items, list) or not items:
        raise LlmError("大模型没有返回 molecules 列表。")
    out, errors, seen = [], [], set()
    for i, item in enumerate(items):
        if not isinstance(item, dict):
            errors.append(f"第 {i + 1} 项不是对象")
            continue
        smiles = str(item.get("smiles") or "").strip()
        name = str(item.get("name") or f"mol{i + 1}").strip()[:60] or f"mol{i + 1}"
        if not smiles:
            errors.append(f"{name}: 缺少 SMILES")
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            errors.append(f"{name}: SMILES 无法被 RDKit 解析")
            continue
        canonical = Chem.MolToSmiles(mol)
        if canonical in seen:
            continue
        seen.add(canonical)
        inferred_q, inferred_m = infer_charge_mult(mol)
        try:
            charge = int(item.get("charge", inferred_q))
        except (TypeError, ValueError):
            charge = inferred_q
        try:
            mult = int(item.get("multiplicity", inferred_m))
        except (TypeError, ValueError):
            mult = inferred_m
        if charge != inferred_q:
            # keep the user's/LLM's charge only if it makes electron parity sense
            electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
            if electrons % 2 == 0:
                mult = mult if mult % 2 == 1 else 1
            else:
                mult = mult if mult % 2 == 0 else 2
        n_f = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "F")
        out.append({
            "name": name,
            "smiles": canonical,
            "charge": charge,
            "multiplicity": mult,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "n_F": n_f,
            "natoms": mol.GetNumAtoms() + sum(a.GetTotalNumHs() for a in mol.GetAtoms()),
            "note": str(item.get("note") or "")[:200],
        })
        if len(out) >= max(1, int(max_molecules)):
            break
    if not out:
        raise LlmError("大模型给出的分子全部无法通过 RDKit 校验：" + "；".join(errors[:5]))
    return {"ok": True, "molecules": out, "errors": errors,
            "explanation": str(payload.get("explanation") or "")[:1000]}


def generate_molecules(cfg, prompt: str, f_range: str | None = None,
                       n_molecules: int = 1,
                       max_molecules: int | None = None,
                       extra_constraints: str | None = None) -> dict:
    """Natural language -> validated molecule list.

    * ``f_range`` given (e.g. "0-3"): one molecule per F count (same scaffold).
    * otherwise ``n_molecules`` > 1: that many *distinct* candidate molecules.
    """
    prompt = (prompt or "").strip()
    if not prompt:
        raise LlmError("请先输入分子描述。")
    counts = parse_f_range(f_range)
    try:
        n_want = max(1, int(n_molecules or 1))
    except (TypeError, ValueError) as exc:
        raise LlmError("生成分子数必须是整数") from exc
    limit = int(max_molecules or cfg.get("llm.max_molecules") or 12)
    limit = max(limit, len(counts) or n_want)
    user = prompt
    if counts:
        user += (
            f"\n\n要求：请针对 F（氟）原子个数 = {', '.join(map(str, counts))} "
            f"分别给出一个分子，共 {len(counts)} 个，保持同一骨架，只改变 F 的个数；"
            f"name 依次为「基名-F0」…；这是为了逐步增加 F 个数的系列计算。")
    elif n_want > 1:
        user += (
            f"\n\n要求：请给出 {n_want} 个互不相同的候选分子（同一类骨架的不同取代/异构体"
            f"或不同设计），每个都要有唯一 name（基名-1、基名-2…）和可直接解析的 SMILES；"
            f"这些将各自生成一个 gjf 并分别计算。")
    if extra_constraints:
        user += "\n补充要求：" + str(extra_constraints)
    user += f"\n最多返回 {limit} 个分子；只输出 JSON。"
    result = chat(cfg, [
        {"role": "system", "content": _SYSTEM_PROMPT},
        {"role": "user", "content": user},
    ])
    payload = _extract_json(result["content"])
    validated = validate_molecules(payload, max_molecules=limit)
    got = {m["n_F"] for m in validated["molecules"]}
    validated["missing_f_counts"] = [c for c in counts if c not in got]
    validated["requested_f_counts"] = counts
    validated["requested_n_molecules"] = n_want
    validated["prompt"] = prompt
    validated["model"] = result.get("model")
    validated["usage"] = result.get("usage")
    validated["raw"] = result["content"][:4000]
    return validated


_WORKFLOW_SYSTEM_PROMPT = """你是计算化学流程设计助手，为 Gaussian 16 生成「计算链」JSON。
只输出一个 JSON 对象，不要输出任何解释文字。结构：

{
  "recipe": {
    "name": "简短方案名（中文，如 自定义链-乙醇溶剂化）",
    "description": "一句话说明这个链的用途",
    "solvent": "DMSO" 或 null,
    "core_stages": [
      {"name": "optfreq", "kind": "optfreq",
       "route": "B3LYP/6-311+G(d) empiricaldispersion=gd3bj",
       "phase": "gas", "pop": null},
      {"name": "sp_gas", "kind": "sp",
       "route": "B2PLYP/def2TZVP empiricaldispersion=gd3bj",
       "phase": "gas", "pop": null},
      {"name": "sp_sol", "kind": "sp",
       "route": "M062X/6-31+G(d)", "phase": "solvent",
       "solvent": "DMSO", "pop": "full"}
    ],
    "li": {"on": true, "elements": ["O", "F"], "bsse_levels": ["optfreq", "sp_gas"]},
    "redox": {"on": false, "states": ["cation", "anion"], "targets": ["neutral"]}
  },
  "explanation": "为什么这样设计（简短中文）"
}

规则：
- kind 只能是 "optfreq" / "opt" / "freq" / "sp"（预优化由程序自动处理，不要写 preopt）；
  至少一个 optfreq 或 sp。
- phase 只能是 "gas" 或 "solvent"；溶剂相 stage 必须给 solvent（如 DMSO、水=Water、
  乙腈=Acetonitrile、乙醇=Ethanol、THF）。
- route 必须是合法的 Gaussian 路由行片段；方法/基组用 Gaussian 写法：
  B3LYP、M062X（不能写 M06-2X）、def2TZVP（不能写 def2-TZVP）、6-311+G(d)、
  empiricaldispersion=gd3bj、scrf 不要写（phase=solvent 时程序自动加 PCM）。
- 溶剂相 stage 可设 "pop": "full" 以输出更多轨道信息；否则 null。
- BSSE/counterpoise 只能用于气相 stage（Gaussian 不支持 PCM+counterpoise）。
- 若用户没提 Li+/氧化还原，li.on 与 redox.on 都给 false。
"""


def generate_workflow(cfg, prompt: str) -> dict:
    """Natural language -> calculation-chain recipe (unvalidated structure)."""
    prompt = (prompt or "").strip()
    if not prompt:
        raise LlmError("请先描述需要的计算链。")
    result = chat(cfg, [
        {"role": "system", "content": _WORKFLOW_SYSTEM_PROMPT},
        {"role": "user", "content": prompt + "\n\n只输出 JSON 对象。"},
    ], thinking="enabled", reasoning_effort="high")
    payload = _extract_json(result["content"])
    recipe = payload.get("recipe") if isinstance(payload.get("recipe"), dict) else None
    if recipe is None and isinstance(payload.get("core_stages"), list):
        recipe = payload  # model returned the recipe directly
    if not isinstance(recipe, dict):
        raise LlmError("大模型没有返回合法的计算链 JSON："
                       + str(payload)[:300])
    return {"ok": True, "recipe": recipe,
            "explanation": str(payload.get("explanation") or "")[:1200],
            "model": result.get("model"), "usage": result.get("usage"),
            "raw": result["content"][:4000]}
