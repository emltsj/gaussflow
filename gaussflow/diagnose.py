"""Failure diagnosis + whitelist-only retry strategies.

The Agent is not allowed to invent arbitrary repairs.  Diagnosis maps a failed
job to one of a small set of known classes and, when the class is on the
whitelist and the retry budget is not exhausted, proposes a *specific, bounded*
retry (route patch, geometry restart, memory adjustment).  Everything else is
reported as ``needs_human`` with a short evidence excerpt instead of the full
log.
"""
from __future__ import annotations

import re
from pathlib import Path

# ordered: first match wins
_PATTERNS: list[tuple[str, str, str]] = [
    # (regex, class, human title)
    (r"Optimization stopped|Stationary point NOT found", "opt_not_converged",
     "几何优化未收敛"),
    (r"Convergence failure", "scf_not_converged", "SCF 收敛失败"),
    (r"out of memory|malloc|Out of memory|Insufficient memory", "memory_error",
     "内存不足"),
    (r"Erroneous write|No space left|Disk quota|FileIO|Unable to open|"
     r"read-write|Read error|write error", "io_error", "磁盘/文件读写错误"),
    (r"FormBX had a problem|Checkpoint|chk.*(corrupt|problem)|guess=read",
     "restart_chk_error", "断点(chk)/guess=read 读取失败"),
    (r"QPErr|Illegal IType|syntax error|End of file in|EOF|Bad data|Illegal|"
     r"Unrecognized|unknown keyword|Atomic number out of range|NAtoms",
     "input_error", "输入文件/坐标格式错误"),
    (r"Basis set .* not|basis.*not found|No basis", "input_error", "基组不可用"),
    (r"Error termination via Lnk1e", "gaussian_error", "Gaussian 错误终止"),
    (r"Error termination", "gaussian_error", "Gaussian 错误终止"),
]

_CLASS_TITLES = {
    "opt_not_converged": "几何优化未收敛",
    "scf_not_converged": "SCF 未收敛",
    "memory_error": "内存不足",
    "io_error": "磁盘/文件读写错误",
    "restart_chk_error": "断点文件/guess=read 读取失败",
    "input_error": "输入文件错误",
    "gaussian_error": "Gaussian 错误终止",
    "imaginary_frequencies": "存在虚频",
    "consistency_error": "电荷/自旋/原子数不一致",
    "timeout": "超过墙钟时间上限",
    "killed_or_truncated": "进程被终止或输出被截断",
    "validation_failed": "结果未通过验收规则",
    "none": "无错误",
}

_KNOWN_EVIDENCE = [
    r"Link\s+9999", r"Convergence failure", r"Error termination[^\n]*",
    r"out of memory[^\n]*", r"Erroneous write[^\n]*", r"FormBX[^\n]*",
    r"Atomic number out of range", r"End of file in[^\n]*",
    r"Unknown keyword[^\n]*", r"NAtoms[^\n]*", r"Problem[^\n]{0,120}",
]


def _evidence(parsed: dict, output_path: str | None, limit: int = 8) -> list[str]:
    found: list[str] = []
    for line in parsed.get("error_lines") or []:
        if line not in found:
            found.append(line)
    if output_path and Path(output_path).exists():
        text = Path(output_path).read_text(encoding="utf-8", errors="replace")
        tail = text[-12000:]
        for pat in _KNOWN_EVIDENCE:
            for m in re.finditer(pat, tail):
                ln = m.group(0).strip()
                if ln and ln not in found:
                    found.append(ln)
                if len(found) >= limit:
                    break
            if len(found) >= limit:
                break
        if not found:
            lines = [ln.strip() for ln in tail.splitlines() if ln.strip()]
            found = lines[-4:]
    return found[:limit]


def classify(parsed: dict, job: dict, validation: dict,
             output_path: str | None = None) -> str:
    if job.get("timed_out"):
        return "timeout"
    hay = "\n".join(parsed.get("error_lines") or [])
    if output_path and Path(output_path).exists():
        try:
            hay += "\n" + Path(output_path).read_text(
                encoding="utf-8", errors="replace")[-12000:]
        except OSError:
            pass
    if not parsed.get("normal_termination"):
        for pat, cls, _title in _PATTERNS:
            if re.search(pat, hay, re.I):
                return cls
        if parsed.get("exists"):
            return "killed_or_truncated"
        return "validation_failed"
    # normal termination but rejected by the validator
    if validation.get("overall") == "fail":
        for chk in validation.get("checks") or []:
            if chk.get("status") != "fail":
                continue
            cid = chk.get("id")
            if cid == "opt_converged":
                return "opt_not_converged"
            if cid == "scf_converged":
                return "scf_not_converged"
            if cid in ("charge", "multiplicity", "atom_count"):
                return "consistency_error"
            if cid == "no_imaginary":
                return "imaginary_frequencies"
        return "validation_failed"
    return "none"


def diagnose(parsed: dict, job: dict, validation: dict, cfg,
             output_path: str | None = None) -> dict:
    cls = classify(parsed, job, validation, output_path=output_path)
    attempt = int(job.get("attempt") or 1)
    max_retries = int(cfg.limits.get("max_retries", 2))
    meta = job.get("meta") or {}
    route = (job.get("route") or "").lower()

    retry = {"allowed": False, "strategy": None, "route_patch": None,
             "use_last_geometry": False, "mem": None, "restart_stage": None}
    suggestions: list[str] = []
    needs_human = cls not in ("opt_not_converged", "scf_not_converged",
                              "restart_chk_error", "memory_error",
                              "killed_or_truncated")

    if cls == "opt_not_converged" and attempt <= max_retries:
        patch = "opt=(maxcycles=200,calcfc)"
        if "maxcycles" in route:
            patch = "opt=(maxcycles=400,calcfc)"
        retry.update(allowed=True, strategy="restart_opt", route_patch=patch,
                     use_last_geometry=True)
        suggestions.append("用最后一步几何重启优化，并放宽最大步数/初始 Hessian")
        needs_human = False
    elif cls == "scf_not_converged" and attempt <= max_retries:
        if "xqc" not in route:
            retry.update(allowed=True, strategy="scf_boost",
                         route_patch="scf=(maxcycle=512,xqc)",
                         use_last_geometry=True)
            suggestions.append("提高 SCF 最大循环并启用 xqc 外推")
            needs_human = False
        else:
            suggestions.append("SCF 已使用 xqc 仍失败：需人工检查电荷/自旋/初始猜测")
    elif cls == "restart_chk_error":
        retry.update(allowed=True, strategy="no_guess", use_last_geometry=True,
                     route_patch=None)
        suggestions.append("放弃 guess=read，改用最后一步几何重新计算")
        needs_human = False
    elif cls == "memory_error" and attempt <= max_retries:
        retry.update(allowed=True, strategy="reduce_mem", use_last_geometry=True)
        suggestions.append("降低 %mem 后重试同一任务")
        needs_human = False
    elif cls == "killed_or_truncated" and attempt <= max_retries:
        # Windows-specific transient: Gaussian sometimes exits right after
        # Link 1 (AV / file-handle race).  Manual re-run of the identical
        # input succeeds, so one whitelist retry from the last geometry.
        retry.update(allowed=True, strategy="retry_transient",
                     use_last_geometry=True)
        suggestions.append("输出异常截断且无错误行（疑似启动/IO 瞬时问题），"
                           "使用最近几何重试一次")
        needs_human = False
    elif cls == "imaginary_frequencies":
        suggestions.append("对当前几何做更严格的优化 (opt=calcfc,tight) 后再算频率；"
                           "需要 Agent/人工决定是否追加计算")
        needs_human = True
    elif cls == "timeout":
        suggestions.append("任务超过墙钟上限：可提高 limits.max_walltime_seconds、"
                           "降低体系规模或更换方法，需人工批准")
        needs_human = True
    elif cls == "input_error":
        suggestions.append("输入文件/坐标/基组问题：修正输入后重新 prepare，系统不会自动重试")
        needs_human = True
    elif cls == "io_error":
        suggestions.append("检查磁盘空间与目录权限（含 %chk 路径），系统不会自动重试")
        needs_human = True
    elif cls == "consistency_error":
        suggestions.append("电荷/自旋多重度/原子数与输入不一致：需人工确认分子状态")
        needs_human = True
    elif cls == "killed_or_truncated":
        suggestions.append("输出被截断或进程被终止：检查是否被手动终止/内存不足")
        needs_human = True
    elif cls == "gaussian_error":
        suggestions.append("Gaussian 报错：请人工查看 evidence 中的错误行")
        needs_human = True
    elif cls == "validation_failed":
        suggestions.append("结果未通过验收规则，需要人工审核")
        needs_human = True

    if attempt > max_retries and cls in ("opt_not_converged", "scf_not_converged",
                                         "memory_error"):
        needs_human = True
        retry["allowed"] = False
        suggestions.append(f"已达到最大重试次数 ({max_retries})，转入人工处理")

    return {
        "job_id": job.get("id"),
        "study_id": job.get("study_id"),
        "stage": job.get("stage"),
        "class": cls,
        "title": _CLASS_TITLES.get(cls, cls),
        "summary": "; ".join(suggestions) if suggestions else "无异常",
        "evidence": _evidence(parsed, output_path, limit=8),
        "retry": retry,
        "needs_human": needs_human,
        "attempt": attempt,
        "max_retries": max_retries,
        "suggestions": suggestions,
        "validation_overall": validation.get("overall"),
        "reasons": validation.get("reasons") or [],
    }
