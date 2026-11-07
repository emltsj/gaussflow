"""Central configuration for gaussflow.

Everything is overridable through ``config.json`` (or a custom path via the
``GFLOW_CONFIG`` environment variable).  The configuration is intentionally
small: it pins the Gaussian/GaussView binaries, the working root, and the
hard resource / approval limits that the Agent is *not* allowed to exceed.
"""
from __future__ import annotations

import copy
import json
import os
import platform
from pathlib import Path

DEFAULT_HOME = Path(__file__).resolve().parents[1]


def default_home() -> Path:
    """Resolve the data home *at call time* (packaged exe sets GFLOW_HOME)."""
    return Path(os.environ.get("GFLOW_HOME", str(DEFAULT_HOME)))


DEFAULTS: dict = {
    "config_version": 4,
    "home": str(DEFAULT_HOME),
    # --- external software -------------------------------------------------
    "g16_home": None,         # set by detect / GUI setting / G16_HOME
    "g16_exe": None,          # derived: <g16_home>/g16.exe
    "gview_exe": None,        # derived: <g16_home>/gview.exe
    "formchk_exe": None,      # derived: <g16_home>/formchk.exe
    # --- layout ------------------------------------------------------------
    "work_root": None,        # derived: <home>/runs
    "db_path": None,          # derived: <home>/gaussflow.db
    "log_dir": None,          # derived: <home>/logs
    # --- resources (hard ceiling; Agent may only lower, not raise) ---------
    "resources": {
        "nproc": 8,           # cores per Gaussian job
        "mem": "8GB",         # %mem per Gaussian job
        "max_nproc": 12,
        "max_mem_gb": 24,
    },
    # --- hard limits -------------------------------------------------------
    "limits": {
        "max_walltime_seconds": 7200,      # per job; killed afterwards
        "max_retries": 2,                  # whitelist-only automatic retries
        "max_jobs_per_batch": 60,
        "max_atoms": 200,
        "max_studies_per_batch": 60,
        "poll_seconds": 2.0,               # worker idle poll interval
        "inter_job_delay_seconds": 1.0,    # settle time between Gaussian jobs
    },
    # --- default science ---------------------------------------------------
    "chemistry": {
        "route": "B3LYP/6-31G(d)",         # main level of theory (no stage kw)
        "preopt_route": "HF/3-21G",        # cheap pre-optimisation level
        "preopt": True,                    # GaussView-Clean equivalent first pass
        "workflow": "opt_freq",            # opt | opt_freq | opt_freq_sp | sp
        "sp_route": None,                  # used by opt_freq_sp (default: route)
        "solvent": None,                   # e.g. "water" -> scrf=(smd,solvent=water)
        "extra_keywords": [],
        "use_chk_guess": True,             # freq/sp reuse opt chk via guess=read
        "charge": None,                    # None -> derive from input structure
        "multiplicity": None,              # None -> derive from electron parity
        "title": None,
    },
    # --- route whitelist / approval gate -----------------------------------
    "approval": {
        "allowed_methods": [
            "HF", "RHF", "UHF", "ROHF",
            "B3LYP", "B3PW91", "PBE0", "PBE", "PBEPBE", "BLYP", "BP86",
            "M06", "M062X", "M06-2X", "M06L", "M06-2X-D3",
            "WB97XD", "wB97XD", "wB97X-D", "CAM-B3LYP", "LC-WPBE",
            "B2PLYP", "B2PLYP-D3", "B2PLYP-D3BJ",
            "TPSSTPSS", "TPSS", "B97D", "B97D3", "B3LYP-D3", "B3LYP-D3BJ",
            "MP2", "CCSD", "CCSD(T)", "QCISD", "CIS", "TD-DFT",
        ],
        "allowed_basis": [
            "STO-3G", "3-21G", "6-31G", "6-31G(d)", "6-31G(d,p)",
            "6-31+G(d)", "6-31+G(d,p)", "6-31++G(d,p)", "6-311G(d,p)",
            "6-311+G(d)", "6-311+G(d,p)", "6-311+G(2d,p)",
            "6-311++G(d,p)", "6-311++G(2d,2p)",
            "cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "aug-cc-pVDZ", "aug-cc-pVTZ",
            "def2-SVP", "def2SVP", "def2-SVPP", "def2-TZVP", "def2TZVP",
            "def2-TZVPP", "def2QZVP", "def2-QZVP",
            "LANL2DZ", "SDD", "DGDZVP",
        ],
        # routes containing these need explicit human approval (expensive)
        "expensive_keywords": ["CCSD(T)", "QCISD(T)", "CBS-", "G4", "W1", "MP4"],
        "approval_atoms": 80,
        "approval_walltime_seconds": 14400,
    },
    # --- parsing / validation ---------------------------------------------
    "validation": {
        "require_scf_converged": True,
        "require_opt_converged": True,
        "require_normal_termination": True,
        "imaginary_freq_tolerance": 0,       # 0 for minima; 1 for TS mode
    },
    # --- preset calculation chains ----------------------------------------
    "recipes": {
        "opt_freq": {
            "name": "优化 + 频率",
            "description": "B3LYP/6-31G(d) 气相 opt+freq；开箱即用的通用起点",
            "solvent": None,
            "core_stages": [
                {"name": "optfreq", "kind": "optfreq",
                 "route": "B3LYP/6-31G(d)",
                 "phase": "gas"},
            ],
        },
        "standard": {
            "name": "通用三级流程（opt+freq / 气相 SP / 溶剂 SP）",
            "description": "B3LYP/6-31G(d)：opt+freq → 气相 SP → PCM(Water) SP；"
                           "可作为 Li⁺/BSSE/氧化还原分析的起点，也可用 AI 生成自己的链",
            "solvent": "Water",
            "core_stages": [
                {"name": "optfreq", "kind": "optfreq",
                 "route": "B3LYP/6-31G(d)", "phase": "gas"},
                {"name": "sp_gas", "kind": "sp",
                 "route": "B3LYP/6-31G(d)", "phase": "gas"},
                {"name": "sp_sol", "kind": "sp",
                 "route": "B3LYP/6-31G(d)", "phase": "solvent",
                 "solvent": "Water", "pop": "full"},
            ],
        },
    },
    # --- redox / thermochemistry references -------------------------------
    "redox": {
        "temperature_k": 298.15,
        # absolute potential scale: E(V) = -dG(F) + reference_offset_v
        "reference_offset_v": 0.0,
        "electron_free_energy_ev": 0.0,
    },
    "shermo": {
        "exe": None,
        "enabled": True,
        "command": "\"{shermo}\" \"{out}\"",
        "timeout_seconds": 300,
    },
    # --- DeepSeek LLM (molecule generation; provider fixed to DeepSeek) ----
    "llm": {
        "provider": "deepseek",
        "base_url": "https://api.deepseek.com",
        # V4 model ids: Flash is the fast/cheap default, Pro is the strongest.
        "model": "deepseek-v4-flash",
        "model_options": ["deepseek-v4-flash", "deepseek-v4-pro"],
        "thinking": "disabled",        # disabled | enabled(推理模式)
        "reasoning_effort": "high",    # low | high | max
        "api_key": "",
        "timeout_seconds": 180,
        "max_molecules": 12,
        "max_fluorination_variants": 20,
        "temperature": 0.2,
    },
}

_d = copy.deepcopy(DEFAULTS)


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


class Config:
    """Loaded configuration with derived paths."""

    def __init__(self, data: dict, path: Path):
        self.data = data
        self.path = Path(path)
        self.home = Path(data["home"]).resolve()

        work_root = data.get("work_root") or (self.home / "runs")
        self.work_root = Path(work_root).resolve()
        self.db_path = Path(data.get("db_path") or (self.home / "gaussflow.db")).resolve()
        self.log_dir = Path(data.get("log_dir") or (self.home / "logs")).resolve()
        self.exports_dir = self.home / "exports"

        g16_home_raw = data.get("g16_home")
        self.g16_home = Path(g16_home_raw).resolve() if g16_home_raw else None

        def _exe(key: str, filename: str):
            value = data.get(key)
            if value:
                return Path(value)
            return (self.g16_home / filename) if self.g16_home else None

        self.g16_exe = _exe("g16_exe", "g16.exe")
        self.gview_exe = _exe("gview_exe", "gview.exe")
        self.formchk_exe = _exe("formchk_exe", "formchk.exe")
        # Keep GAUSS_EXEDIR consistent with the configured g16.exe.
        if self.g16_home is None and self.g16_exe:
            self.g16_home = Path(self.g16_exe).resolve().parent

        self.ensure_dirs()

    # ------------------------------------------------------------------
    @classmethod
    def load(cls, path: str | os.PathLike | None = None) -> "Config":
        env_home = os.environ.get("GFLOW_HOME")
        if path is None:
            path = Path(os.environ.get("GFLOW_CONFIG")
                        or (default_home() / "config.json"))
        path = Path(path)
        exists = path.exists()
        migrated = False
        if exists:
            try:
                user = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:  # pragma: no cover
                raise RuntimeError(f"无法读取配置文件 {path}: {exc}") from exc
            # migrate older configs: built-in policy sections (method/basis
            # whitelist, recipes, validation rules) must follow the code, not a
            # stale list from an earlier version.
            if int(user.get("config_version") or 1) < int(DEFAULTS["config_version"]):
                for key in ("approval", "recipes", "validation", "shermo", "redox"):
                    user.pop(key, None)
                user["config_version"] = DEFAULTS["config_version"]
                migrated = True
            # upgrade retired model names to the current DeepSeek V4 ids
            llm = user.get("llm")
            if isinstance(llm, dict):
                if llm.get("model") in (None, "", "deepseek-chat",
                                        "deepseek-reasoner", "deepseek-flash"):
                    llm["model"] = "deepseek-v4-flash"
                    migrated = True
                if "thinking" not in llm:
                    llm["thinking"] = "disabled"
                    migrated = True
                options = llm.get("model_options")
                if (not options
                        or set(options) & {"deepseek-chat", "deepseek-reasoner",
                                           "deepseek-flash"}):
                    llm["model_options"] = ["deepseek-v4-flash", "deepseek-v4-pro"]
                    migrated = True
                user["llm"] = llm
        else:
            user = {}
        data = _deep_merge(DEFAULTS, user)
        # A config file defines its own home: this keeps temp/test configs and
        # custom GFLOW_CONFIG files self-contained instead of falling back to
        # the default project data directory.
        if not user.get("home"):
            data["home"] = str(path.resolve().parent)
        if env_home:
            data["home"] = env_home
        data["home"] = str(Path(data["home"]).resolve())
        cfg = cls(data, path)
        if not exists or migrated:
            cfg.save()
        return cfg

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(
            json.dumps(self.data, ensure_ascii=False, indent=2), encoding="utf-8"
        )

    def ensure_dirs(self) -> None:
        for d in (self.home, self.work_root, self.log_dir, self.exports_dir):
            d.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    def get(self, dotted: str, default=None):
        """Read a nested value with a dotted key, e.g. ``resources.nproc``."""
        cur = self.data
        for part in dotted.split("."):
            if not isinstance(cur, dict) or part not in cur:
                return default
            cur = cur[part]
        return cur

    def set(self, dotted: str, value) -> None:
        parts = dotted.split(".")
        cur = self.data
        for part in parts[:-1]:
            cur = cur.setdefault(part, {})
        cur[parts[-1]] = value

    # convenience typed views ------------------------------------------
    @property
    def resources(self) -> dict:
        return self.data["resources"]

    @property
    def limits(self) -> dict:
        return self.data["limits"]

    @property
    def chemistry(self) -> dict:
        return self.data["chemistry"]

    @property
    def approval(self) -> dict:
        return self.data["approval"]

    @property
    def validation(self) -> dict:
        return self.data["validation"]

    def require_binaries(self) -> list[str]:
        """Return a list of human-readable problems with configured binaries."""
        problems = []
        for label, path in (
            ("g16.exe", self.g16_exe),
            ("gview.exe", self.gview_exe),
        ):
            if not path or not Path(path).exists():
                problems.append(f"{label} 不存在: {path or '（未配置）'}")
        return problems

    def describe(self) -> dict:
        return {
            "home": str(self.home),
            "work_root": str(self.work_root),
            "db_path": str(self.db_path),
            "g16_exe": str(self.g16_exe),
            "gview_exe": str(self.gview_exe),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "resources": self.resources,
            "limits": self.limits,
        }



