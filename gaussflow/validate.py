"""Scientific acceptance rules.

A zero exit code or the presence of "Normal termination" is *not* enough.
Each stage is validated against the checks listed in the task specification:

* normal termination, no error termination, output complete
* SCF converged
* geometry optimisation really converged (4 criteria, stationary point)
* frequency job: no imaginary frequencies for a minimum
* charge / multiplicity / atom count consistency with the prepared input
* final energy taken from the final (highest-level) job of the chain
* thermochemistry present (for freq jobs) -- warning, not hard failure
* warnings are surfaced, not silently dropped
"""
from __future__ import annotations


def _check(cid: str, label: str, status: str, detail: str = "") -> dict:
    return {"id": cid, "label": label, "status": status, "detail": detail}


def validate(parsed: dict, job: dict, cfg) -> dict:
    """Return a validation report for one finished Gaussian job."""
    meta = job.get("meta") or {}
    stage = job.get("stage") or "sp"
    checks: list[dict] = []
    reasons: list[str] = []

    def add(cid, label, status, detail=""):
        checks.append(_check(cid, label, status, detail))
        if status == "fail":
            reasons.append(f"{label}: {detail}" if detail else label)

    # --- completeness --------------------------------------------------
    add("output_exists", "输出文件存在",
        "pass" if parsed.get("exists") else "fail",
        "" if parsed.get("exists") else "未找到 Gaussian 输出文件")

    add("normal_termination", "正常结束(Normal termination)",
        "pass" if parsed.get("normal_termination") else "fail",
        "" if parsed.get("normal_termination") else "输出中没有 Normal termination")

    add("no_error_termination", "无错误终止",
        "pass" if not parsed.get("error_termination") else "fail",
        "; ".join(parsed.get("error_lines") or [])[:300])

    # --- identity / consistency ---------------------------------------
    exp_natoms = meta.get("expected_natoms")
    if exp_natoms is not None:
        got = parsed.get("natoms")
        add("atom_count", "原子数一致",
            "pass" if got == exp_natoms else "fail",
            f"输入 {exp_natoms} / 输出 {got}")
    exp_charge = meta.get("expected_charge")
    exp_mult = meta.get("expected_multiplicity")
    if exp_charge is not None and parsed.get("charge") is not None:
        ok = parsed["charge"] == exp_charge
        add("charge", "电荷一致", "pass" if ok else "fail",
            f"输入 {exp_charge} / 输出 {parsed['charge']}")
    if exp_mult is not None and parsed.get("multiplicity") is not None:
        ok = parsed["multiplicity"] == exp_mult
        add("multiplicity", "自旋多重度一致", "pass" if ok else "fail",
            f"输入 {exp_mult} / 输出 {parsed['multiplicity']}")

    # --- SCF -----------------------------------------------------------
    scf_ok = bool(parsed.get("scf_converged")) and parsed.get("final_energy_hartree") is not None
    scf_status = "pass" if scf_ok else "fail"
    if cfg.validation.get("require_scf_converged", True):
        add("scf_converged", "SCF 收敛", scf_status,
            "" if scf_ok else "未找到收敛的 SCF Done 或出现 Convergence failure")
    else:
        add("scf_converged", "SCF 收敛", scf_status if scf_ok else "warn",
            "" if scf_ok else "SCF 未收敛（配置未强制）")

    # --- geometry optimisation ----------------------------------------
    opt = parsed.get("opt")
    expect_opt = meta.get("expect_opt_converged")
    if expect_opt is None:
        expect_opt = stage in ("preopt", "opt") or bool(parsed.get("is_opt"))
    if expect_opt and cfg.validation.get("require_opt_converged", True):
        if not parsed.get("is_opt"):
            add("opt_converged", "几何优化收敛", "fail", "输出中没有优化过程")
        elif not opt:
            add("opt_converged", "几何优化收敛", "fail", "优化信息缺失")
        else:
            ok = bool(opt.get("converged"))
            detail = ""
            if not ok:
                bad = [k for k, v in (opt.get("criteria") or {}).items()
                       if not v.get("converged")]
                detail = ("未达到收敛标准: " + ", ".join(bad)) if bad else \
                    f"completed={opt.get('completed')} stationary={opt.get('stationary_point')} stopped={opt.get('stopped')}"
            add("opt_converged", "几何优化收敛", "pass" if ok else "fail", detail)

    # --- frequencies ---------------------------------------------------
    expect_freq = meta.get("expect_freq")
    if expect_freq is None:
        expect_freq = stage == "freq"
    if expect_freq:
        freq = parsed.get("freq") or {}
        has = bool(freq.get("has_freq"))
        add("freq_present", "振动频率存在", "pass" if has else "fail",
            "" if has else "输出中没有频率信息")
        if has:
            tol = 1 if meta.get("ts_mode") else int(
                cfg.validation.get("imaginary_freq_tolerance", 0))
            n_imag = int(freq.get("n_imag") or 0)
            ok = n_imag <= tol
            detail = ""
            if not ok:
                detail = f"虚频 {n_imag} 个: {freq.get('imag_freqs')[:6]}"
            add("no_imaginary", "虚频检查", "pass" if ok else "fail", detail)
            if freq.get("zpve_hartree") is None and freq.get("enthalpy_hartree") is None:
                add("thermochem", "热化学数据", "warn", "未解析到 ZPE/焓/自由能")

    # --- BSSE / counterpoise ------------------------------------------
    if meta.get("expect_bsse"):
        cp = parsed.get("counterpoise") or {}
        ok = cp.get("bsse_energy_hartree") is not None and \
            cp.get("corrected_energy_hartree") is not None
        add("bsse", "BSSE/counterpoise 解析", "pass" if ok else "fail",
            "" if ok else "未解析到 Counterpoise BSSE/corrected energy")

    # --- energy --------------------------------------------------------
    e = parsed.get("final_energy_hartree")
    add("final_energy", "最终能量存在", "pass" if e is not None else "fail",
        "" if e is not None else "未提取到最终能量")

    # --- warnings ------------------------------------------------------
    warns = parsed.get("warnings") or []
    if warns:
        add("warnings", "Gaussian 警告", "warn", "; ".join(warns[:5])[:400])

    overall = "fail" if any(c["status"] == "fail" for c in checks) else (
        "warn" if any(c["status"] == "warn" for c in checks) else "pass")
    return {"overall": overall, "stage": stage, "checks": checks,
            "reasons": reasons}
