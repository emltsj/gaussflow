"""Pipeline orchestration: study creation, stage chaining, retries, summaries.

Stage graph (default ``opt_freq`` + preopt):

    preopt (HF/3-21G opt)  ->  opt (production level)  ->  freq (production level)
                                     |                          |
                                  job.chk                   %oldchk + guess=read

Every stage is a separate Gaussian job in its own directory; the chain is
advanced only after the previous stage passes validation.
"""
from __future__ import annotations

import json
import random
from pathlib import Path

from . import parse as parser
from .diagnose import diagnose
from .gjf import (build_counterpoise_gjf, build_gjf, compose_route,
                  describe_input, normalize_route, split_route,
                  strip_stage_keywords, validate_route)
from .store import Store, new_id, now
from .structure import (add_ion, coords_from_mol, formula_of_coords,
                        infer_charge_mult, list_coordination_sites,
                        molecule_from_file, molecule_from_smiles, xyz_block)
from .validate import validate

STAGE_SEQ = {"preopt": 0, "opt": 1, "freq": 2, "sp": 3}

WORKFLOWS = {
    "opt": ["opt"],
    "opt_freq": ["opt", "freq"],
    "opt_freq_sp": ["opt", "freq", "sp"],
    "sp": ["sp"],
    "freq": ["freq"],
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _study_dir(cfg, study_id: str) -> Path:
    d = cfg.work_root / study_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def _patch_meta(store: Store, study_id: str, patch: dict) -> dict:
    study = store.get_study(study_id)
    meta = dict(study.get("meta") or {})
    meta.update(patch)
    store.update_study(study_id, meta=meta)
    return meta


def _patch_job_meta(store: Store, job_id: str, patch: dict) -> dict:
    job = store.get_job(job_id)
    meta = dict(job.get("meta") or {})
    meta.update(patch)
    store.update_job(job_id, meta=meta)
    return meta


def stage_sequence(study: dict) -> list[str]:
    meta = study.get("meta") or {}
    specs = meta.get("stage_specs")
    if specs:
        return [s["name"] for s in specs]
    wf = study.get("workflow") or "opt_freq"
    seq = list(WORKFLOWS.get(wf, ["opt"]))
    if study.get("preopt") and seq and seq[0] == "opt":
        seq = ["preopt"] + seq
    return seq


def next_stage(study: dict, current: str) -> str | None:
    seq = stage_sequence(study)
    if current not in seq:
        return None
    i = seq.index(current)
    return seq[i + 1] if i + 1 < len(seq) else None


def _coords_from_geometry(geometry) -> list[tuple[str, float, float, float]] | None:
    if not geometry:
        return None
    try:
        return [(g["symbol"], float(g["x"]), float(g["y"]), float(g["z"]))
                for g in geometry]
    except (KeyError, TypeError, ValueError):
        return None


def _coords_from_study(study: dict) -> list[tuple[str, float, float, float]]:
    meta = study.get("meta") or {}
    coords = meta.get("initial_coords") or []
    return [(c["symbol"], float(c["x"]), float(c["y"]), float(c["z"])) for c in coords]


def _previous_stage_job(store: Store, study: dict, stage: str) -> dict | None:
    """Last *completed* job of the stage immediately preceding ``stage``."""
    seq = stage_sequence(study)
    if stage not in seq:
        return None
    idx = seq.index(stage)
    if idx == 0:
        return None
    prev_stage = seq[idx - 1]
    jobs = [j for j in store.study_jobs(study["id"])
            if j["stage"] == prev_stage and j["status"] == "completed"]
    if not jobs:
        return None
    return sorted(jobs, key=lambda j: (j.get("finished_at") or "", j["created_at"]))[-1]


def _result_of(store: Store, job: dict | None) -> dict:
    if not job:
        return {}
    res = store.get_result(job["id"]) or {}
    return res.get("parsed") or {}


def create_study(store: Store, cfg, *,
                 name: str,
                 smiles: str | None = None,
                 file: str | None = None,
                 route: str | None = None,
                 workflow: str | None = None,
                 charge: int | None = None,
                 multiplicity: int | None = None,
                 solvent: str | None = None,
                 preopt: bool | None = None,
                 title: str | None = None,
                 ion: str | None = None,
                 seed: int = 42,
                 n_conformers: int = 1,
                 sp_route: str | None = None,
                 batch_id: str | None = None,
                 variant: dict | None = None) -> dict:
    """Prepare one study; no Gaussian job is submitted here."""
    chem = cfg.chemistry
    if not smiles and not file:
        raise ValueError("必须提供 SMILES 或结构文件")
    if smiles and file:
        raise ValueError("SMILES 与结构文件只能二选一")

    if smiles:
        base_mol = molecule_from_smiles(smiles, seed=seed, n_conformers=n_conformers)
        src = {"smiles": smiles}
    else:
        base_mol, src = molecule_from_file(file, seed=seed)
        if base_mol.GetNumConformers() == 0:
            from rdkit import Chem
            base_mol = molecule_from_smiles(Chem.MolToSmiles(base_mol), seed=seed)

    # optional ion / variant
    ion_charge = 0
    mol = base_mol
    variant = dict(variant or {})
    if ion:
        mol = add_ion(base_mol, ion, seed=seed)
        ion_charge = _ion_charge(ion)
        variant.update({"type": "ion", "ion": ion})

    # charge / multiplicity
    inferred_charge, inferred_mult = infer_charge_mult(mol)
    if charge is None and multiplicity is None:
        final_charge, final_mult = inferred_charge, inferred_mult
    elif charge is not None and multiplicity is None:
        final_charge = int(charge) + ion_charge
        final_mult = 1 if (sum(a.GetAtomicNum() for a in mol.GetAtoms()) - final_charge) % 2 == 0 else 2
    elif charge is None and multiplicity is not None:
        final_charge, final_mult = inferred_charge, int(multiplicity)
    else:
        final_charge = int(charge) + ion_charge
        final_mult = int(multiplicity)

    coords = coords_from_mol(mol)
    natoms = len(coords)
    if natoms > int(cfg.limits.get("max_atoms", 200)):
        raise ValueError(f"原子数 {natoms} 超过上限 {cfg.limits['max_atoms']}")

    base_route = route or chem.get("route") or "B3LYP/6-31G(d)"
    raw_tokens = split_route(base_route)
    route_tokens, removed = strip_stage_keywords(raw_tokens)
    route_check = validate_route(route_tokens, cfg)
    sp_route_tokens = None
    if sp_route:
        sp_route_tokens, _removed = strip_stage_keywords(split_route(sp_route))
        sp_check = validate_route(sp_route_tokens, cfg)
        if sp_check["problems"]:
            raise ValueError(f"单点能路由不合法: {sp_check['problems']}")

    wf = workflow or chem.get("workflow") or "opt_freq"
    if wf not in WORKFLOWS:
        raise ValueError(f"未知 workflow: {wf}（可选 {', '.join(WORKFLOWS)}）")
    use_preopt = chem.get("preopt", True) if preopt is None else bool(preopt)
    if wf == "sp" or wf == "freq":
        use_preopt = False

    solv = solvent if solvent is not None else chem.get("solvent")

    approval_reasons: list[str] = []
    approval_reasons.extend(route_check["problems"])
    if route_check.get("needs_approval"):
        approval_reasons.extend(route_check.get("warnings") or [])
    if charge is not None or multiplicity is not None:
        if (final_charge != inferred_charge) or (final_mult != inferred_mult):
            approval_reasons.append(
                f"手动覆盖电荷/自旋多重度（推断 {inferred_charge}/{inferred_mult} → "
                f"指定 {final_charge}/{final_mult}）需要人工批准")
    if natoms > int(cfg.approval.get("approval_atoms", 80)):
        approval_reasons.append(f"原子数 {natoms} 大于 {cfg.approval['approval_atoms']}")
    route_upper = " ".join(route_tokens).upper()
    for kw in cfg.approval.get("expensive_keywords", []):
        if kw.upper() in route_upper:
            approval_reasons.append(f"昂贵方法关键字: {kw}")
    needs_approval = bool(approval_reasons) or not route_check["ok"]

    # first-approximation route preview
    preview = {}
    for st in ([s for s in ["preopt"] if use_preopt] + WORKFLOWS[wf]):
        toks = route_tokens
        preview[st] = compose_route(toks, st, solv,
                                    cfg.chemistry.get("extra_keywords") or [],
                                    None if st != "preopt" else None)
    if use_preopt:
        preview["preopt"] = compose_route(
            split_route(chem.get("preopt_route") or "HF/3-21G"), "opt", solv)

    study = store.create_study(
        name=name, smiles=smiles, charge=final_charge, multiplicity=final_mult,
        route=" ".join(route_tokens), workflow=wf, solvent=solv, preopt=use_preopt,
        status="needs_approval" if needs_approval else "prepared",
        spec={"source": src, "file": file, "route_input": base_route,
              "removed_stage_keywords": removed, "seed": seed,
              "n_conformers": n_conformers, "title": title or name},
        meta={
            "initial_coords": [{"symbol": s, "x": x, "y": y, "z": z}
                               for s, x, y, z in coords],
            "route_tokens": route_tokens,
            "route_input": base_route,
            "route_check": route_check,
            "sp_route_tokens": sp_route_tokens,
            "sp_route_input": sp_route,
            "workflow": wf,
            "preopt": use_preopt,
            "solvent": solv,
            "title": title or name,
            "natoms": natoms,
            "formula": formula_of_coords(coords),
            "variant": variant,
            "ion": ion,
            "batch_id": batch_id,
            "needs_approval": needs_approval,
            "approval_reasons": approval_reasons,
            "route_preview": preview,
            "created_at": now(),
        },
        batch_id=batch_id,
    )
    d = _study_dir(cfg, study["id"])
    (d / "initial.xyz").write_text(
        xyz_block(coords, f"{name} ({study['id']})"), encoding="utf-8")
    first_stage = (["preopt"] if use_preopt else []) + WORKFLOWS[wf][:1]
    first = first_stage[0]
    first_tokens = (split_route(chem.get("preopt_route") or "HF/3-21G")
                    if first == "preopt" else route_tokens)
    try:
        preview_route = compose_route(first_tokens, first, solv)
        (d / "initial.gjf").write_text(build_gjf(
            route_line=preview_route, title=title or name,
            charge=final_charge, multiplicity=final_mult, coords=coords,
            mem=cfg.resources["mem"], nproc=int(cfg.resources["nproc"]),
        ), encoding="ascii", errors="replace")
    except Exception:
        pass
    (d / "study.json").write_text(
        json.dumps({**study, "meta": study["meta"]}, ensure_ascii=False, indent=2,
                   default=str), encoding="utf-8")
    store.add_event(f"study prepared: {name} ({natoms} atoms, {wf})",
                    study_id=study["id"])
    return study


def _ion_charge(ion: str) -> int:
    from rdkit import Chem
    from .structure import COMMON_IONS
    smiles = COMMON_IONS.get(ion, ion)
    m = Chem.MolFromSmiles(smiles)
    return int(Chem.GetFormalCharge(m)) if m is not None else 0


# ---------------------------------------------------------------------------
# job creation
# ---------------------------------------------------------------------------
def create_job_for_stage(store: Store, cfg, study: dict, stage: str, *,
                         attempt: int = 1,
                         retry_of: str | None = None,
                         route_patch: str | None = None,
                         use_last_geometry: bool = False,
                         mem_override: str | None = None,
                         parent_job: str | None = None,
                         lane: int | None = None) -> dict:
    meta = study.get("meta") or {}
    if lane is None:
        lane = int(meta.get("lane") or 0)
    if meta.get("stage_specs"):
        return _create_recipe_job(
            store, cfg, study, stage, attempt=attempt, retry_of=retry_of,
            route_patch=route_patch, use_last_geometry=use_last_geometry,
            mem_override=mem_override, parent_job=parent_job, lane=lane)
    seq = STAGE_SEQ.get(stage, 9)
    study_dir = _study_dir(cfg, study["id"])
    job_dir = study_dir / f"{seq:02d}_{stage}_a{attempt}"
    n = 1
    while job_dir.exists():
        n += 1
        job_dir = study_dir / f"{seq:02d}_{stage}_a{attempt}_{n}"
    job_dir.mkdir(parents=True, exist_ok=True)

    # ---- geometry source ------------------------------------------------
    sources: list[str] = []
    if stage == "preopt":
        coords = _coords_from_study(study)
        sources.append("initial")
    elif use_last_geometry and retry_of:
        parsed = _result_of(store, store.get_job(retry_of))
        coords = _coords_from_geometry(parsed.get("final_geometry"))
        if coords:
            sources.append("failed_output_last_geometry")
        else:
            coords = None
    else:
        coords = None

    if coords is None and stage != "preopt":
        prev = _previous_stage_job(store, study, stage)
        coords = _coords_from_geometry(_result_of(store, prev).get("final_geometry"))
        if coords:
            sources.append(f"previous_stage:{prev['stage']}")
    if coords is None:
        coords = _coords_from_study(study)
        sources.append("initial")

    # ---- route ----------------------------------------------------------
    base_tokens = list(meta.get("route_tokens") or ["B3LYP/6-31G(d)"])
    stage_tokens = base_tokens
    if stage == "sp" and meta.get("sp_route_tokens"):
        stage_tokens = list(meta["sp_route_tokens"])
    if stage == "preopt":
        stage_tokens = split_route(cfg.chemistry.get("preopt_route") or "HF/3-21G")

    extra = list(cfg.chemistry.get("extra_keywords") or [])
    oldchk = None
    oldchk_src = None
    if stage in ("freq", "sp") and bool(cfg.chemistry.get("use_chk_guess", True)):
        prev = _previous_stage_job(store, study, stage)
        same_level = stage_tokens == base_tokens
        if prev and same_level and prev.get("chk_path") and Path(prev["chk_path"]).exists():
            oldchk_src = Path(prev["chk_path"])
            oldchk = "prev.chk"
            extra.append("guess=read")

    route_line = compose_route(
        stage_tokens, stage,
        meta.get("solvent"),
        extra,
        route_patch,
    )

    # ---- input ----------------------------------------------------------
    mem = mem_override or cfg.resources["mem"]
    nproc = int(cfg.resources["nproc"])
    input_text = build_gjf(
        route_line=route_line,
        title=meta.get("title") or study["name"],
        charge=int(study["charge"]),
        multiplicity=int(study["multiplicity"]),
        coords=coords,
        chk="job.chk",
        oldchk=oldchk,
        mem=mem,
        nproc=nproc,
        geom_check=False,
    )
    input_path = job_dir / "input.gjf"
    input_path.write_text(input_text, encoding="ascii", errors="replace")
    if oldchk_src:
        import shutil
        shutil.copy2(oldchk_src, job_dir / oldchk)

    job_meta = {
        "expected_natoms": len(coords),
        "expected_charge": int(study["charge"]),
        "expected_multiplicity": int(study["multiplicity"]),
        "geom_sources": sources,
        "oldchk": oldchk,
        "retry_route_patch": route_patch,
        "route_tokens": stage_tokens,
        "nproc": nproc,
        "mem": mem,
        "title": meta.get("title") or study["name"],
        "study_name": study["name"],
        "natoms": len(coords),
        "formula": formula_of_coords(coords),
        "attempt": attempt,
        "ts_mode": bool(meta.get("ts_mode")),
    }
    (job_dir / "meta.json").write_text(
        json.dumps(job_meta, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8")

    job = store.create_job(
        study_id=study["id"], stage=stage, seq=seq, attempt=attempt,
        status="pending", dir=str(job_dir), route=route_line,
        charge=int(study["charge"]), multiplicity=int(study["multiplicity"]),
        input_path=str(input_path), output_path=str(job_dir / "input.out"),
        chk_path=str(job_dir / "job.chk"), parent_job=parent_job,
        retry_of=retry_of, meta=job_meta, lane=lane,
    )
    store.add_event(
        f"job created: {stage} attempt {attempt} ({route_line})",
        study_id=study["id"], job_id=job["id"])
    return job


# ---------------------------------------------------------------------------
# recipe studies
# ---------------------------------------------------------------------------
def _compose_recipe_route(spec: dict, route_patch: str | None = None) -> str:
    """Build the final Gaussian route for one recipe stage spec."""
    tokens = split_route(spec.get("route") or "")
    kind = spec.get("kind")
    if kind == "optfreq":
        tokens += ["opt", "freq"]
    elif kind in ("opt", "preopt"):
        tokens += ["opt"]
    elif kind == "sp":
        pass
    elif kind == "bsse":
        tokens += ["counterpoise=2"]
    else:
        raise ValueError(f"未知 stage kind: {kind}")
    if spec.get("phase") == "solvent" and spec.get("solvent"):
        tokens.append(f"scrf=(pcm,solvent={spec['solvent']})")
    if spec.get("pop"):
        tokens.append(f"pop={spec['pop']}")
    if route_patch:
        low = route_patch.lower()
        if low.startswith(("opt=", "freq=")):
            # replace the stage keyword, e.g. opt -> opt=(maxcycles=200,calcfc)
            target = low.split("=")[0]
            for i, tok in enumerate(tokens):
                if tok.lower() == target:
                    tokens[i] = route_patch
                    break
            else:
                tokens.append(route_patch)
        else:
            tokens.append(route_patch)
    return "#p " + " ".join(tokens)


def _create_recipe_job(store: Store, cfg, study: dict, stage: str, *,
                       attempt: int = 1,
                       retry_of: str | None = None,
                       route_patch: str | None = None,
                       use_last_geometry: bool = False,
                       mem_override: str | None = None,
                       parent_job: str | None = None,
                       lane: int = 0) -> dict:
    meta = study.get("meta") or {}
    specs = meta.get("stage_specs") or []
    spec = next((s for s in specs if s.get("name") == stage), None)
    if spec is None:
        raise ValueError(f"stage '{stage}' 不在 study 的 stage_specs 中")
    seq = next(i for i, s in enumerate(specs) if s.get("name") == stage)
    study_dir = _study_dir(cfg, study["id"])
    job_dir = study_dir / f"{seq:02d}_{stage}_a{attempt}"
    n = 1
    while job_dir.exists():
        n += 1
        job_dir = study_dir / f"{seq:02d}_{stage}_a{attempt}_{n}"
    job_dir.mkdir(parents=True, exist_ok=True)

    # ---- geometry source ------------------------------------------------
    sources: list[str] = []
    coords = None
    if use_last_geometry and retry_of:
        parsed = _result_of(store, store.get_job(retry_of))
        coords = _coords_from_geometry(parsed.get("final_geometry"))
        if coords:
            sources.append("failed_output_last_geometry")
    if coords is None:
        prev = _previous_stage_job(store, study, stage)
        coords = _coords_from_geometry(_result_of(store, prev).get("final_geometry"))
        if coords:
            sources.append(f"previous_stage:{prev['stage']}")
    if coords is None:
        coords = _coords_from_study(study)
        sources.append("initial")

    mem = mem_override or cfg.resources["mem"]
    nproc = int(cfg.resources["nproc"])
    route_line = _compose_recipe_route(spec, route_patch)
    title = f"{meta.get('title') or study['name']} | {stage}"

    if spec.get("kind") == "bsse":
        fragments = spec.get("fragments") or meta.get("fragments")
        if not fragments:
            raise ValueError("BSSE 任务缺少 fragment 定义")
        if len(coords) != sum(int(f["n"]) for f in fragments):
            raise ValueError(
                f"BSSE 原子数不匹配: 坐标 {len(coords)}, 片段 "
                f"{[f['n'] for f in fragments]}")
        input_text = build_counterpoise_gjf(
            route_line=route_line, title=title, coords=coords,
            fragments=fragments,
            total_charge=int(study["charge"]),
            total_multiplicity=int(study["multiplicity"]),
            chk="job.chk", mem=mem, nproc=nproc)
        expected_charge = None
        expected_mult = None
    else:
        input_text = build_gjf(
            route_line=route_line, title=title,
            charge=int(study["charge"]), multiplicity=int(study["multiplicity"]),
            coords=coords, chk="job.chk", mem=mem, nproc=nproc)
        expected_charge = int(study["charge"])
        expected_mult = int(study["multiplicity"])

    input_path = job_dir / "input.gjf"
    input_path.write_text(input_text, encoding="ascii", errors="replace")
    job_meta = {
        "expected_natoms": len(coords),
        "expected_charge": expected_charge,
        "expected_multiplicity": expected_mult,
        "geom_sources": sources,
        "kind": spec.get("kind"),
        "stage_spec": spec,
        "expect_opt_converged": spec.get("kind") in ("optfreq", "opt"),
        "expect_freq": spec.get("kind") == "optfreq",
        "expect_bsse": spec.get("kind") == "bsse",
        "phase": spec.get("phase"),
        "solvent": spec.get("solvent"),
        "retry_route_patch": route_patch,
        "route_tokens": split_route(spec.get("route") or ""),
        "nproc": nproc, "mem": mem,
        "title": title, "study_name": study["name"],
        "natoms": len(coords),
        "formula": formula_of_coords(coords),
        "attempt": attempt, "ts_mode": False,
    }
    (job_dir / "meta.json").write_text(
        json.dumps(job_meta, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8")
    job = store.create_job(
        study_id=study["id"], stage=stage, seq=seq, attempt=attempt,
        status="pending", dir=str(job_dir), route=route_line,
        charge=int(study["charge"]) if expected_charge is not None else None,
        multiplicity=int(study["multiplicity"]) if expected_mult is not None else None,
        input_path=str(input_path), output_path=str(job_dir / "input.out"),
        chk_path=str(job_dir / "job.chk"), parent_job=parent_job,
        retry_of=retry_of, meta=job_meta, lane=lane)
    store.add_event(f"job created: {stage} attempt {attempt} ({route_line})",
                    study_id=study["id"], job_id=job["id"])
    return job


def create_recipe_system(store: Store, cfg, *,
                         name: str,
                         smiles: str | None = None,
                         file: str | None = None,
                         recipe_name: str = "opt_freq",
                         options: dict | None = None,
                         seed: int = 42,
                         n_conformers: int = 1,
                         charge: int | None = None,
                         multiplicity: int | None = None) -> dict:
    """Create the whole calculation system (batch of studies)."""
    from .recipes import plan_studies

    if not smiles and not file:
        raise ValueError("必须提供 SMILES 或结构文件")
    if smiles:
        base_mol = molecule_from_smiles(smiles, seed=seed, n_conformers=n_conformers)
        src = {"smiles": smiles}
    else:
        base_mol, src = molecule_from_file(file, seed=seed)
        if base_mol.GetNumConformers() == 0:
            from rdkit import Chem
            base_mol = molecule_from_smiles(Chem.MolToSmiles(base_mol), seed=seed)

    inf_charge, inf_mult = infer_charge_mult(base_mol)
    base_charge = inf_charge if charge is None else int(charge)
    base_mult = inf_mult if multiplicity is None else int(multiplicity)
    total_z = sum(a.GetAtomicNum() for a in base_mol.GetAtoms())
    natoms = base_mol.GetNumAtoms()
    if natoms > int(cfg.limits.get("max_atoms", 200)):
        raise ValueError(f"原子数 {natoms} 超过上限 {cfg.limits['max_atoms']}")

    li_cfg = (cfg.get("recipes") or {}).get(recipe_name, {}).get("li") or {}
    elements = list((options or {}).get("li_elements") or li_cfg.get("elements", ["O", "F"]))
    sites = list_coordination_sites(base_mol, elements)

    plan = plan_studies(
        cfg, recipe_name, name=name, charge=base_charge,
        multiplicity=base_mult, total_z=total_z, natoms=natoms,
        sites=sites, options=options)

    batch_id = new_id("bt_")
    created: list[dict] = []
    for item in plan["studies"]:
        role = item["role"]
        li_site = item.get("li_site")
        if role == "li_ref":
            # bare Li+ reference: single atom at the origin
            from .structure import mol_from_coords
            coords = [("Li", 0.0, 0.0, 0.0)]
            mol = mol_from_coords(coords, charge=1)
        else:
            mol = base_mol
            if li_site:
                mol = add_ion(base_mol, "Li+", anchor_index=li_site["atom_index"],
                              dist=float(li_site.get("dist") or 2.0), seed=seed)
            coords = coords_from_mol(mol)
        stage_specs = item["stage_specs"]
        # validate every route of the recipe once
        problems, warnings = [], []
        for spec in stage_specs:
            chk = validate_route(split_route(spec.get("route") or ""), cfg)
            for item in (chk.get("problems") or []):
                if item not in problems:
                    problems.append(item)
            for item in (chk.get("warnings") or []):
                if item not in warnings:
                    warnings.append(item)
        needs_approval = bool(problems)
        if natoms > int(cfg.approval.get("approval_atoms", 80)):
            needs_approval = True
        if charge is not None and base_charge != inf_charge:
            needs_approval = True
        approval_reasons = list(problems)
        if natoms > int(cfg.approval.get("approval_atoms", 80)):
            approval_reasons.append(f"原子数 {natoms} 大于 {cfg.approval['approval_atoms']}")

        meta = {
            "initial_coords": [{"symbol": s, "x": x, "y": y, "z": z}
                               for s, x, y, z in coords],
            "natoms": len(coords),
            "formula": formula_of_coords(coords),
            "recipe": recipe_name,
            "recipe_name": plan["recipe_name"],
            "role": role,
            "stage_specs": stage_specs,
            "stage_names": [s["name"] for s in stage_specs],
            "li_site": li_site,
            "fragments": item.get("fragments"),
            "redox": item.get("redox"),
            "variant": item.get("variant"),
            "title": item["name"],
            "batch_id": batch_id,
            "base_charge": base_charge,
            "base_multiplicity": base_mult,
            "total_z": total_z,
            "needs_approval": needs_approval,
            "approval_reasons": approval_reasons,
            "route_check": {"problems": problems, "warnings": warnings},
            "notes": plan.get("notes") or [],
            "solvent": plan.get("solvent"),
        }
        study = store.create_study(
            name=item["name"], smiles=smiles, charge=item["charge"],
            multiplicity=item["multiplicity"],
            route=stage_specs[0].get("route"),
            workflow=f"recipe:{recipe_name}", solvent=plan.get("solvent"),
            preopt=False,
            status="needs_approval" if needs_approval else "prepared",
            spec={"source": src, "recipe": recipe_name, "role": role,
                  "variant": item.get("variant")},
            meta=meta, batch_id=batch_id)
        d = _study_dir(cfg, study["id"])
        (d / "initial.xyz").write_text(
            xyz_block(coords, f"{study['name']} ({study['id']})"), encoding="utf-8")
        try:
            first_route = _compose_recipe_route(stage_specs[0])
            (d / "initial.gjf").write_text(build_gjf(
                route_line=first_route, title=item["name"],
                charge=int(item["charge"]), multiplicity=int(item["multiplicity"]),
                coords=coords, mem=cfg.resources["mem"],
                nproc=int(cfg.resources["nproc"])), encoding="ascii", errors="replace")
        except Exception:
            pass
        created.append(study)
    store.add_event(
        f"recipe system created: {name} / {recipe_name} "
        f"({len(created)} studies, {len(sites)} sites)",
        study_id=None)
    return {
        "ok": True,
        "batch_id": batch_id,
        "name": name,
        "recipe": recipe_name,
        "recipe_name": plan["recipe_name"],
        "natoms": natoms,
        "base_charge": base_charge,
        "base_multiplicity": base_mult,
        "li_sites": sites,
        "n_studies": len(created),
        "studies": [{"study_id": s["id"], "name": s["name"],
                     "role": (s.get("meta") or {}).get("role"),
                     "charge": s["charge"], "multiplicity": s["multiplicity"],
                     "status": s["status"],
                     "stages": (s.get("meta") or {}).get("stage_names")}
                    for s in created],
    }


def submit_study(store: Store, cfg, study_id: str, approve: bool = False) -> dict:
    """Submit a study's first-stage job (asynchronous)."""
    study = store.get_study(study_id)
    if study is None:
        return {"ok": False, "error": f"未知 study: {study_id}"}
    if study["status"] in ("submitted", "running", "completed"):
        return {"ok": False, "error": f"study 状态为 {study['status']}，不能重复提交"}
    meta = study.get("meta") or {}
    if meta.get("needs_approval") and not approve:
        store.update_study(study_id, status="needs_approval")
        return {"ok": True, "status": "needs_approval", "study_id": study_id,
                "reasons": meta.get("approval_reasons") or [],
                "message": "该任务需要 --approve 显式批准后才能提交"}
    seq = stage_sequence(study)
    first = seq[0]
    lane = store.next_lane()
    _patch_meta(store, study_id, {"lane": lane})
    job = create_job_for_stage(store, cfg, study, first, lane=lane)
    store.update_study(study_id, status="submitted")
    store.add_event(f"study submitted (first stage {first}, lane {lane})",
                    study_id=study_id)
    return {"ok": True, "status": "submitted", "study_id": study_id,
            "jobs": [job["id"]], "stages": seq}


# ---------------------------------------------------------------------------
# input from existing .gjf files (multi-select) and AI-generated molecules
# ---------------------------------------------------------------------------
def gjf_study_spec(path) -> dict:
    """Parse an existing .gjf into the fields needed to queue it as a study."""
    from .structure import read_gjf
    p = Path(path)
    parsed = read_gjf(p)
    tokens = split_route(parsed["route"])
    base_tokens, removed = strip_stage_keywords(tokens)
    low = " ".join(t.lower() for t in removed)
    if "opt" in low and "freq" in low:
        workflow = "opt_freq"
    elif "opt" in low:
        workflow = "opt"
    elif "freq" in low:
        workflow = "freq"
    else:
        workflow = "sp"
    return {
        "path": str(p),
        "title": parsed["title"],
        "route_raw": parsed["route"],
        "route_tokens": base_tokens,
        "removed_stage_keywords": removed,
        "charge": int(parsed["charge"]),
        "multiplicity": int(parsed["multiplicity"]),
        "coords": parsed["coords"],
        "natoms": len(parsed["coords"]),
        "formula": formula_of_coords(parsed["coords"]),
        "workflow": workflow,
    }


def create_batch_from_gjfs(store: Store, cfg, paths: list[str], *,
                           name_prefix: str | None = None,
                           submit: bool = False,
                           approve: bool = False) -> dict:
    """Queue one study per selected .gjf file (multi-select support)."""
    paths = [str(p) for p in (paths or []) if str(p).strip()]
    if not paths:
        raise ValueError("没有选择 .gjf 文件")
    if len(paths) > int(cfg.limits.get("max_jobs_per_batch", 60)):
        raise ValueError(f"一次最多 {cfg.limits['max_jobs_per_batch']} 个 gjf 文件")
    batch_id = new_id("bt_")
    studies, specs, used_names = [], [], set()
    for idx, p in enumerate(paths):
        spec = gjf_study_spec(p)
        base = Path(p).stem or f"gjf{idx + 1}"
        name = f"{name_prefix}-{base}" if name_prefix else base
        n = 2
        while name in used_names:
            name = f"{base}-{n}"
            n += 1
        used_names.add(name)
        route_check = validate_route(spec["route_tokens"], cfg)
        needs_approval = bool(route_check["problems"]) or \
            spec["natoms"] > int(cfg.approval.get("approval_atoms", 80))
        study = store.create_study(
            name=name, smiles=None, charge=spec["charge"],
            multiplicity=spec["multiplicity"],
            route=" ".join(spec["route_tokens"]),
            workflow=spec["workflow"], solvent=None, preopt=False,
            status="needs_approval" if needs_approval else "prepared",
            spec={"source": {"type": "gjf", "path": str(p),
                             "raw_route": spec["route_raw"]}},
            meta={
                "initial_coords": [{"symbol": s, "x": x, "y": y, "z": z}
                                   for s, x, y, z in spec["coords"]],
                "natoms": spec["natoms"],
                "formula": spec["formula"],
                "title": spec["title"] or name,
                "route_tokens": spec["route_tokens"],
                "removed_stage_keywords": spec["removed_stage_keywords"],
                "workflow": spec["workflow"],
                "preopt": False,
                "variant": {"type": "gjf", "label": Path(p).stem},
                "needs_approval": needs_approval,
                "approval_reasons": list(route_check["problems"]),
                "route_check": route_check,
                "notes": [],
                "batch_id": batch_id,
            },
            batch_id=batch_id)
        d = _study_dir(cfg, study["id"])
        (d / "initial.xyz").write_text(
            xyz_block(spec["coords"], f"{name} ({study['id']})"), encoding="utf-8")
        try:
            (d / "initial.gjf").write_text(
                Path(p).read_text(encoding="utf-8", errors="replace"), encoding="utf-8")
        except OSError:
            pass
        studies.append(study)
        specs.append({k: v for k, v in spec.items() if k != "coords"})
    if submit:
        for s in studies:
            submit_study(store, cfg, s["id"], approve=approve)
    store.add_event(f"batch created from {len(studies)} gjf files", study_id=None)
    return {
        "ok": True,
        "batch_id": batch_id,
        "n_studies": len(studies),
        "studies": [{"study_id": s["id"], "name": s["name"],
                     "charge": s["charge"], "multiplicity": s["multiplicity"],
                     "workflow": (s.get("meta") or {}).get("workflow"),
                     "formula": (s.get("meta") or {}).get("formula"),
                     "status": s["status"]} for s in studies],
        "specs": specs,
        "submitted": bool(submit),
    }


def export_gjf_previews(cfg, molecules: list[dict], *,
                        recipe_name: str = "opt_freq",
                        name_prefix: str = "AI") -> dict:
    """Write one .gjf per LLM-generated molecule (preview / manual reuse)."""
    import re as _re

    from .recipes import get_recipe
    from .structure import molecule_from_smiles

    recipe = get_recipe(cfg, recipe_name)
    core = recipe.get("core_stages") or [{"route": "B3LYP/6-311+G(d)",
                                          "kind": "optfreq"}]
    first = core[0]
    tokens = split_route(first.get("route") or "B3LYP/6-311+G(d)")
    stage_kw = ""
    if first.get("kind") == "optfreq":
        stage_kw = " opt freq"
    elif first.get("kind") in ("opt", "preopt"):
        stage_kw = " opt"
    elif first.get("kind") == "freq":
        stage_kw = " freq"
    if first.get("phase") == "solvent" and first.get("solvent"):
        stage_kw += f" scrf=(pcm,solvent={first['solvent']})"
    route_line = "#p " + " ".join(tokens) + stage_kw
    out_dir = cfg.home / "generated_gjf" / new_id("gen_")
    out_dir.mkdir(parents=True, exist_ok=True)
    files = []
    for i, m in enumerate(molecules):
        mol = molecule_from_smiles(m["smiles"], seed=42 + i)
        coords = coords_from_mol(mol)
        title = f"{name_prefix}-{m.get('name') or i + 1}"
        text = build_gjf(route_line=route_line, title=title,
                         charge=int(m.get("charge") or 0),
                         multiplicity=int(m.get("multiplicity") or 1),
                         coords=coords, mem=cfg.resources["mem"],
                         nproc=int(cfg.resources["nproc"]))
        safe = _re.sub(r"[^0-9A-Za-z_.\-]+", "_", f"{i + 1:02d}_{m.get('name') or 'mol'}")
        path = out_dir / f"{safe}.gjf"
        path.write_text(text, encoding="ascii", errors="replace")
        files.append({
            "name": title, "smiles": m["smiles"],
            "formula": m.get("formula"), "n_F": m.get("n_F"),
            "charge": m.get("charge"), "multiplicity": m.get("multiplicity"),
            "gjf": str(path),
        })
    (out_dir / "index.json").write_text(
        json.dumps({"molecules": molecules, "files": files, "route": route_line},
                   ensure_ascii=False, indent=2), encoding="utf-8")
    return {"ok": True, "gjf_dir": str(out_dir), "files": files,
            "route": route_line, "index": str(out_dir / "index.json")}


def create_ai_systems(store: Store, cfg, molecules: list[dict], *,
                      recipe_name: str = "opt_freq",
                      options: dict | None = None,
                      name_prefix: str = "AI",
                      submit: bool = True,
                      approve: bool = False) -> dict:
    """Create one calculation system per AI-generated molecule, then queue."""
    if not molecules:
        raise ValueError("没有可计算的分子")
    batches = []
    for i, m in enumerate(molecules):
        name = f"{name_prefix}-{m.get('name') or i + 1}"
        res = create_recipe_system(
            store, cfg, name=name, smiles=m["smiles"],
            recipe_name=recipe_name, options=options,
            charge=m.get("charge"), multiplicity=m.get("multiplicity"))
        item = {
            "batch_id": res["batch_id"],
            "name": name,
            "smiles": m["smiles"],
            "formula": m.get("formula"),
            "n_F": m.get("n_F"),
            "studies": [s["study_id"] for s in res["studies"]],
        }
        if submit:
            item["submit"] = submit_batch(store, cfg, res["batch_id"], approve=approve)
        batches.append(item)
    return {"ok": True, "n_molecules": len(molecules),
            "batches": batches, "submitted": bool(submit)}


# ---------------------------------------------------------------------------
# user cancellation (force stop)
# ---------------------------------------------------------------------------
def cancel_jobs(store: Store, cfg, *, batch_ids: list[str] | None = None,
                all_jobs: bool = False) -> dict:
    """Force-stop running Gaussian jobs and cancel the queue.

    Running jobs: kill the g16.exe process tree immediately and mark the job
    ``cancelled`` (no automatic retry).  Pending jobs: mark ``cancelled``.
    Studies in the affected scope are marked ``cancelled`` as well.
    """
    from .runner import kill_pid, request_cancel

    scope: set[str] | None = None
    if batch_ids and not all_jobs:
        scope = set()
        for bid in batch_ids:
            scope |= {s["id"] for s in store.list_studies(batch_id=bid, limit=2000)}

    jobs = store.list_jobs(limit=5000)
    running, pending = [], []
    for job in jobs:
        if job["status"] not in ("running", "pending"):
            continue
        if scope is not None and job["study_id"] not in scope:
            continue
        if job["status"] == "pending":
            store.update_job(job["id"], status="cancelled",
                             error_class="stopped_by_user", finished_at=now())
            pending.append(job["id"])
            continue
        pid = job.get("pid")
        killed = request_cancel(job["id"])
        if pid and int(pid) > 0 and not killed:
            kill_pid(int(pid))
        store.update_job(job["id"], status="cancelled", pid=None,
                         error_class="stopped_by_user", finished_at=now())
        running.append(job["id"])
        store.add_event(f"force-killed Gaussian pid={pid} for {job['stage']}",
                        study_id=job["study_id"], job_id=job["id"], level="warning")

    study_ids = scope if scope is not None else {j["study_id"] for j in jobs}
    for sid in study_ids:
        s = store.get_study(sid)
        if s and s["status"] in ("prepared", "submitted", "running",
                                 "needs_approval"):
            store.update_study(sid, status="cancelled")

    store.add_event(
        f"user stop: {len(running)} running + {len(pending)} pending cancelled",
        level="warning")
    return {"ok": True,
            "running_cancelled": running,
            "pending_cancelled": pending,
            "cancelled_count": len(running) + len(pending)}


def purge_cancelled(store: Store, cfg, *, batch_ids: list[str] | None = None,
                    delete_files: bool = False) -> dict:
    """Remove cancelled studies (and their jobs/results) from the task list.

    By default the files under ``runs/`` are kept as audit records; pass
    ``delete_files=True`` to delete the cancelled job directories too.
    """
    import shutil

    studies = store.list_studies(limit=5000)
    scope = set(batch_ids) if batch_ids else None
    removed = []
    for s in studies:
        if s["status"] != "cancelled":
            continue
        if scope is not None and s.get("batch_id") not in scope:
            continue
        sid = s["id"]
        for table in ("results", "jobs", "events"):
            store.conn.execute(f"delete from {table} where study_id=?", (sid,))
        store.conn.execute("delete from studies where id=?", (sid,))
        if delete_files:
            d = cfg.work_root / sid
            if d.exists():
                shutil.rmtree(d, ignore_errors=True)
        removed.append(sid)
    store.conn.commit()
    store.add_event(f"purged {len(removed)} cancelled studies from task list",
                    level="info")
    return {"ok": True, "purged": removed, "count": len(removed)}


# ---------------------------------------------------------------------------
# job completion handling
# ---------------------------------------------------------------------------
def after_job(store: Store, cfg, job: dict, outcome: dict) -> dict:
    """Parse/validate a finished job, then advance, retry or fail the study."""
    parsed = parser.parse_output(outcome["output_path"])
    parsed["walltime_s"] = outcome.get("walltime")
    parsed["exit_code"] = outcome.get("exit_code")
    # optional Shermo thermochemistry for frequency jobs
    if (job.get("meta") or {}).get("expect_freq") and cfg.get("shermo.enabled", True):
        try:
            from .shermo import run_shermo
            parsed["shermo"] = run_shermo(
                cfg, outcome["output_path"],
                energy_hartree=parsed.get("final_energy_hartree"))
        except Exception as exc:  # pragma: no cover - external tool
            parsed["shermo"] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
    validation = validate(parsed, job, cfg)
    store.save_result(job["id"], job["study_id"], job["stage"], parsed, validation)
    store.update_job(job["id"], status="completed", exit_code=outcome.get("exit_code"),
                     timed_out=int(bool(outcome.get("timed_out"))),
                     walltime=outcome.get("walltime"),
                     output_path=outcome.get("output_path"),
                     finished_at=now())
    if validation["overall"] == "fail" or outcome.get("timed_out"):
        store.update_job(job["id"], status="failed")

    study = store.get_study(job["study_id"])
    job_status = ("failed" if (validation["overall"] == "fail"
                               or outcome.get("timed_out")) else "completed")
    _log_job_result(cfg, job, parsed, validation,
                    batch_id=(study or {}).get("batch_id"), status=job_status)
    if validation["overall"] in ("pass", "warn") and not outcome.get("timed_out"):
        nxt = next_stage(study, job["stage"])
        if nxt:
            new_job = create_job_for_stage(store, cfg, study, nxt,
                                           parent_job=job["id"])
            store.update_study(job["study_id"], status="running")
            store.add_event(f"{job['stage']} ok -> next stage {nxt}",
                            study_id=job["study_id"], job_id=job["id"])
            return {"action": "next_stage", "next_job": new_job["id"],
                    "validation": validation}
        finalize_study(store, cfg, job["study_id"])
        store.add_event("study completed", study_id=job["study_id"], job_id=job["id"])
        return {"action": "study_completed", "validation": validation}

    # failure path ------------------------------------------------------
    diag = diagnose(parsed, job, validation, cfg,
                    output_path=outcome.get("output_path"))
    _patch_job_meta(store, job["id"], {"diagnosis": diag})
    store.update_job(job["id"], error_class=diag["class"])
    store.add_event(
        f"{job['stage']} failed: {diag['class']} -- {diag['title']}",
        study_id=job["study_id"], job_id=job["id"], level="error")

    r = diag.get("retry") or {}
    attempt = int(job.get("attempt") or 1)
    if r.get("allowed") and attempt <= int(cfg.limits.get("max_retries", 2)):
        mem = None
        if r.get("strategy") == "reduce_mem":
            mem = _halve_mem(job.get("meta", {}).get("mem") or cfg.resources["mem"])
        new_job = create_job_for_stage(
            store, cfg, study, job["stage"],
            attempt=attempt + 1,
            retry_of=job["id"],
            route_patch=r.get("route_patch"),
            use_last_geometry=bool(r.get("use_last_geometry")),
            mem_override=mem,
            parent_job=job["id"],
        )
        store.add_event(
            f"auto-retry {job['stage']} attempt {attempt + 1} "
            f"(strategy={r['strategy']})",
            study_id=job["study_id"], job_id=new_job["id"])
        store.update_study(job["study_id"], status="running")
        return {"action": "retry", "next_job": new_job["id"],
                "diagnosis": diag, "validation": validation}

    _patch_meta(store, job["study_id"], {"last_diagnosis": diag})
    store.update_study(job["study_id"], status="failed")
    return {"action": "study_failed", "diagnosis": diag, "validation": validation}


def _halve_mem(mem: str) -> str:
    import re
    m = re.match(r"(\d+(?:\.\d+)?)\s*(GB|MB|GW|MW)?", str(mem), re.I)
    if not m:
        return "2GB"
    val = float(m.group(1)) / 2.0
    unit = (m.group(2) or "GB").upper()
    return f"{max(0.5, val):g}{unit}"


# ---------------------------------------------------------------------------
# real-time result log (crash safety)
# ---------------------------------------------------------------------------
def _log_job_result(cfg, job: dict, parsed: dict, validation: dict,
                    batch_id: str | None = None,
                    status: str | None = None) -> None:
    """Append one compact job result to the real-time JSONL log."""
    try:
        from .resultlog import append_record
        freq = parsed.get("freq") or {}
        orbitals = parsed.get("orbitals") or {}
        append_record(cfg, {
            "type": "job_result",
            "batch_id": batch_id,
            "study_id": job.get("study_id"),
            "job_id": job.get("id"),
            "stage": job.get("stage"),
            "route": job.get("route"),
            "status": status or job.get("status"),
            "validation": (validation or {}).get("overall"),
            "energy_hartree": parsed.get("final_energy_hartree"),
            "n_imag": freq.get("n_imag"),
            "homo_ev": orbitals.get("homo_ev"),
            "lumo_ev": orbitals.get("lumo_ev"),
            "opt_converged": (parsed.get("opt") or {}).get("converged"),
            "normal_termination": parsed.get("normal_termination"),
            "counterpoise": parsed.get("counterpoise"),
        })
    except Exception:
        pass


def _log_study_result(cfg, summary: dict) -> None:
    """Append the completed-study summary to the real-time JSONL log."""
    try:
        from .resultlog import append_record
        clean = {k: v for k, v in summary.items() if k != "final_geometry"}
        append_record(cfg, {
            "type": "study_result",
            "batch_id": summary.get("batch_id"),
            "study_id": summary.get("study_id"),
            "name": summary.get("name"),
            "smiles": summary.get("smiles"),
            "role": (summary.get("variant") or {}).get("role"),
            "summary": clean,
        })
    except Exception:
        pass


# ---------------------------------------------------------------------------
# AI-generated calculation chain (validated before registration)
# ---------------------------------------------------------------------------
AI_STAGE_KINDS = ("optfreq", "opt", "freq", "sp")


def _validate_ai_stage(cfg, stage: dict, index: int) -> tuple[dict, list[str]]:
    problems: list[str] = []
    kind = str(stage.get("kind") or "").strip().lower()
    name = str(stage.get("name") or f"stage{index}").strip()
    route = str(stage.get("route") or "").strip()
    phase = str(stage.get("phase") or "gas").strip().lower()
    if kind not in AI_STAGE_KINDS:
        problems.append(f"{name}: kind 只支持 {', '.join(AI_STAGE_KINDS)}（预优化由程序自动处理）")
    if phase not in ("gas", "solvent"):
        problems.append(f"{name}: phase 只支持 gas/solvent")
    if not route:
        problems.append(f"{name}: 缺少 route")
    else:
        route = normalize_route(route)
        try:
            tokens = split_route(route)
            base, _removed = strip_stage_keywords(tokens)
            route = " ".join(base)
            if any(t.lower().startswith(("scrf", "smd")) for t in base):
                problems.append(
                    f"{name}: 不要在 route 里写 scrf/smd；phase=solvent 时程序会自动加 PCM")
            check = validate_route(base, cfg)
            problems.extend(f"{name}: {p}" for p in check.get("problems") or [])
            problems.extend(f"{name}: {w}" for w in check.get("warnings") or [])
        except Exception as exc:
            problems.append(f"{name}: 路由解析失败 {type(exc).__name__}: {exc}")
    spec = {"name": name, "kind": kind, "route": route, "phase": phase}
    pop = stage.get("pop")
    if pop:
        spec["pop"] = str(pop)
    if phase == "solvent":
        spec["solvent"] = str(stage.get("solvent") or "").strip() or None
    return spec, problems


def register_ai_recipe(cfg, recipe: dict, key: str | None = None) -> dict:
    """Validate an LLM-generated calculation chain and store it as a recipe."""
    stages = recipe.get("core_stages") or recipe.get("stages") or []
    if not isinstance(stages, list) or not stages:
        return {"ok": False, "problems": ["core_stages 为空"], "key": None}
    if len(stages) > 8:
        return {"ok": False, "problems": ["计算链最多 8 个阶段"], "key": None}
    specs, problems = [], []
    for i, st in enumerate(stages, 1):
        if not isinstance(st, dict):
            problems.append(f"stage{i}: 不是对象")
            continue
        spec, p = _validate_ai_stage(cfg, st, i)
        specs.append(spec)
        problems.extend(p)
    names = [s["name"] for s in specs]
    if len(set(names)) != len(names):
        problems.append("阶段名重复：" + ", ".join(names))
    solvent = (recipe.get("solvent") or "").strip() or None
    if not solvent:
        solvent = next((s.get("solvent") for s in specs if s.get("solvent")), None)
    for spec in specs:
        if spec.get("phase") == "solvent" and not spec.get("solvent"):
            problems.append(f"{spec['name']}: 溶剂相 stage 必须给 solvent 名称")
    if problems:
        return {"ok": False, "problems": problems, "key": None}
    if not key:
        base = str(recipe.get("name") or "ai_workflow").strip()
        slug = "".join(c if (c.isascii() and c.isalnum()) else "_"
                       for c in base)[:40].strip("_")
        key = f"ai_{slug}" if slug else f"ai_{new_id('')}"
    existing = cfg.data.setdefault("recipes", {})
    if key in existing:
        key = f"{key}_{new_id()[:6]}"
    clean = {
        "name": str(recipe.get("name") or key)[:80],
        "description": str(recipe.get("description") or "AI 生成的计算链")[:300],
        "solvent": solvent,
        "core_stages": specs,
        "li": {
            "on": bool((recipe.get("li") or {}).get("on", False)),
            "elements": list((recipe.get("li") or {}).get("elements")
                             or ["O", "F"]),
            "bsse_levels": list((recipe.get("li") or {}).get("bsse_levels")
                                 or ["optfreq", "sp_gas"]),
        },
        "redox": {
            "on": bool((recipe.get("redox") or {}).get("on", False)),
            "states": list((recipe.get("redox") or {}).get("states")
                           or ["cation", "anion"]),
            "targets": list((recipe.get("redox") or {}).get("targets")
                            or ["neutral"]),
        },
        "ai_generated": True,
    }
    existing[key] = clean
    cfg.save()
    return {"ok": True, "problems": [], "key": key, "recipe": clean}


# ---------------------------------------------------------------------------
# summaries / finalisation
# ---------------------------------------------------------------------------
def finalize_study(store: Store, cfg, study_id: str) -> dict:
    summary = study_summary(store, cfg, study_id)
    summary["status"] = "completed"
    study_dir = _study_dir(cfg, study_id)
    full = dict(summary)
    (study_dir / "result.json").write_text(
        json.dumps(full, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8")
    meta = {k: v for k, v in summary.items() if k != "final_geometry"}
    _patch_meta(store, study_id, {"result": meta,
                                 "completed_at": now()})
    store.update_study(study_id, status="completed")
    _log_study_result(cfg, summary)
    return summary


def job_summary_rows(store: Store, study_id: str) -> list[dict]:
    rows = []
    for job in store.study_jobs(study_id):
        res = store.get_result(job["id"]) or {}
        parsed = res.get("parsed") or {}
        validation = res.get("validation") or {}
        diag = (job.get("meta") or {}).get("diagnosis") or {}
        rows.append({
            "job_id": job["id"],
            "stage": job["stage"],
            "attempt": job.get("attempt"),
            "status": job["status"],
            "route": job.get("route"),
            "walltime_s": job.get("walltime"),
            "exit_code": job.get("exit_code"),
            "timed_out": bool(job.get("timed_out")),
            "error_class": job.get("error_class"),
            "validation": validation.get("overall"),
            "normal_termination": parsed.get("normal_termination"),
            "energy_hartree": parsed.get("final_energy_hartree"),
            "n_imag": (parsed.get("freq") or {}).get("n_imag"),
            "opt_converged": (parsed.get("opt") or {}).get("converged"),
            "output": job.get("output_path"),
        })
    return rows


def _final_job(store: Store, study: dict) -> dict | None:
    jobs = store.study_jobs(study["id"])
    seq = stage_sequence(study)
    if not seq:
        return None
    last_stage = seq[-1]
    # prefer the completed last-stage job that passed validation
    for job in reversed([j for j in jobs if j["stage"] == last_stage]):
        if job["status"] == "completed":
            res = store.get_result(job["id"]) or {}
            if (res.get("validation") or {}).get("overall") in ("pass", "warn"):
                return job
    # otherwise: the last job that produced any result (honest failure summary)
    for job in reversed(jobs):
        if store.get_result(job["id"]):
            return job
    return None


def _stage_results(store: Store, study_id: str) -> dict:
    """Per-stage structured facts (energy/thermo/orbitals/BSSE) of a study."""
    out: dict = {}
    for job in store.study_jobs(study_id):
        res = store.get_result(job["id"]) or {}
        parsed = res.get("parsed") or {}
        validation = res.get("validation") or {}
        spec = (job.get("meta") or {}).get("stage_spec") or {}
        freq = parsed.get("freq") or {}
        orbitals = parsed.get("orbitals") or {}
        out[job["stage"]] = {
            "job_id": job["id"], "status": job["status"],
            "kind": spec.get("kind"), "phase": spec.get("phase"),
            "route": job.get("route"),
            "energy_hartree": parsed.get("final_energy_hartree"),
            "energy_ev": parsed.get("final_energy_ev"),
            "zpe_hartree": freq.get("zpve_hartree"),
            "enthalpy_hartree": freq.get("enthalpy_hartree"),
            "freeenergy_hartree": freq.get("freeenergy_hartree"),
            "entropy_cal_mol_k": freq.get("entropy_cal_mol_k"),
            "homo_ev": orbitals.get("homo_ev"),
            "lumo_ev": orbitals.get("lumo_ev"),
            "gap_ev": orbitals.get("gap_ev"),
            "counterpoise": parsed.get("counterpoise"),
            "n_imag": freq.get("n_imag"),
            "opt_converged": (parsed.get("opt") or {}).get("converged"),
            "validation": validation.get("overall"),
            "walltime_s": job.get("walltime"),
        }
    return out


def study_summary(store: Store, cfg, study_id: str) -> dict:
    study = store.get_study(study_id)
    if study is None:
        return {"ok": False, "error": f"未知 study: {study_id}"}
    meta = study.get("meta") or {}
    jobs = store.study_jobs(study_id)
    if not jobs and isinstance(meta.get("result"), dict):
        # recovered from the real-time result log: show the stored summary
        out = dict(meta["result"])
        out["ok"] = True
        out["recovered"] = True
        out.setdefault("study_id", study_id)
        out.setdefault("name", study["name"])
        out.setdefault("batch_id", study.get("batch_id"))
        return out
    final_job = _final_job(store, study)
    parsed = _result_of(store, final_job)
    validation = {}
    if final_job:
        validation = (store.get_result(final_job["id"]) or {}).get("validation") or {}
    freq = parsed.get("freq") or {}
    opt = parsed.get("opt") or {}
    wall = sum(float(j.get("walltime") or 0) for j in jobs)
    variant = meta.get("variant") or {}
    stage_results = _stage_results(store, study_id)
    orbitals = next((v for v in stage_results.values()
                     if v.get("homo_ev") is not None), {})
    counterpoise = {k: v["counterpoise"] for k, v in stage_results.items()
                    if v.get("counterpoise")}
    return {
        "ok": True,
        "study_id": study_id,
        "name": study["name"],
        "batch_id": study.get("batch_id"),
        "status": study["status"],
        "variant": variant,
        "smiles": study.get("smiles"),
        "formula": meta.get("formula"),
        "natoms": meta.get("natoms"),
        "charge": study.get("charge"),
        "multiplicity": study.get("multiplicity"),
        "route": study.get("route"),
        "workflow": study.get("workflow"),
        "solvent": study.get("solvent"),
        "preopt": bool(study.get("preopt")),
        "final_stage": final_job["stage"] if final_job else None,
        "final_job_id": final_job["id"] if final_job else None,
        "energy_hartree": parsed.get("final_energy_hartree"),
        "energy_kcal_mol": parsed.get("final_energy_kcal_mol"),
        "energy_ev": parsed.get("final_energy_ev"),
        "zpe_hartree": freq.get("zpve_hartree"),
        "enthalpy_hartree": freq.get("enthalpy_hartree"),
        "freeenergy_hartree": freq.get("freeenergy_hartree"),
        "entropy_cal_mol_k": freq.get("entropy_cal_mol_k"),
        "n_imag": freq.get("n_imag"),
        "imag_freqs": freq.get("imag_freqs"),
        "lowest_modes": freq.get("lowest_modes"),
        "dipole_debye": parsed.get("dipole_debye"),
        "opt_converged": opt.get("converged"),
        "opt_steps": opt.get("nsteps"),
        "scf_converged": parsed.get("scf_converged"),
        "scf_cycles": parsed.get("scf_cycles"),
        "normal_termination": parsed.get("normal_termination"),
        "validation_overall": validation.get("overall"),
        "validation_reasons": validation.get("reasons"),
        "warnings": parsed.get("warnings"),
        "walltime_total_s": round(wall, 1),
        "n_jobs": len(jobs),
        "created_at": study.get("created_at"),
        "updated_at": study.get("updated_at"),
        "needs_approval": meta.get("needs_approval"),
        "approval_reasons": meta.get("approval_reasons"),
        "last_diagnosis": meta.get("last_diagnosis"),
        "final_geometry": parsed.get("final_geometry"),
        "jobs": job_summary_rows(store, study_id),
        "stage_results": stage_results,
        "orbitals": orbitals,
        "counterpoise": counterpoise,
        "role": meta.get("role"),
        "recipe": meta.get("recipe"),
        "li_site": meta.get("li_site"),
        "redox": meta.get("redox"),
        "g16_version": parsed.get("version"),
    }


# ---------------------------------------------------------------------------
# batch / iteration
# ---------------------------------------------------------------------------
def plan_variants(n: int, ions: list[str], seed: int = 42,
                  include_base: bool = True) -> list[dict]:
    rng = random.Random(seed)
    variants: list[dict] = []
    if include_base:
        variants.append({"type": "base", "label": "base"})
    if not ions:
        if not include_base:
            return []
        while len(variants) < n:
            variants.append({"type": "base", "label": "base"})
        return variants[:n]
    for ion in ions:
        if len(variants) >= n:
            break
        variants.append({"type": "ion", "ion": ion, "label": ion})
    while len(variants) < n:
        ion = rng.choice(ions)
        variants.append({"type": "ion", "ion": ion, "label": ion})
    return variants[:n]


def create_batch(store: Store, cfg, *, base_name: str,
                 smiles: str | None = None, file: str | None = None,
                 n: int = 3, ions: list[str] | None = None, seed: int = 42,
                 route: str | None = None, workflow: str | None = None,
                 solvent: str | None = None, preopt: bool | None = None,
                 sp_route: str | None = None,
                 include_base: bool = True) -> dict:
    ions = ions or []
    n = int(n)
    max_n = int(cfg.limits.get("max_studies_per_batch", 60))
    if n > max_n:
        raise ValueError(f"变体数 {n} 超过上限 {max_n}")
    batch_id = new_id("bt_")
    variants = plan_variants(n, ions, seed=seed, include_base=include_base)
    studies = []
    for i, v in enumerate(variants):
        label = v.get("label") or v.get("ion") or "base"
        study = create_study(
            store, cfg,
            name=f"{base_name}-{label}",
            smiles=smiles, file=file, route=route, workflow=workflow,
            solvent=solvent, preopt=preopt, sp_route=sp_route,
            ion=v.get("ion"), seed=seed + i, batch_id=batch_id,
            variant=v,
        )
        studies.append(study)
    store.add_event(f"batch created: {batch_id} ({len(studies)} variants)",
                    study_id=None)
    return {"ok": True, "batch_id": batch_id,
            "studies": [{"study_id": s["id"], "name": s["name"],
                         "variant": (s.get("meta") or {}).get("variant"),
                         "status": s["status"]} for s in studies],
            "count": len(studies)}


def submit_batch(store: Store, cfg, batch_id: str, approve: bool = False) -> dict:
    studies = store.list_studies(batch_id=batch_id)
    if not studies:
        return {"ok": False, "error": f"未知 batch: {batch_id}"}
    results = []
    for s in studies:
        results.append(submit_study(store, cfg, s["id"], approve=approve))
    return {"ok": True, "batch_id": batch_id, "submitted": results}
