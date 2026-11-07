"""Preset calculation chains ("recipes") for the Gaussian workflow.

A recipe is a system built from several Gaussian studies:

    neutral molecule            : opt+freq → gas SP → solvent SP
    Li⁺@O / Li⁺@F complexes     : opt+freq → gas SP → solvent SP → BSSE SP(×levels)
    Li⁺ reference (single ion)  : gas/solvent SP only
    cation / anion (redox)      : opt+freq → gas SP → solvent SP
                                  (optionally also for the Li⁺ complexes)

``plan_studies`` produces the study skeleton; :mod:`gaussflow.pipeline`
materialises it in the database and :mod:`gaussflow.analysis` later combines
the energies into binding energies / redox potentials.
"""
from __future__ import annotations

import copy


# ---------------------------------------------------------------------------
def list_recipes(cfg) -> dict[str, str]:
    return {k: v.get("name", k) for k, v in (cfg.get("recipes") or {}).items()}


def get_recipe(cfg, name: str) -> dict:
    recipes = cfg.get("recipes") or {}
    if name not in recipes:
        raise ValueError(f"未知计算方案: {name}（可选: {', '.join(recipes)}）")
    return copy.deepcopy(recipes[name])


def mult_for_charge(total_z: int, charge: int) -> int:
    """Closed-shell multiplicity guess from electron parity."""
    electrons = int(total_z) - int(charge)
    return 1 if electrons % 2 == 0 else 2


# ---------------------------------------------------------------------------
def _level_stage(core_stages: list[dict], level: str) -> dict | None:
    for stage in core_stages:
        if stage["name"] == level:
            return stage
    return None


def _bsse_stages(core_stages: list[dict], levels: list[str]) -> tuple[list[dict], list[str]]:
    """BSSE single points at the requested levels.

    Gaussian 16 refuses ``scrf`` together with ``counterpoise``
    ("No SCRF with counterpoise/fragments."), so solvent-phase BSSE is skipped
    with an explicit note instead of silently producing a failed job.
    """
    out: list[dict] = []
    notes: list[str] = []
    for level in levels:
        src = _level_stage(core_stages, level)
        if src is None:
            continue
        if src.get("phase") == "solvent":
            notes.append(
                f"Gaussian 16 不支持 PCM 与 counterpoise 同时使用，已跳过 "
                f"{level} 层级的 BSSE（该层级只给未校正结合能）。")
            continue
        view = copy.deepcopy(src)
        view["name"] = f"bsse_{level}"
        view["kind"] = "bsse"
        view["level"] = level
        view.pop("pop", None)
        out.append(view)
    return out, notes


def plan_studies(cfg, recipe_name: str, *,
                 name: str,
                 charge: int,
                 multiplicity: int,
                 total_z: int,
                 natoms: int,
                 sites: list[dict],
                 options: dict | None = None) -> dict:
    """Return the study plan for one molecular system (no DB / no Gaussian).

    Parameters
    ----------
    sites : list from :func:`gaussflow.structure.list_coordination_sites`
    options : overrides
        ``li`` (bool), ``li_elements`` (list), ``redox`` (bool),
        ``redox_states`` (list), ``redox_targets`` (list), ``preopt`` (bool),
        ``solvent`` (str or None), ``bsse_levels`` (list)
    """
    recipe = get_recipe(cfg, recipe_name)
    opts = {
        "li": bool((recipe.get("li") or {}).get("on", False)),
        "li_elements": list((recipe.get("li") or {}).get("elements", ["O", "F"])),
        "redox": bool((recipe.get("redox") or {}).get("on", False)),
        "redox_states": list((recipe.get("redox") or {}).get("states", ["cation", "anion"])),
        "redox_targets": list((recipe.get("redox") or {}).get("targets", ["neutral"])),
        "preopt": True,
        "solvent": recipe.get("solvent"),
        "bsse_levels": list((recipe.get("li") or {}).get(
            "bsse_levels", ["optfreq", "sp_gas", "sp_sol"])),
        "max_sites_per_element": int((recipe.get("li") or {}).get(
            "max_sites_per_element", 4)),
        "li_dist": dict((recipe.get("li") or {}).get("dist", {"O": 2.05, "F": 1.90})),
    }
    opts.update({k: v for k, v in (options or {}).items() if v is not None})

    core = copy.deepcopy(recipe["core_stages"])
    if opts["solvent"]:
        for st in core:
            if st.get("phase") == "solvent":
                st["solvent"] = opts["solvent"]
    if opts["preopt"]:
        preopt_route = cfg.chemistry.get("preopt_route") or "HF/3-21G"
        core.insert(0, {"name": "preopt", "kind": "opt", "route": preopt_route,
                        "phase": "gas"})

    studies: list[dict] = []

    def add(name_, role, charge_, mult_, stages, variant, extra=None):
        item = {
            "name": name_, "role": role, "charge": int(charge_),
            "multiplicity": int(mult_), "stage_specs": copy.deepcopy(stages),
            "variant": variant, "natoms": natoms,
        }
        item.update(extra or {})
        studies.append(item)

    # 1) neutral molecule ------------------------------------------------
    add(f"{name}-neutral", "neutral", charge, multiplicity, core,
        {"type": "neutral", "label": "neutral"})

    # 2) Li+ complexes ----------------------------------------------------
    li_sites = []
    li_z = total_z + 3          # Li has Z=3
    notes: list[str] = []
    if opts["li"]:
        wanted = [e.upper() for e in opts["li_elements"]]
        per_element: dict[str, int] = {}
        for site in sites:
            el = site["element"].upper()
            if el not in wanted:
                continue
            if per_element.get(el, 0) >= opts["max_sites_per_element"]:
                continue
            per_element[el] = per_element.get(el, 0) + 1
            li_sites.append(site)
        bsse_stages, bsse_notes = _bsse_stages(core, opts["bsse_levels"])
        notes.extend(bsse_notes)
        for site in li_sites:
            el = site["element"]
            site_no = site["idx_in_element"] + 1
            li_charge = charge + 1
            li_mult = mult_for_charge(li_z, li_charge)
            stages = copy.deepcopy(core) + copy.deepcopy(bsse_stages)
            fragments = [
                {"n": natoms, "q": charge, "m": multiplicity},
                {"n": 1, "q": 1, "m": 1},
            ]
            for st in stages:
                if st["kind"] == "bsse":
                    st["fragments"] = fragments
            add(f"{name}-Li{el}-{site_no}", f"li_{el}{site_no}", li_charge,
                li_mult, stages,
                {"type": "li_complex", "label": f"Li⁺@{el}{site_no}", "element": el},
                {"li_site": {"element": el, "atom_index": site["atom_index"],
                             "site_no": site_no,
                             "dist": opts["li_dist"].get(el, 2.0)},
                 "fragments": fragments})

    # 3) redox species ----------------------------------------------------
    if opts["redox"]:
        redox_core = copy.deepcopy(core)
        targets = []
        if "neutral" in opts["redox_targets"]:
            targets.append(("neutral", charge, multiplicity,
                            {"type": "neutral", "label": "neutral"}, None))
        if "li" in opts["redox_targets"]:
            for site in li_sites:
                el = site["element"]
                site_no = site["idx_in_element"] + 1
                targets.append((f"Li{el}{site_no}", charge + 1,
                                mult_for_charge(li_z, charge + 1),
                                {"type": "li_complex", "label": f"Li⁺@{el}{site_no}",
                                 "element": el}, site))
        state_delta = {"cation": +1, "anion": -1}
        for state in opts["redox_states"]:
            dq = state_delta.get(state)
            if dq is None:
                continue
            for tag, base_q, base_m, variant, site in targets:
                q = base_q + dq
                z = total_z if site is None else li_z
                m = mult_for_charge(z, q)
                role = f"redox_{state}" if tag == "neutral" else f"{tag}_redox_{state}"
                extra = {"redox": {"state": state, "parent": tag}}
                if site:
                    extra["li_site"] = {"element": site["element"],
                                        "atom_index": site["atom_index"],
                                        "site_no": site["idx_in_element"] + 1,
                                        "dist": opts["li_dist"].get(site["element"], 2.0)}
                add(f"{name}-{tag}-{state}", role, q, m,
                    copy.deepcopy(redox_core),
                    {**variant, "state": state,
                     "label": f"{variant['label']} ({state})"},
                    extra)

    # 4) Li+ single-ion reference ----------------------------------------
    if opts["li"]:
        ref_stages = []
        level_names = {"optfreq": "opt", "sp_gas": "gas", "sp_sol": "sol"}
        for level in ("optfreq", "sp_gas", "sp_sol"):
            src = _level_stage(core, level)
            if src is None:
                continue
            ref_stages.append({
                "name": f"li_{level_names.get(level, level)}", "kind": "sp",
                "route": src["route"], "phase": src.get("phase", "gas"),
                "solvent": src.get("solvent"),
            })
        if ref_stages:
            add(f"{name}-Li+ref", "li_ref", 1, 1, ref_stages,
                {"type": "li_ref", "label": "Li⁺ ref"}, {"natoms": 1})

    return {
        "recipe": recipe_name,
        "recipe_name": recipe.get("name", recipe_name),
        "solvent": opts["solvent"],
        "stage_names": [s["name"] for s in core],
        "li_enabled": opts["li"],
        "li_sites": li_sites,
        "redox_enabled": opts["redox"],
        "studies": studies,
        "n_studies": len(studies),
        "notes": notes,
    }


def route_of_stage(study_meta: dict, stage_name: str) -> str | None:
    for stage in (study_meta.get("stage_specs") or []):
        if stage["name"] == stage_name:
            return stage.get("route")
    return None
