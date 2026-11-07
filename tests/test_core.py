"""Unit tests for the deterministic core (no Gaussian needed except fixtures).

Run from the project root:  python -m unittest discover -s tests -v
"""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gaussflow.config import Config                      # noqa: E402
from gaussflow.diagnose import diagnose                  # noqa: E402
from gaussflow.gjf import (build_counterpoise_gjf, build_gjf,  # noqa: E402
                           compose_route, split_route,
                           strip_stage_keywords, validate_route)
from gaussflow.parse import parse_output                 # noqa: E402
from gaussflow.store import Store                        # noqa: E402
from gaussflow.structure import (add_ion, coords_from_mol,  # noqa: E402
                                 formula_of_coords, infer_charge_mult,
                                 molecule_from_smiles)
from gaussflow.validate import validate                  # noqa: E402

DATA = Path(__file__).resolve().parent / "data"


def _tmp_config() -> Config:
    tmp = Path(tempfile.mkdtemp(prefix="gflow_test_"))
    return Config.load(tmp / "config.json")


class TestRoute(unittest.TestCase):
    def test_split_and_strip(self):
        toks = split_route("#p B3LYP/6-31G(d) opt freq")
        self.assertEqual(toks[0], "B3LYP/6-31G(d)")
        base, removed = strip_stage_keywords(toks)
        self.assertEqual(base, ["B3LYP/6-31G(d)"])
        self.assertIn("opt", removed)

    def test_compose(self):
        route = compose_route(["B3LYP/6-31G(d)"], "freq", "water",
                              ["guess=read"])
        self.assertIn("freq", route)
        self.assertIn("scrf=(smd,solvent=water)", route)
        self.assertIn("guess=read", route)

    def test_validation_ok(self):
        cfg = _tmp_config()
        self.assertTrue(validate_route(["B3LYP/6-31G(d)"], cfg)["ok"])
        self.assertTrue(validate_route(["HF/3-21G"], cfg)["ok"])

    def test_validation_rejects_unknown_method(self):
        cfg = _tmp_config()
        res = validate_route(["FOO/6-31G(d)"], cfg)
        self.assertFalse(res["ok"])
        self.assertTrue(res["problems"])

    def test_validation_unknown_basis_needs_approval(self):
        cfg = _tmp_config()
        res = validate_route(["B3LYP/not-a-basis"], cfg)
        self.assertTrue(res["ok"])
        self.assertTrue(res["needs_approval"])

    def test_retry_patch_replaces_stage(self):
        route = compose_route(["B3LYP/6-31G(d)"], "opt", None, None,
                              "opt=(maxcycles=200,calcfc)")
        self.assertIn("opt=(maxcycles=200,calcfc)", route)
        self.assertNotIn(" opt opt=", route)
        # scf patch is appended, not replacing the stage keyword
        route2 = compose_route(["B3LYP/6-31G(d)"], "opt", None, None,
                               "scf=(maxcycle=512,xqc)")
        self.assertIn(" opt ", route2 + " ")
        self.assertIn("scf=(maxcycle=512,xqc)", route2)


class TestGjf(unittest.TestCase):
    def test_build(self):
        text = build_gjf(route_line="#p HF/3-21G opt", title="water",
                         charge=0, multiplicity=1,
                         coords=[("O", 0.0, 0.0, 0.1173),
                                 ("H", 0.0, 0.7572, -0.4692),
                                 ("H", 0.0, -0.7572, -0.4692)],
                         mem="4GB", nproc=4)
        self.assertIn("%chk=job.chk", text)
        self.assertIn("%mem=4GB", text)
        self.assertIn("%nproc=4", text)
        self.assertIn("#p HF/3-21G opt", text)
        self.assertIn("0 1", text)
        self.assertIn("O ", text)
        self.assertTrue(text.endswith("\n\n"))

    def test_oldchk(self):
        text = build_gjf(route_line="#p HF/3-21G freq guess=read", title="t",
                         charge=0, multiplicity=1,
                         coords=[("O", 0.0, 0.0, 0.0)], oldchk="prev.chk")
        self.assertIn("%oldchk=prev.chk", text)


class TestStructure(unittest.TestCase):
    def test_ethanol_3d(self):
        mol = molecule_from_smiles("CCO", seed=42)
        coords = coords_from_mol(mol)
        self.assertEqual(len(coords), 9)
        self.assertEqual(formula_of_coords(coords), "C2H6O")
        self.assertEqual(infer_charge_mult(mol), (0, 1))

    def test_ions(self):
        mol = molecule_from_smiles("CCO", seed=42)
        na = add_ion(mol, "Na+")
        self.assertEqual(len(coords_from_mol(na)), 10)
        self.assertEqual(infer_charge_mult(na), (1, 1))
        cl = add_ion(mol, "Cl-")
        self.assertEqual(infer_charge_mult(cl), (-1, 1))

    def test_water(self):
        mol = molecule_from_smiles("O", seed=1)
        self.assertEqual(len(coords_from_mol(mol)), 3)


class TestParse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.water = parse_output(DATA / "water_hf_sto3g_sp.out")
        cls.ethanol = parse_output(DATA / "ethanol_hf_321g_opt.out")
        cls.pcm = parse_output(DATA / "water_pcm.out")
        cls.cp = parse_output(DATA / "cp_li_water.out")

    def test_water(self):
        p = self.water
        self.assertTrue(p["normal_termination"])
        self.assertEqual(p["natoms"], 3)
        self.assertAlmostEqual(p["final_energy_hartree"], -74.9630231626, places=6)
        self.assertTrue(p["scf_converged"])
        self.assertEqual(len(p["final_geometry"]), 3)
        self.assertEqual(p["final_geometry"][0]["symbol"], "O")
        self.assertIsNotNone(p["final_energy_kcal_mol"])

    def test_ethanol_opt(self):
        p = self.ethanol
        self.assertTrue(p["normal_termination"])
        self.assertEqual(p["natoms"], 9)
        self.assertTrue(p["is_opt"])
        self.assertTrue(p["opt"]["converged"])
        self.assertGreaterEqual(p["opt"]["nsteps"], 1)
        self.assertEqual(len(p["final_geometry"]), 9)

    def test_pcm_orbitals(self):
        p = self.pcm
        self.assertTrue(p["normal_termination"])
        orb = p["orbitals"]
        self.assertIsNotNone(orb["homo_ev"])
        self.assertIsNotNone(orb["lumo_ev"])
        self.assertGreater(orb["lumo_ev"], orb["homo_ev"])
        self.assertAlmostEqual(orb["gap_ev"], orb["lumo_ev"] - orb["homo_ev"], places=9)
        # water: 5 occupied MOs in 6-31+G(d)
        self.assertEqual(orb["n_occupied_basis"], 5)
        self.assertIn("Polarizable Continuum Model", (DATA / "water_pcm.out").read_text(
            encoding="utf-8", errors="replace"))

    def test_counterpoise(self):
        cp = self.cp["counterpoise"]
        self.assertIsNotNone(cp)
        self.assertIsNotNone(cp["bsse_energy_hartree"])
        self.assertIsNotNone(cp["corrected_energy_hartree"])
        self.assertGreater(cp["bsse_energy_hartree"], 0)
        # Gaussian convention: corrected = complex + BSSE (corrected > complex)
        self.assertAlmostEqual(
            cp["corrected_energy_hartree"],
            cp["complex_energy_hartree"] + cp["bsse_energy_hartree"], places=6)
        # raw + BSSE = corrected (in kcal/mol as printed)
        if cp["complexation_raw_kcal"] is not None:
            raw = cp["complexation_raw_kcal"]
            corr = cp["complexation_corrected_kcal"]
            bsse_kcal = cp["bsse_energy_hartree"] * 627.5094740631
            self.assertAlmostEqual(raw + bsse_kcal, corr, delta=0.05)


class TestCounterpoiseGjf(unittest.TestCase):
    def test_total_first_and_fragment_labels(self):
        text = build_counterpoise_gjf(
            route_line="#p HF/3-21G counterpoise=2", title="t",
            coords=[("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.0),
                    ("Li", 0.0, 0.0, -1.85)],
            fragments=[{"n": 2, "q": 0, "m": 1}, {"n": 1, "q": 1, "m": 1}],
            total_charge=1, total_multiplicity=1)
        lines = text.splitlines()
        self.assertIn("1 1 0 1 1 1", lines)
        self.assertIn("O(Fragment=1)", text)
        self.assertIn("Li(Fragment=2)", text)


class TestRecipes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = _tmp_config()

    def test_plan_water_with_li_and_redox(self):
        from gaussflow.recipes import plan_studies
        from gaussflow.structure import (list_coordination_sites,
                                         molecule_from_smiles)
        mol = molecule_from_smiles("O", seed=42)
        sites = list_coordination_sites(mol, ["O", "F"])
        plan = plan_studies(
            self.cfg, "standard", name="w", charge=0, multiplicity=1,
            total_z=10, natoms=3, sites=sites,
            options={"li": True, "li_elements": ["O"], "redox": True,
                     "redox_states": ["cation", "anion"],
                     "redox_targets": ["neutral"], "preopt": True,
                     "solvent": "DMSO"})
        roles = [s["role"] for s in plan["studies"]]
        self.assertIn("neutral", roles)
        self.assertIn("li_O1", roles)
        self.assertIn("redox_cation", roles)
        self.assertIn("redox_anion", roles)
        self.assertIn("li_ref", roles)
        neutral = next(s for s in plan["studies"] if s["role"] == "neutral")
        self.assertEqual([x["name"] for x in neutral["stage_specs"]],
                         ["preopt", "optfreq", "sp_gas", "sp_sol"])
        li = next(s for s in plan["studies"] if s["role"] == "li_O1")
        # water + Li+: 13 electrons, charge +1 -> singlet
        self.assertEqual((li["charge"], li["multiplicity"]), (1, 1))
        self.assertTrue(any(s["kind"] == "bsse" for s in li["stage_specs"]))
        cat = next(s for s in plan["studies"] if s["role"] == "redox_cation")
        self.assertEqual((cat["charge"], cat["multiplicity"]), (1, 2))

    def test_route_composition(self):
        from gaussflow.pipeline import _compose_recipe_route
        spec = {"name": "optfreq", "kind": "optfreq",
                "route": "B3LYP/6-311+G(d) empiricaldispersion=gd3bj",
                "phase": "gas"}
        route = _compose_recipe_route(spec)
        self.assertIn("B3LYP/6-311+G(d)", route)
        self.assertIn("opt freq", route)
        # retry patch replaces, not duplicates, the stage keyword
        patched = _compose_recipe_route(spec, "opt=(maxcycles=200,calcfc)")
        self.assertIn("opt=(maxcycles=200,calcfc)", patched)
        self.assertNotIn(" opt freq opt=(maxcycles", patched)
        spec_sol = {"name": "sp_sol", "kind": "sp", "route": "M062X/6-31+G(d)",
                    "phase": "solvent", "solvent": "DMSO", "pop": "full"}
        route_sol = _compose_recipe_route(spec_sol)
        self.assertIn("scrf=(pcm,solvent=DMSO)", route_sol)
        self.assertIn("pop=full", route_sol)

    def test_default_recipes_are_valid_routes(self):
        from gaussflow.pipeline import _compose_recipe_route
        from gaussflow.gjf import split_route, validate_route
        recipes = self.cfg.get("recipes") or {}
        for key, recipe in recipes.items():
            for stage in recipe.get("core_stages", []):
                check = validate_route(split_route(stage["route"]), self.cfg)
                self.assertTrue(check["ok"], f"{key}/{stage['name']}: {check['problems']}")
                _compose_recipe_route(stage)


class TestShermo(unittest.TestCase):
    def test_parse_block(self):
        from gaussflow.shermo import parse_shermo_output
        text = (
            " ========== Total ==========\n"
            " Total S:      224.240 J/mol/K      53.595 cal/mol/K    -TS:    18.758 kcal/mol\n"
            " Sum of electronic energy and ZPE, namely U/H/G at 0 K:      -114.5231595 a.u.\n"
            " Sum of electronic energy and thermal correction to U:       -114.5197478 a.u.\n"
            " Sum of electronic energy and thermal correction to H:       -114.5186394 a.u.\n"
            " Sum of electronic energy and thermal correction to G:       -114.5485323 a.u.\n"
        )
        out = parse_shermo_output(text)
        self.assertTrue(out["ok"])
        self.assertAlmostEqual(out["G_hartree"], -114.5485323)
        self.assertAlmostEqual(out["H_hartree"], -114.5186394)
        self.assertAlmostEqual(out["S_cal_mol_k"], 53.595)

    def test_not_installed_is_graceful(self):
        cfg = _tmp_config()
        from gaussflow.shermo import run_shermo
        res = run_shermo(cfg, DATA / "water_pcm.out")
        self.assertFalse(res.get("ok"))
        self.assertIn("error", res)


class TestLocate(unittest.TestCase):
    def test_detect_all_structure(self):
        from gaussflow.locate import detect_all
        found = detect_all()
        for key in ("gaussian", "gaussview", "shermo"):
            self.assertIn(key, found)
            self.assertIn("best", found[key])
            self.assertIsInstance(found[key]["candidates"], list)


class TestValidate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = _tmp_config()
        cls.water = parse_output(DATA / "water_hf_sto3g_sp.out")

    def test_pass(self):
        job = {"stage": "sp",
               "meta": {"expected_natoms": 3, "expected_charge": 0,
                        "expected_multiplicity": 1}}
        rep = validate(self.water, job, self.cfg)
        # HF/STO-3G carries a "small basis set" warning -> overall may be warn
        self.assertNotEqual(rep["overall"], "fail", rep["reasons"])
        self.assertTrue(all(c["status"] != "fail" for c in rep["checks"]))

    def test_atom_mismatch(self):
        job = {"stage": "sp",
               "meta": {"expected_natoms": 4, "expected_charge": 0,
                        "expected_multiplicity": 1}}
        rep = validate(self.water, job, self.cfg)
        self.assertEqual(rep["overall"], "fail")
        self.assertTrue(any(c["id"] == "atom_count" for c in rep["checks"]))

    def test_imaginary_freq_fails(self):
        parsed = {
            "exists": True, "normal_termination": True, "error_termination": False,
            "error_lines": [], "natoms": 3, "charge": 0, "multiplicity": 1,
            "scf_converged": True, "final_energy_hartree": -75.0,
            "opt": None, "is_opt": False,
            "freq": {"has_freq": True, "n_imag": 2, "imag_freqs": [-120.0, -80.0],
                     "zpve_hartree": 0.02, "enthalpy_hartree": -74.9,
                     "freeenergy_hartree": -74.95},
            "warnings": [],
        }
        job = {"stage": "freq", "meta": {"expected_natoms": 3}}
        rep = validate(parsed, job, self.cfg)
        self.assertEqual(rep["overall"], "fail")
        self.assertTrue(any(c["id"] == "no_imaginary" for c in rep["checks"]))


class TestDiagnose(unittest.TestCase):
    def test_opt_not_converged_retry(self):
        cfg = _tmp_config()
        parsed = {"exists": True, "normal_termination": False,
                  "error_termination": True,
                  "error_lines": ["Optimization stopped.",
                                  "Error termination via Lnk1e in D:\\g16w\\l9999.exe"],
                  "final_energy_hartree": -100.0, "final_geometry": []}
        job = {"id": "j1", "stage": "opt", "attempt": 1,
               "meta": {"mem": "4GB"}, "timed_out": 0}
        validation = {"overall": "fail", "checks": [], "reasons": []}
        diag = diagnose(parsed, job, validation, cfg)
        self.assertEqual(diag["class"], "opt_not_converged")
        self.assertTrue(diag["retry"]["allowed"])
        self.assertIsNotNone(diag["retry"]["route_patch"])

    def test_timeout_needs_human(self):
        cfg = _tmp_config()
        parsed = {"exists": False, "normal_termination": False,
                  "error_termination": True, "error_lines": []}
        job = {"id": "j2", "stage": "opt", "attempt": 1, "timed_out": 1,
               "meta": {}}
        diag = diagnose(parsed, job, {"overall": "fail", "checks": [],
                                      "reasons": []}, cfg)
        self.assertEqual(diag["class"], "timeout")
        self.assertTrue(diag["needs_human"])
        self.assertFalse(diag["retry"]["allowed"])

    def test_transient_retry(self):
        cfg = _tmp_config()
        parsed = {"exists": True, "normal_termination": False,
                  "error_termination": True, "error_lines": [],
                  "final_energy_hartree": None}
        job = {"id": "j3", "stage": "opt", "attempt": 1, "meta": {}}
        diag = diagnose(parsed, job, {"overall": "fail", "checks": [],
                                      "reasons": []}, cfg)
        self.assertEqual(diag["class"], "killed_or_truncated")
        self.assertTrue(diag["retry"]["allowed"])


class TestStore(unittest.TestCase):
    def test_roundtrip(self):
        tmp = Path(tempfile.mkdtemp(prefix="gflow_store_"))
        try:
            store = Store(tmp / "t.db")
            study = store.create_study(name="s", smiles="O", charge=0,
                                       multiplicity=1, route="HF/3-21G",
                                       workflow="opt_freq", status="prepared",
                                       meta={"initial_coords": [], "x": 1})
            self.assertEqual(study["meta"]["x"], 1)
            job = store.create_job(study_id=study["id"], stage="opt", seq=1,
                                   route="#p HF/3-21G opt", charge=0,
                                   multiplicity=1, meta={"expected_natoms": 3})
            claimed = store.claim_next_job()
            self.assertEqual(claimed["id"], job["id"])
            self.assertEqual(store.claim_next_job(), None)
            store.save_result(job["id"], study["id"], "opt",
                              {"final_energy_hartree": -75.0},
                              {"overall": "pass"})
            res = store.get_result(job["id"])
            self.assertEqual(res["parsed"]["final_energy_hartree"], -75.0)
            store.close()
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
