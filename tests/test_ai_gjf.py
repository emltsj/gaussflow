"""Tests for DeepSeek AI generation helpers and multi-.gjf batch importing."""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gaussflow.config import Config                        # noqa: E402
from gaussflow import fluorinate, llm, pipeline, resultlog  # noqa: E402
from gaussflow import analysis                              # noqa: E402
from gaussflow.llm import (LlmError, _extract_json, chat,   # noqa: E402
                           parse_f_range, validate_molecules)
from gaussflow.pipeline import (create_batch_from_gjfs,     # noqa: E402
                                gjf_study_spec)
from gaussflow.store import Store                          # noqa: E402


def _tmp_config() -> Config:
    tmp = Path(tempfile.mkdtemp(prefix="gflow_ai_test_"))
    return Config.load(tmp / "config.json")


WATER_GJF = """/%chk=job.chk
%mem=1GB
%nproc=2
#p B3LYP/6-31G(d) opt freq

water test

0 1
O   0.000000   0.000000   0.117300
H   0.000000   0.757200  -0.469200
H   0.000000  -0.757200  -0.469200
"""


class TestFfRange(unittest.TestCase):
    def test_ranges(self):
        self.assertEqual(parse_f_range("0-3"), [0, 1, 2, 3])
        self.assertEqual(parse_f_range("3"), [3])
        self.assertEqual(parse_f_range("0,2,4"), [0, 2, 4])
        self.assertEqual(parse_f_range("4-2"), [2, 3, 4])
        self.assertEqual(parse_f_range(""), [])
        self.assertEqual(parse_f_range(None), [])


class TestValidateMolecules(unittest.TestCase):
    def test_valid_and_invalid(self):
        payload = {"molecules": [
            {"name": "TEP-F0", "smiles": "CCOP(=O)(OCC)OCC", "charge": 0,
             "multiplicity": 1},
            {"name": "bad", "smiles": "not-a-smiles", "charge": 0},
            {"name": "dup", "smiles": "CCOP(=O)(OCC)OCC", "charge": 0},
            {"name": "anion", "smiles": "CC(=O)[O-]", "charge": -1,
             "multiplicity": 1},
        ], "explanation": "test"}
        out = validate_molecules(payload, max_molecules=10)
        names = [m["name"] for m in out["molecules"]]
        self.assertIn("TEP-F0", names)
        self.assertIn("anion", names)
        self.assertNotIn("bad", names)
        self.assertEqual(len([n for n in names if n == "dup"]), 0)  # deduped
        self.assertEqual(out["molecules"][1]["charge"], -1)
        self.assertGreaterEqual(len(out["errors"]), 1)

    def test_missing_field_raises(self):
        with self.assertRaises(LlmError):
            validate_molecules({"nope": 1})

    def test_extract_json_fenced(self):
        text = '```json\n{"molecules": [{"name": "a", "smiles": "O"}]}\n```'
        payload = _extract_json(text)
        self.assertIn("molecules", payload)

    def test_chat_without_key_raises(self):
        cfg = _tmp_config()
        cfg.set("llm.api_key", "")
        cfg.save()
        with self.assertRaises(LlmError):
            chat(cfg, [{"role": "user", "content": "hi"}])

    def test_config_llm_defaults(self):
        cfg = _tmp_config()
        self.assertEqual(cfg.get("llm.provider"), "deepseek")
        self.assertIn("deepseek", (cfg.get("llm.base_url") or ""))
        self.assertEqual(cfg.get("llm.model"), "deepseek-v4-flash")
        self.assertEqual(cfg.get("llm.model_options"),
                         ["deepseek-v4-flash", "deepseek-v4-pro"])
        self.assertEqual(Path(cfg.home), Path(cfg.path).resolve().parent)

    def test_model_aliases(self):
        self.assertEqual(llm.resolve_model("v4pro"), "deepseek-v4-pro")
        self.assertEqual(llm.resolve_model("V4 Flash"), "deepseek-v4-flash")
        self.assertEqual(llm.resolve_model("deepseek-chat"), "deepseek-v4-flash")
        self.assertEqual(llm.model_label("v4pro"), llm.MODEL_LABELS["deepseek-v4-pro"])
        self.assertEqual(llm.model_from_label("V4 Pro（最强 / 较慢）"),
                         "deepseek-v4-pro")

    def test_old_config_migration(self):
        tmp = Path(tempfile.mkdtemp(prefix="gflow_cfg_"))
        try:
            path = tmp / "config.json"
            path.write_text(json.dumps({
                "config_version": 3,
                "llm": {"model": "deepseek-chat", "api_key": ""}}),
                encoding="utf-8")
            cfg = Config.load(path)
            self.assertEqual(cfg.get("llm.model"), "deepseek-v4-flash")
            self.assertEqual(cfg.get("llm.model_options"),
                             ["deepseek-v4-flash", "deepseek-v4-pro"])
            self.assertEqual(Path(cfg.home), tmp.resolve())
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestGjfBatch(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="gflow_gjf_"))
        self.cfg = Config.load(self.tmp / "config.json")
        self.store = Store(self.cfg.db_path)
        self.gjf1 = self.tmp / "water1.gjf"
        self.gjf2 = self.tmp / "water2.gjf"
        self.gjf1.write_text(WATER_GJF, encoding="utf-8")
        self.gjf2.write_text(WATER_GJF, encoding="utf-8")

    def tearDown(self):
        self.store.close()
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_gjf_study_spec(self):
        spec = gjf_study_spec(self.gjf1)
        self.assertEqual(spec["workflow"], "opt_freq")
        self.assertEqual(spec["charge"], 0)
        self.assertEqual(spec["multiplicity"], 1)
        self.assertEqual(spec["natoms"], 3)
        self.assertIn("B3LYP/6-31G(d)", " ".join(spec["route_tokens"]))

    def test_create_batch_from_gjfs(self):
        res = create_batch_from_gjfs(
            self.store, self.cfg, [str(self.gjf1), str(self.gjf2)],
            name_prefix="manual", submit=False)
        self.assertTrue(res["ok"])
        self.assertEqual(res["n_studies"], 2)
        studies = self.store.list_studies(batch_id=res["batch_id"])
        self.assertEqual(len(studies), 2)
        self.assertTrue(all(s["workflow"] == "opt_freq" for s in studies))
        # each study has a copied initial.gjf
        for s in studies:
            self.assertTrue((self.cfg.work_root / s["id"] / "initial.gjf").exists())


class TestPurgeCancelled(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="gflow_purge_"))
        self.cfg = Config.load(self.tmp / "config.json")
        self.store = Store(self.cfg.db_path)

    def tearDown(self):
        self.store.close()
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_purge_removes_only_cancelled(self):
        from gaussflow.pipeline import purge_cancelled
        keep = self.store.create_study(name="keep", charge=0, multiplicity=1,
                                       route="HF/3-21G", workflow="sp",
                                       status="completed")
        gone = self.store.create_study(name="gone", charge=0, multiplicity=1,
                                       route="HF/3-21G", workflow="sp",
                                       status="cancelled")
        self.store.create_job(study_id=gone["id"], stage="sp", seq=0,
                              status="cancelled")
        res = purge_cancelled(self.store, self.cfg, batch_ids=None)
        self.assertEqual(res["count"], 1)
        self.assertIsNotNone(self.store.get_study(keep["id"]))
        self.assertIsNone(self.store.get_study(gone["id"]))


class TestFluorination(unittest.TestCase):
    def test_ethanol_sites_are_labelled(self):
        sites = fluorinate.describe_fluorination_sites("CCO")
        self.assertEqual(len(sites), 3)
        kinds = {s["site"]: s["site_type"] for s in sites}
        self.assertEqual(kinds["C1"], "CH3 端基")
        self.assertEqual(kinds["C2"], "CH2")
        self.assertEqual(kinds["O3"], "OH")
        self.assertTrue(any("CH3" in s["label"] for s in sites))
        self.assertTrue(any("CH2" in s["label"] for s in sites))

    def test_propane_terminal_and_middle_carbon(self):
        sites = fluorinate.describe_fluorination_sites("CCC")
        # C1 and C3 are symmetry-equivalent, so one terminal + one CH2 remain
        self.assertEqual(len(sites), 2)
        self.assertEqual({s["site_type"] for s in sites}, {"CH3 端基", "CH2"})

    def test_heteroaromatic_site_label(self):
        sites = fluorinate.describe_fluorination_sites("c1cc[nH]c1")
        self.assertTrue(any("芳环 NH" in s["site_type"] for s in sites))
        self.assertTrue(any("芳环 CH" in s["site_type"] for s in sites))

    def test_multi_f_is_deterministic_and_labelled(self):
        first = fluorinate.enumerate_fluorinated_variants("CCO", n_f=2)
        second = fluorinate.enumerate_fluorinated_variants("CCO", n_f=2)
        self.assertEqual([v["smiles"] for v in first],
                         [v["smiles"] for v in second])
        self.assertEqual(len(first), 5)
        self.assertTrue(all(v["site_counts"] for v in first))
        self.assertTrue(any("×" in v["label"] or "+" in v["label"] for v in first))
        for variant in first:
            mol = Chem.MolFromSmiles(variant["smiles"])
            self.assertEqual(sum(1 for a in mol.GetAtoms()
                                 if a.GetSymbol() == "F"), 2)

    def test_zero_to_three_series(self):
        base = {"name": "EC", "smiles": "CCO", "formula": "C2H6O",
                "charge": 0, "multiplicity": 1, "natoms": 9}
        molecules, notes = fluorinate.molecules_for_f_counts(
            base, [0, 1, 2], max_variants=20)
        self.assertEqual(notes, [])
        self.assertEqual(len(molecules), 1 + 3 + 5)
        by_f = {}
        for mol in molecules:
            by_f.setdefault(mol["n_F"], []).append(mol)
        self.assertEqual(len(by_f[0]), 1)
        self.assertEqual(by_f[0][0]["smiles"], "CCO")
        self.assertEqual(len(by_f[1]), 3)
        self.assertEqual(len(by_f[2]), 5)
        self.assertEqual(len({m["name"] for m in molecules}), len(molecules))

    def test_too_many_f_raises(self):
        with self.assertRaises(ValueError):
            fluorinate.enumerate_fluorinated_variants("C", n_f=5)


class TestAiGjfExport(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="gflow_aigjf_"))
        self.cfg = Config.load(self.tmp / "config.json")
        self.store = Store(self.cfg.db_path)

    def tearDown(self):
        self.store.close()
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_one_gjf_per_molecule(self):
        molecules = [
            {"name": "a", "smiles": "CCO", "formula": "C2H6O", "n_F": 0,
             "charge": 0, "multiplicity": 1},
            {"name": "b", "smiles": "CCF", "formula": "C2H5F", "n_F": 1,
             "charge": 0, "multiplicity": 1},
            {"name": "c", "smiles": "FC(F)F", "formula": "CHF3", "n_F": 3,
             "charge": 0, "multiplicity": 1},
        ]
        res = pipeline.export_gjf_previews(
            self.cfg, molecules, recipe_name="opt_freq", name_prefix="AI")
        self.assertTrue(res["ok"])
        self.assertEqual(len(res["files"]), 3)
        for item in res["files"]:
            path = Path(item["gjf"])
            self.assertTrue(path.exists())
            self.assertIn("#p", path.read_text(encoding="ascii"))
        index = json.loads(Path(res["index"]).read_text(encoding="utf-8"))
        self.assertEqual(len(index["files"]), 3)


class TestResultLog(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="gflow_log_"))
        self.cfg = Config.load(self.tmp / "config.json")
        self.store = Store(self.cfg.db_path)

    def tearDown(self):
        self.store.close()
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _summary(self, sid="st_test_1"):
        return {
            "ok": True, "study_id": sid, "name": "ethanol-neutral",
            "batch_id": "bt_test", "status": "completed", "role": "neutral",
            "smiles": "CCO", "formula": "C2H6O", "natoms": 9,
            "charge": 0, "multiplicity": 1, "route": "B3LYP/6-31G(d)",
            "workflow": "opt_freq", "solvent": None,
            "variant": {"type": "neutral", "label": "neutral"},
            "energy_hartree": -154.0,
            "stage_results": {"optfreq": {
                "energy_hartree": -154.0, "freeenergy_hartree": -153.9,
                "homo_ev": -7.0, "lumo_ev": 1.0, "gap_ev": 8.0,
                "zpe_hartree": 0.08, "enthalpy_hartree": -153.95,
                "entropy_cal_mol_k": 60.0}},
            "orbitals": {"homo_ev": -7.0, "lumo_ev": 1.0, "gap_ev": 8.0},
            "n_imag": 0, "validation_overall": "pass",
        }

    def test_append_and_recover(self):
        pipeline._log_study_result(self.cfg, self._summary())
        missing = resultlog.missing_studies(self.store, self.cfg)
        self.assertEqual(missing["count"], 1)
        self.assertEqual(resultlog.summary(self.cfg)["counts"]["study_result"], 1)
        res = resultlog.recover(self.store, self.cfg)
        self.assertEqual(res["count"], 1)
        study = self.store.get_study("st_test_1")
        self.assertIsNotNone(study)
        self.assertEqual(study["batch_id"], "bt_test")
        summary = pipeline.study_summary(self.store, self.cfg, "st_test_1")
        self.assertTrue(summary.get("recovered"))
        self.assertEqual(summary["energy_hartree"], -154.0)
        self.assertEqual(summary["status"], "completed")
        self.assertIn("optfreq", summary["stage_results"])
        # idempotent: a second recovery skips the existing study
        again = resultlog.recover(self.store, self.cfg)
        self.assertEqual(again["count"], 0)
        self.assertEqual(again["skipped"], 1)
        self.assertEqual(resultlog.missing_studies(self.store, self.cfg)["count"], 0)
        # the result page aggregation must also work from the recovered summary
        page = analysis.system_summary(self.store, self.cfg, "bt_test")
        self.assertTrue(page["ok"])
        self.assertEqual(page["n_studies"], 1)
        self.assertEqual(len(page["orbitals"]), 1)
        self.assertEqual(page["orbitals"][0]["HOMO_eV"], -7.0)


class TestAiRecipe(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="gflow_airecipe_"))
        self.cfg = Config.load(self.tmp / "config.json")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _recipe(self, route="B3LYP/6-311+G(d) empiricaldispersion=gd3bj"):
        return {
            "name": "AI 测试链",
            "description": "单测生成的计算链",
            "solvent": "DMSO",
            "core_stages": [
                {"name": "optfreq", "kind": "optfreq", "route": route,
                 "phase": "gas"},
                {"name": "sp_sol", "kind": "sp", "route": "M06-2X/6-31+G(d)",
                 "phase": "solvent", "solvent": "DMSO", "pop": "full"},
            ],
            "li": {"on": False},
            "redox": {"on": False},
        }

    def test_valid_recipe_is_registered_and_normalised(self):
        res = pipeline.register_ai_recipe(self.cfg, self._recipe())
        self.assertTrue(res["ok"], res.get("problems"))
        self.assertIn("ai_", res["key"])
        routes = [s["route"] for s in res["recipe"]["core_stages"]]
        self.assertIn("M062X/6-31+G(d)", routes)
        self.assertIn(res["key"], (self.cfg.get("recipes") or {}))
        saved = json.loads(Path(self.cfg.path).read_text(encoding="utf-8"))
        self.assertIn(res["key"], saved["recipes"])

    def test_bad_basis_is_rejected(self):
        res = pipeline.register_ai_recipe(
            self.cfg, self._recipe(route="B3LYP/not-a-real-basis"))
        self.assertFalse(res["ok"])
        self.assertTrue(res["problems"])

    def test_scrf_inside_route_is_rejected(self):
        res = pipeline.register_ai_recipe(
            self.cfg, self._recipe(
                route="B3LYP/6-31G(d) scrf=(pcm,solvent=DMSO)"))
        self.assertFalse(res["ok"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
