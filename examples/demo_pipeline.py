"""End-to-end demo of the gaussflow tool interface.

Runs the whole chain through the same JSON CLI the Agent uses:
  doctor -> iterate (base + Na+/Cl- variants) -> submit -> poll -> export

Usage:
    cd gaussflow
    python examples/demo_pipeline.py            # full demo (real Gaussian jobs)
    python examples/demo_pipeline.py --dry-run  # just prepare, do not submit
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable


def gflow(*args: str) -> dict:
    cmd = [PY, str(ROOT / "gflow.py"), *args, "--json"]
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                          errors="replace", cwd=str(ROOT))
    text = proc.stdout.strip()
    try:
        data = json.loads(text) if text else {}
    except json.JSONDecodeError:
        raise RuntimeError(f"gflow {' '.join(args)} 输出不是 JSON:\n{text[:800]}")
    if proc.returncode not in (0, 3) and not data.get("ok", False):
        raise RuntimeError(f"gflow {' '.join(args)} 失败: {data}")
    return data


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles", default="O", help="基准分子（默认水）")
    ap.add_argument("--name", default="water-ion-demo")
    ap.add_argument("--route", default="B3LYP/6-31G(d)")
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--ions", default="Na+,Cl-")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    print("== 1) doctor ==")
    doc = gflow("doctor", "--test-run")
    print(f"   g16 ok: {doc['checks'].get('g16_test', {}).get('ok')} "
          f"energy={doc['checks'].get('g16_test', {}).get('energy_hartree')}")

    print("== 2) iterate ==")
    batch = gflow("iterate", "--smiles", args.smiles, "--name", args.name,
                  "--n", str(args.n), "--ions", args.ions,
                  "--route", args.route, "--workflow", "opt_freq",
                  "--seed", str(args.seed))
    bid = batch["batch_id"]
    print(f"   batch={bid} variants={batch['count']}")
    for s in batch["studies"]:
        print(f"   - {s['study_id']} {s['name']} variant={s['variant']}")
    if args.dry_run:
        print("dry-run: 不提交计算。可用 gflow submit-batch --batch", bid, "提交")
        return 0

    print("== 3) submit ==")
    sub = gflow("submit-batch", "--batch", bid)
    print(f"   submitted, worker={sub.get('worker')}")

    print("== 4) poll ==")
    deadline = time.time() + 3600
    last = ""
    while time.time() < deadline:
        st = gflow("status", "--batch", bid)
        states = " ".join(f"{s['variant']}:{s['status']}" for s in st["studies"])
        if states != last:
            print("   ", states)
            last = states
        if all(s["status"] in ("completed", "failed") for s in st["studies"]):
            break
        time.sleep(10)

    print("== 5) results ==")
    for s in gflow("status", "--batch", bid)["studies"]:
        r = gflow("result", "--study", s["study_id"])
        print(f"   {r['name']:<28} {r['status']:<10} "
              f"E={r.get('energy_hartree')} Ha  NImag={r.get('n_imag')}  "
              f"valid={r.get('validation_overall')}")

    print("== 6) export ==")
    exp = gflow("export", "--batch", bid)
    print("   out_dir:", exp.get("out_dir"))
    print("   csv:", exp.get("results_csv"))
    print("   xlsx:", exp.get("summary_xlsx"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
