"""Constrained-axis test: re-run 5 actinyl complexes with O=M=O held at 180 deg.

conda activate env-gxtb && python xtb-complex-gxtb-actinyl/constrain-test/run_constrained.py [--angle 180.0]

Same start structures as the parent run (../starts/), seeds 1-4, same charge/uhf, gas phase.
Per complex an xcontrol file is written with the verified 1-based atom indices:

    $constrain
      force constant=1.0
      angle: <O1>,<M>,<O2>,<target>
    $end

Runs in runs/<name>_s<seed>/ (and runs_178/<...> for the 178 deg retry). Parent folder is read-only.
"""
import argparse
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
CASES = ["17_U6", "259_U6", "280_Np5", "12_Pu6", "3_U6"]
SEEDS = range(1, 5)
NPROC = 9


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    el = np.array([l.split()[0] for l in lines[2:2 + n]])
    X = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
    return el, X


def indices(xyz, row):
    """Verify the builder's atom order in this file: metal at n_heavy, then two oxo O (0-based)."""
    el, X = read_xyz(xyz)
    L = int(row.n_heavy)
    assert el[L] == row.element, f"{xyz}: expected {row.element} at index {L}, found {el[L]}"
    o1, o2 = L + 1, L + 2
    assert el[o1] == "O" and el[o2] == "O", f"{xyz}: oxo positions hold {el[o1]},{el[o2]}"
    d = [np.linalg.norm(X[o] - X[L]) for o in (o1, o2)]
    assert all(1.5 < x < 2.1 for x in d), f"{xyz}: oxo distances {d} do not look like M=O"
    return L, o1, o2, d


def run(task):
    name, seed, row, angle, out_root = task
    d = out_root / f"{name}_s{seed:02d}"
    d.mkdir(parents=True, exist_ok=True)
    src = PARENT / "starts" / f"{name}_s{seed:02d}.xyz"
    L, o1, o2, dist = indices(src, row)
    (d / "in.xyz").write_text(src.read_text())
    # xcontrol uses 1-based atom numbering
    (d / "constrain.inp").write_text(
        "$constrain\n  force constant=1.0\n"
        f"  angle: {o1 + 1},{L + 1},{o2 + 1},{angle:.1f}\n$end\n")
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_STACKSIZE="4G")
    cmd = ["xtb", "in.xyz", "--gxtb", "--opt", "--chrg", str(int(row.chrg)),
           "--uhf", str(int(row.uhf)), "--input", "constrain.inp"]
    t0 = time.time()
    with open(d / "xtb.log", "w") as log:
        rc = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
    wall = time.time() - t0
    txt = (d / "xtb.log").read_text(errors="replace")
    ok = "GEOMETRY OPTIMIZATION CONVERGED" in txt
    print(f"{out_root.name}/{name}_s{seed:02d}: rc={rc} converged={ok} {wall:.0f}s "
          f"(constraint {o1 + 1},{L + 1},{o2 + 1} -> {angle:.1f}, start M=O {dist[0]:.2f}/{dist[1]:.2f})",
          flush=True)
    return {"name": name, "seed": seed, "angle_target": angle, "rc": rc, "converged": ok,
            "wall_s": wall, "dir": str(d.relative_to(HERE))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--angle", type=float, default=180.0)
    ap.add_argument("--out", default="runs")
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    if "env-gxtb" not in sys.executable:
        sys.exit("activate env-gxtb (the g-xTB binary must be on PATH)")

    lst = pd.read_csv(PARENT / "list.csv").set_index("name")
    cases = args.only or CASES
    tasks = [(n, s, lst.loc[n], args.angle, HERE / args.out) for n in cases for s in SEEDS]
    with ThreadPoolExecutor(NPROC) as ex:
        rows = list(ex.map(run, tasks))
    out = HERE / f"run_log_{args.out}.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\n{sum(r['converged'] for r in rows)}/{len(rows)} converged -> {out}")


if __name__ == "__main__":
    main()
