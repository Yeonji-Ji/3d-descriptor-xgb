"""Compare the constrained runs with the parent unconstrained runs of the same seeds.

conda activate env-pd && python xtb-complex-gxtb-actinyl/constrain-test/compare.py

Geometry and log parsing are reused from the parent analyze.py (which itself reuses
xtb-complex-gxtb/qc.py), so the numbers are defined exactly as in the parent QC.
Writes results.csv (one row per run) and prints the per-complex comparison.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(PARENT))
import analyze as pa  # noqa: E402  (geometry, base.read_xyz, base.parse_log)

SEEDS = range(1, 5)


def one(run_dir, row, mode, seed, target):
    out = {"name": row["name"], "element": row["element"], "ox": row["ox"], "seed": seed,
           "mode": mode, "angle_target": target}
    log = run_dir / "xtb.log"
    if not log.exists():
        return {**out, "converged": False, "note": "no run"}
    out.update(pa.base.parse_log(log, row["element"]))
    opt = run_dir / "xtbopt.xyz"
    if opt.exists():
        el, X = pa.base.read_xyz(opt)
        g = pa.geometry(el, X, row)
        out.update({k: g[k] for k in ("oxo_angle", "oxo_d1", "oxo_d2", "eq_CN", "donor_oop_max",
                                      "water_oop_max", "M_donor_min", "n_water_bound")})
        out["spin_excess"] = out["M_spin"] - row["uhf"] if out.get("M_spin") == out.get("M_spin") else np.nan
    return out


def main():
    lst = pd.read_csv(PARENT / "list.csv").set_index("name")
    rows = []
    for d in sorted(HERE.glob("runs*/*_s??")):
        name, seed = d.name.rsplit("_s", 1)
        target = 178.0 if d.parent.name.endswith("178") else 180.0
        rows.append(one(d, {**lst.loc[name].to_dict(), "name": name}, f"constrained_{target:.0f}",
                        int(seed), target))
        if target == 180.0:   # matching unconstrained run from the parent set
            rows.append(one(PARENT / "runs" / d.name, {**lst.loc[name].to_dict(), "name": name},
                            "unconstrained", int(seed), np.nan))
    df = pd.DataFrame(rows).sort_values(["name", "seed", "mode"])
    df.to_csv(HERE / "results.csv", index=False)

    pd.set_option("display.width", 250)
    cols = ["oxo_angle", "oxo_d1", "oxo_d2", "eq_CN", "donor_oop_max", "M_q", "spin_excess",
            "final_grad", "energy_Eh", "converged"]
    print(df[["name", "seed", "mode"] + cols].round(3).to_string(index=False))

    print("\nper complex, lowest-energy converged run of each mode:")
    best = (df[df.converged & df.energy_Eh.notna()].sort_values("energy_Eh")
              .groupby(["name", "mode"], as_index=False).first())
    print(best[["name", "mode", "seed"] + cols].round(3).to_string(index=False))
    piv = best.pivot(index="name", columns="mode", values=["oxo_angle", "oxo_d1", "final_grad", "energy_Eh"])
    if ("energy_Eh", "constrained_180") in piv and ("energy_Eh", "unconstrained") in piv:
        piv[("dE_kcal", "")] = (piv[("energy_Eh", "constrained_180")]
                                - piv[("energy_Eh", "unconstrained")]) * 627.509
    print("\n", piv.round(3).to_string())


if __name__ == "__main__":
    main()
