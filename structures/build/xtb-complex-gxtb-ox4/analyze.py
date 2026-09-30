"""Steps 4-5: per-seed energies, best seed per combination, QC of the best structure.

conda activate env-pd && python xtb-complex-gxtb-ox4/analyze.py

Writes seeds.csv (one row per seed), structures.csv (one row per combination, QC of its best seed),
best/<idx>_<El>.xyz and prints which combinations need extra seeds (best seed in the 15-20 window).
QC rules and parsing are reused from xtb-complex-gxtb/qc.py; CN rule from xtb-complex-1/pd.
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent / "xtb-complex-gxtb"))
import qc as base  # noqa: E402  (read_xyz, geometry, parse_log, wall_times, SHORT_BOND)

SPIN_TOL = 0.2
LATE_WINDOW = (15, 20)   # a best seed in this window means the search may not have converged


class Row:
    """Attribute view of a list.csv record, for base.geometry()."""
    def __init__(self, d):
        self.__dict__.update(d)


def wall_times():
    """status.txt lines: <name> uhf=<n> chrg=<q> rc=<rc> wall_s=<s> skipped=<0|1>"""
    import re
    w = {}
    f = HERE / "runs" / "status.txt"
    if f.exists():
        # keep the real timing: a resumed pass re-logs skipped runs with wall_s=0
        for m in re.finditer(r"^(\S+) .*? wall_s=(\d+) skipped=(\d)", f.read_text(), re.M):
            if m.group(3) == "0":
                w[m.group(1)] = int(m.group(2))
    return w


def main():
    lst = pd.read_csv(HERE / "list.csv")
    wall = wall_times()
    seeds, best_rows = [], []
    (HERE / "best").mkdir(exist_ok=True)

    for rec in lst.to_dict("records"):
        r = Row(rec)
        runs = sorted((HERE / "runs").glob(f"{r.name}_s*"))
        for d in runs:
            seed = int(d.name.rsplit("_s", 1)[1])
            row = {"name": r.name, "idx": r.idx, "element": r.element, "seed": seed,
                   "wall_s": wall.get(d.name)}
            if (d / "xtb.log").exists():
                row.update(base.parse_log(d / "xtb.log", r.element))
            else:
                row.update({"converged": False, "energy_Eh": np.nan})
            seeds.append(row)
        s = pd.DataFrame([x for x in seeds if x["name"] == r.name])
        ok = s[s.converged & s.energy_Eh.notna()].sort_values("energy_Eh")
        row = {k: rec[k] for k in ("name", "idx", "element", "Z", "uhf", "chrg", "SMILES", "family",
                                   "splits", "denticity", "n_water", "n_atoms")}
        row.update({"n_seeds": len(s), "n_converged": len(ok),
                    "wall_s_total": float(s.wall_s.dropna().sum())})
        flags = []
        if len(ok):
            b = ok.iloc[0]
            row.update({"best_seed": int(b.seed), "energy_Eh": float(b.energy_Eh),
                        "gap_to_second_kcal": (float(ok.iloc[1].energy_Eh - b.energy_Eh) * 627.509
                                               if len(ok) > 1 else np.nan),
                        "M_q": b.M_q, "M_spin": b.M_spin, "M_spin_f": b.M_spin_f,
                        "spin_excess": float(b.M_spin) - r.uhf,
                        "spin_ok": abs(b.M_spin - r.uhf) <= SPIN_TOL})
            src = HERE / "runs" / f"{r.name}_s{int(b.seed):02d}" / "xtbopt.xyz"
            shutil.copy(src, HERE / "best" / f"{r.name}.xyz")
            el0, X0 = base.read_xyz(HERE / "starts" / f"{r.name}_s{int(b.seed):02d}.xyz")
            el, X = base.read_xyz(src)
            g0, g = base.geometry(el0, X0, r), base.geometry(el, X, r)
            row.update({f"start_{k}": v for k, v in g0.items()})
            row.update(g)
            if g["CN"] != g0["CN"]:
                flags.append("cn_change")
            if g["M_donor_min"] < base.SHORT_BOND:
                flags.append("short_bond")
            if r.n_water and g["n_water_bound"] < r.n_water:
                flags.append("water_left")
            if not row["spin_ok"]:
                flags.append("spin_mismatch")
            if LATE_WINDOW[0] <= b.seed <= LATE_WINDOW[1]:
                flags.append("late_best_seed")
        else:
            flags.append("missing")
            if s.get("scf_fail", pd.Series(dtype=bool)).any():
                flags.append("scf_fail")
            if len(s) and not s.converged.any():
                flags.append("not_converged")
        row["n_scf_fail"] = int(s.get("scf_fail", pd.Series(dtype=bool)).sum())
        row["n_not_converged"] = int((~s.converged).sum()) if len(s) else 0
        row["flags"] = ";".join(flags)
        row["ok"] = not [f for f in flags if f != "late_best_seed"]
        best_rows.append(row)

    sd = pd.DataFrame(seeds)
    sd.to_csv(HERE / "seeds.csv", index=False)
    df = pd.DataFrame(best_rows)
    df.to_csv(HERE / "structures.csv", index=False)

    print(f"combinations {len(df)}; with a converged structure {int(df.n_converged.gt(0).sum())}; "
          f"clean {int(df.ok.sum())}")
    print(f"seed runs {len(sd)}; converged {int(sd.converged.sum())} "
          f"({sd.converged.mean():.1%}); compute {sd.wall_s.dropna().sum()/3600:.1f} h")
    print("\nflags:")
    print(df[df['flags'] != ""]['flags'].value_counts().to_string() or "none")
    late = df[df['flags'].str.contains("late_best_seed")]
    print(f"\nbest seed in {LATE_WINDOW[0]}-{LATE_WINDOW[1]} (needs 10 extra seeds): {len(late)}")
    if len(late):
        print(late[["name", "best_seed", "n_converged", "gap_to_second_kcal"]].to_string(index=False))
    print("\nby element:")
    print(df.groupby("element").agg(combinations=("name", "size"), clean=("ok", "sum"),
                                    mean_converged_seeds=("n_converged", "mean")).round(2).to_string())
    print("\ngap to second best (kcal/mol):")
    print(df.gap_to_second_kcal.describe().round(2).to_string())


if __name__ == "__main__":
    main()
