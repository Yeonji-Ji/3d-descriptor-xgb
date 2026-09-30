"""Steps 4-5: per-seed energies, best seed per combination, actinyl QC.

conda activate env-pd && python xtb-complex-gxtb-actinyl/analyze.py

Writes seeds.csv, structures.csv, best/<idx>_<El><ox>.xyz, and prints which combinations still have
their best seed in 15-20 (they get one extra round of seeds 21-30).

Flags: the +4 set's rules (scf_fail, not_converged, cn_change, short_bond, water_left, missing,
spin_mismatch, late_best_seed) plus the actinyl checks:
  oxo_angle     O=M=O below 170 deg
  oxo_dist      an M=O outside 1.70-1.90 A
  eq_out_of_plane   an equatorial donor more than 25 deg off the equatorial plane
  eq_CN_range   equatorial CN outside 4-6
Atom order (from build_list.py): ligand heavy, metal, oxo O x2, water O, then H.
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent / "xtb-complex-gxtb"))
import qc as base  # noqa: E402  (read_xyz, parse_log, SHORT_BOND)

SPIN_TOL, LATE_WINDOW = 0.2, (15, 20)
OXO_ANGLE_MIN, OXO_RANGE, EQ_DEV_MAX, EQ_CN_RANGE = 170.0, (1.70, 1.90), 25.0, (4, 6)
NO_CN, S_EDGE = 3.2, 3.3          # same coordination rule as xtb-complex-1/pd


def geometry(el, X, r):
    L = int(r["n_heavy"])
    assert el[L] == r["element"], f"{r['name']}: metal not at index {L}"
    M = X[L]
    oxo = [L + 1, L + 2]
    wat = list(range(L + 3, L + 3 + int(r["n_water_eq"])))
    don = [int(x) for x in r["donors"].split()]
    assert all(el[i] == "O" for i in oxo + wat), "oxo/water indices do not hold O"
    d = np.linalg.norm(X - M, axis=1)

    v = [(X[i] - M) / d[i] for i in oxo]
    ang = float(np.degrees(np.arccos(np.clip(np.dot(v[0], v[1]), -1, 1))))
    axis = v[0] - v[1]
    axis /= np.linalg.norm(axis)

    def oop(i):   # angle away from the equatorial plane
        return float(abs(90 - np.degrees(np.arccos(np.clip(abs(np.dot((X[i] - M) / d[i], axis)), 0, 1)))))

    eq = [i for i in range(len(el)) if i not in oxo and i != L
          and ((el[i] in ("N", "O") and d[i] < NO_CN) or (el[i] == "S" and d[i] < S_EDGE))]
    eq_O = [i for i in eq if el[i] == "O"]
    eq_N = [i for i in eq if el[i] == "N"]
    return {"oxo_d1": float(d[oxo[0]]), "oxo_d2": float(d[oxo[1]]), "oxo_angle": ang,
            "eq_CN": len(eq), "eq_CN_O": len(eq_O), "eq_CN_N": len(eq_N),
            "M_O_mean": float(d[eq_O].mean()) if eq_O else np.nan,
            "M_N_mean": float(d[eq_N].mean()) if eq_N else np.nan,
            "M_donor_min": float(d[don].min()), "M_donor_max": float(d[don].max()),
            "donor_oop_max": max((oop(i) for i in don), default=np.nan),
            "n_donor_bound": int(sum(1 for i in don if d[i] < (S_EDGE if el[i] == "S" else NO_CN))),
            "n_water_bound": int(sum(1 for i in wat if d[i] < NO_CN)),
            "water_oop_max": max((oop(i) for i in wat), default=np.nan),
            "M_water_min": float(d[wat].min()) if wat else np.nan,
            "M_water_max": float(d[wat].max()) if wat else np.nan}


def wall_times():
    import re
    w = {}
    f = HERE / "runs" / "status.txt"
    if f.exists():
        for m in re.finditer(r"^(\S+) .*? wall_s=(\d+) skipped=(\d)", f.read_text(), re.M):
            if m.group(3) == "0":
                w[m.group(1)] = int(m.group(2))
    return w


def main():
    lst = pd.read_csv(HERE / "list.csv")
    wall = wall_times()
    (HERE / "best").mkdir(exist_ok=True)
    seeds, rows = [], []

    for r in lst.to_dict("records"):
        mine = []
        for d in sorted((HERE / "runs").glob(f"{r['name']}_s*")):
            row = {"name": r["name"], "idx": r["idx"], "element": r["element"], "ox": r["ox"],
                   "seed": int(d.name.rsplit("_s", 1)[1]), "wall_s": wall.get(d.name)}
            row.update(base.parse_log(d / "xtb.log", r["element"]) if (d / "xtb.log").exists()
                       else {"converged": False, "energy_Eh": np.nan})
            mine.append(row)
        seeds += mine
        s = pd.DataFrame(mine)
        ok = s[s.converged & s.energy_Eh.notna()].sort_values("energy_Eh") if len(s) else s
        row = {k: r[k] for k in ("name", "idx", "element", "ox", "Z", "uhf", "chrg", "SMILES",
                                 "family", "splits", "denticity", "n_water_eq", "n_atoms",
                                 "flag_build")}
        row.update({"n_seeds": len(s), "n_converged": len(ok),
                    "wall_s_total": float(s.wall_s.dropna().sum()) if len(s) else 0.0,
                    "n_scf_fail": int(s.scf_fail.sum()) if len(s) else 0})
        flags = []
        if len(ok):
            b = ok.iloc[0]
            row.update({"best_seed": int(b.seed), "energy_Eh": float(b.energy_Eh),
                        "gap_to_second_kcal": (float(ok.iloc[1].energy_Eh - b.energy_Eh) * 627.509
                                               if len(ok) > 1 else np.nan),
                        "M_q": b.M_q, "M_spin": b.M_spin, "spin_excess": float(b.M_spin) - r["uhf"]})
            src = HERE / "runs" / f"{r['name']}_s{int(b.seed):02d}" / "xtbopt.xyz"
            shutil.copy(src, HERE / "best" / f"{r['name']}.xyz")
            g0 = geometry(*base.read_xyz(HERE / "starts" / f"{r['name']}_s{int(b.seed):02d}.xyz"), r)
            g = geometry(*base.read_xyz(src), r)
            row.update({f"start_{k}": v for k, v in g0.items()})
            row.update(g)
            if g["eq_CN"] != g0["eq_CN"]:
                flags.append("cn_change")
            if g["M_donor_min"] < base.SHORT_BOND:
                flags.append("short_bond")
            if r["n_water_eq"] and g["n_water_bound"] < r["n_water_eq"]:
                flags.append("water_left")
            if abs(row["spin_excess"]) > SPIN_TOL:
                flags.append("spin_mismatch")
            if g["oxo_angle"] < OXO_ANGLE_MIN:
                flags.append("oxo_angle")
            if not all(OXO_RANGE[0] <= g[k] <= OXO_RANGE[1] for k in ("oxo_d1", "oxo_d2")):
                flags.append("oxo_dist")
            if g["donor_oop_max"] > EQ_DEV_MAX:
                flags.append("eq_out_of_plane")
            if not EQ_CN_RANGE[0] <= g["eq_CN"] <= EQ_CN_RANGE[1]:
                flags.append("eq_CN_range")
            if LATE_WINDOW[0] <= b.seed <= LATE_WINDOW[1]:
                flags.append("late_best_seed")
        else:
            flags.append("missing")
            if len(s) and s.scf_fail.any():
                flags.append("scf_fail")
            if len(s) and not s.converged.any():
                flags.append("not_converged")
        row["flags"] = ";".join(flags)
        row["ok"] = not [f for f in flags if f != "late_best_seed"]
        rows.append(row)

    sd = pd.DataFrame(seeds)
    sd.to_csv(HERE / "seeds.csv", index=False)
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "structures.csv", index=False)

    print(f"combinations {len(df)}; with a structure {int(df.n_converged.gt(0).sum())}; clean {int(df.ok.sum())}")
    if len(sd):
        print(f"seed runs {len(sd)}; converged {int(sd.converged.sum())} ({sd.converged.mean():.1%}); "
              f"wall {sd.wall_s.dropna().sum()/3600:.1f} h")
    print("\nflag counts:")
    for f in ("scf_fail", "not_converged", "missing", "cn_change", "short_bond", "water_left",
              "spin_mismatch", "oxo_angle", "oxo_dist", "eq_out_of_plane", "eq_CN_range", "late_best_seed"):
        n = int(df["flags"].str.contains(f).sum())
        if n:
            print(f"  {f:16s} {n}")
    late = df[df["flags"].str.contains("late_best_seed")]
    print(f"\nbest seed in {LATE_WINDOW}: {len(late)}")
    if len(late):
        print(late[["name", "best_seed", "n_converged", "gap_to_second_kcal"]].to_string(index=False))
    print("\nby metal:")
    print(df.groupby(["element", "ox", "uhf"]).agg(
        n=("name", "size"), clean=("ok", "sum"), oxo_angle=("oxo_angle", "mean"),
        oxo_d=("oxo_d1", "mean"), eq_CN=("eq_CN", "mean"), spin_excess=("spin_excess", "mean"),
        M_q=("M_q", "mean")).round(3).to_string())


if __name__ == "__main__":
    main()
