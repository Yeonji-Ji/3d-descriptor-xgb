"""Step 3: QC every g-xTB run and write structures.csv (one row per combination).

conda activate env-pd && python xtb-complex-gxtb/qc.py

Flags (a run can carry several):
  scf_fail        the log reports an SCF that did not converge
  not_converged   no "GEOMETRY OPTIMIZATION CONVERGED" in the log (or the run never produced xtbopt.xyz)
  cn_change       coordination number differs from the start structure
  short_bond      a metal-donor distance below 2.0 A
  water_left      a water O further than 3.2 A from the metal
  missing         no run folder / no output yet
Coordination rule as in xtb-complex-1/pd: N/O within 3.2 A, S within 3.3 A.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
XTB = HERE.parent / "xtb-complex-1"
sys.dont_write_bytecode = True
sys.path.insert(0, str(XTB / "pd"))
import compute_pd as cpd  # noqa: E402  (NO_CN = 3.2, S_EDGE = 3.3)
sys.path.insert(0, str(HERE))
from build_list import SYMBOL  # noqa: E402

Z_OF = {v: k for k, v in SYMBOL.items()}
SHORT_BOND = 2.0


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    el = np.array([l.split()[0] for l in lines[2:2 + n]])
    X = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
    return el, X


def geometry(el, X, r):
    L = int(r.n_heavy)
    assert el[L] == r.element, f"{r.name}: metal not at index {L}"
    d = np.linalg.norm(X - X[L], axis=1)
    isN, isO, isS = el == "N", el == "O", el == "S"
    cN, cO, cS = isN & (d < cpd.NO_CN), isO & (d < cpd.NO_CN), isS & (d < cpd.S_EDGE)
    don = np.array([int(x) for x in r.donors.split()])
    wat = np.arange(L + 1, L + 1 + int(r.n_water))
    return {"CN": int(cN.sum() + cO.sum() + cS.sum()),
            "M_O_mean": float(d[cO].mean()) if cO.any() else np.nan,
            "M_N_mean": float(d[cN].mean()) if cN.any() else np.nan,
            "M_donor_min": float(d[don].min()), "M_donor_max": float(d[don].max()),
            "n_donor_bound": int((d[don] < np.where(el[don] == "S", cpd.S_EDGE, cpd.NO_CN)).sum()),
            "n_water_bound": int((d[wat] < cpd.NO_CN).sum()) if len(wat) else 0,
            "M_water_min": float(d[wat].min()) if len(wat) else np.nan,
            "M_water_max": float(d[wat].max()) if len(wat) else np.nan}


def parse_log(path, element):
    txt = Path(path).read_text(errors="replace")
    e = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", txt)
    g = re.findall(r"::\s+gradient norm\s+(\S+)\s+Eh/a0", txt)
    out = {"converged": "GEOMETRY OPTIMIZATION CONVERGED" in txt,
           "scf_fail": "SCF not converged" in txt,
           "error": (re.findall(r"^\s*-1-\s*(.+?)\s*$", txt, re.M) or [""])[0][:80],
           "opt_cycles": len(re.findall(r"CYCLE\s+\d+\s", txt)),
           "final_grad": float(g[-1]) if g else np.nan,
           "energy_Eh": float(e[-1]) if e else np.nan}
    pat = rf"^\s*\d+\s+{Z_OF[element]} {element}\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)"
    for head, keys in (("Atomic charges and shell populations", ("M_q", "M_pf")),
                       ("Atomic and shell magnetization", ("M_spin", "M_spin_f"))):
        m = re.search(pat, txt.rsplit(head, 1)[1], re.M) if head in txt else None
        out[keys[0]], out[keys[1]] = (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)
    return out


def wall_times():
    """Last status line per run (the file is appended across restarts)."""
    w = {}
    f = HERE / "runs" / "status.txt"
    if f.exists():
        for ln in f.read_text().split("\n"):
            m = re.match(r"(\S+) uhf=(\d+)(?: chrg=(-?\d+))? rc=(\d+) wall_s=(\d+) skipped=(\d)", ln)
            if m:
                w[m.group(1)] = {"rc": int(m.group(4)), "wall_s": int(m.group(5)),
                                 "skipped": int(m.group(6))}
    return w


def main():
    lst = pd.read_csv(HERE / "list.csv")
    wall = wall_times()
    rows = []
    for r in lst.itertuples():
        d = HERE / "runs" / r.name
        row = {k: getattr(r, k) for k in ("name", "idx", "element", "Z", "uhf", "SMILES", "family",
                                          "group", "splits", "denticity", "n_water", "n_atoms",
                                          "gfn2_source", "gfn2_seed")}
        row.update(wall.get(r.name, {"rc": None, "wall_s": None, "skipped": None}))
        el0, X0 = read_xyz(HERE / "starts" / f"{r.name}.xyz")
        g0 = geometry(el0, X0, r)
        row.update({f"start_{k}": v for k, v in g0.items()})
        flags = []
        if not (d / "xtb.log").exists():
            flags.append("missing")
            row["converged"] = False
        else:
            row.update(parse_log(d / "xtb.log", r.element))
            if row["scf_fail"]:
                flags.append("scf_fail")
            if not row["converged"]:
                flags.append("not_converged")
            if (d / "xtbopt.xyz").exists():
                el, X = read_xyz(d / "xtbopt.xyz")
                g = geometry(el, X, r)
                row.update(g)
                if g["CN"] != g0["CN"]:
                    flags.append("cn_change")
                if g["M_donor_min"] < SHORT_BOND:
                    flags.append("short_bond")
                if r.n_water and g["n_water_bound"] < r.n_water:
                    flags.append("water_left")
            else:
                flags.append("missing")
        row["flags"] = ";".join(flags)
        row["ok"] = not flags
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(HERE / "structures.csv", index=False)
    done = df[df.wall_s.notna()]
    print(f"combinations {len(df)}; with a result {(~df['flags'].str.contains('missing')).sum()}; "
          f"clean {int(df.ok.sum())}; flagged {int((~df.ok).sum())}")
    print("flag counts:")
    for f in ("scf_fail", "not_converged", "cn_change", "short_bond", "water_left", "missing"):
        n = int(df['flags'].str.contains(f).sum())
        if n:
            print(f"  {f:15s} {n}")
    bad = df[~df.ok]
    if len(bad):
        print("\nflagged by element:")
        print(pd.crosstab(bad.element, bad['flags']).to_string())
        print("\nsuccess rate by element:")
        print(df.groupby("element").ok.agg(["sum", "size"]).assign(
            rate=lambda t: (t["sum"] / t["size"]).round(3)).to_string())
    if len(done):
        print(f"\nwall time: total {done.wall_s.sum() / 3600:.1f} h, median {done.wall_s.median():.0f} s, "
              f"max {done.wall_s.max():.0f} s")


if __name__ == "__main__":
    main()
