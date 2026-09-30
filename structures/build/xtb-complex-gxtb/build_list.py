"""Step 1: list the unique (SMILES, atomic number) combinations with Oxidation_state == 3 and
write one g-xTB start structure per combination.

conda activate env-pd && python xtb-complex-gxtb/build_list.py

Start structure = lowest-energy GFN2 seed of that ligand in xtb-complex-1/ (best_seed in the group's
summary.csv), copied unchanged except that the Eu label becomes the target element. Waters are kept.
xtb-complex-1/ is only read. Outputs: starts/<idx>_<El>.xyz, list.csv, unmatched.csv.
"""
import glob
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "orig-data").is_dir())   # project root
XTB = HERE.parent / "xtb-complex-1"
ORIG = ROOT / "orig-data" / "0.1-1"

SYMBOL = {57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb",
          66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 89: "Ac", 90: "Th", 91: "Pa",
          92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 98: "Cf"}
# Unpaired electrons of the +3 ion: f^n has n unpaired up to the half-filled shell (n <= 7) and
# 14 - n beyond it (Hund). The first run used n throughout, which forced the surplus unpaired
# electrons onto the ligands for Tb-Yb and Cf; see README "Problems found".
UHF = {"La": 0, "Ce": 1, "Pr": 2, "Nd": 3, "Pm": 4, "Sm": 5, "Eu": 6, "Gd": 7, "Tb": 6, "Dy": 5,
       "Ho": 4, "Er": 3, "Tm": 2, "Yb": 1, "Lu": 0, "Ac": 0, "Th": 1, "Pa": 2, "U": 3,
       "Np": 4, "Pu": 5, "Am": 6, "Cm": 7, "Cf": 5}


def gfn2_best_seeds():
    """idx -> (source, group, seed, E_min, n_converged, xyz path) for the lowest-energy seed."""
    best = {}
    for f in glob.glob(str(XTB / "results*/G*/summary.csv")):
        if "validation" in f:
            continue
        src = Path(f).parts[-3]
        for r in pd.read_csv(f).itertuples():
            if pd.isna(r.best_seed):
                continue
            xyz = XTB / src / r.group / str(r.idx) / f"seed_{int(r.best_seed):02d}_opt.xyz"
            if not xyz.exists():
                continue
            if r.idx not in best or r.E_min_Eh < best[r.idx]["E_min_Eh"]:
                best[r.idx] = {"source": src, "group": r.group, "seed": int(r.best_seed),
                               "E_min_Eh": float(r.E_min_Eh), "n_converged": int(r.n_converged),
                               "xyz": xyz}
    return best


def relabel(src_xyz, dst_xyz, element):
    lines = src_xyz.read_text().splitlines()
    n = int(lines[0].split()[0])
    body, seen = [], 0
    for l in lines[2:2 + n]:
        if l.split()[0] == "Eu":
            seen += 1
            l = f"{element:<2s}" + l[2:]   # "Eu" is 2 chars; pad 1-letter symbols to keep columns
        body.append(l)
    assert seen == 1, f"{src_xyz}: expected exactly one Eu, found {seen}"
    comment = lines[1] + ("" if element == "Eu" else f" | Eu->{element}")
    dst_xyz.write_text("\n".join([lines[0], comment] + body) + "\n")


def main():
    (HERE / "starts").mkdir(parents=True, exist_ok=True)
    rows = []
    for name in ("trainVal", "test"):
        d = pd.read_csv(ORIG / f"{name}_dataset.csv", usecols=["SMILES", "Atomic_number", "Oxidation_state"])
        d = d[d.Oxidation_state == 3]
        rows.append(d.assign(split=name))
    df = pd.concat(rows, ignore_index=True)
    combos = (df.groupby(["SMILES", "Atomic_number"], as_index=False)
                .agg(n_rows=("split", "size"), splits=("split", lambda s: "+".join(sorted(set(s))))))
    print(f"trivalent rows {len(df)}, unique ligands {df.SMILES.nunique()}, "
          f"unique (SMILES, Z) combinations {len(combos)}")

    dt = pd.read_csv(XTB / "results" / "donor_table.csv")
    smi2idx = dict(zip(dt.SMILES, dt.idx))
    best = gfn2_best_seeds()

    out, bad = [], []
    for c in combos.itertuples():
        el = SYMBOL.get(c.Atomic_number)
        idx = smi2idx.get(c.SMILES)
        if el is None or idx is None or idx not in best:
            bad.append({"SMILES": c.SMILES, "Atomic_number": c.Atomic_number, "n_rows": c.n_rows,
                        "reason": "unknown atomic number" if el is None else
                                  ("SMILES not in donor_table" if idx is None else
                                   "no converged GFN2 seed for this ligand")})
            continue
        b, d = best[idx], dt[dt.idx == idx].iloc[0]
        name = f"{idx}_{el}"
        relabel(b["xyz"], HERE / "starts" / f"{name}.xyz", el)
        out.append({"name": name, "idx": idx, "element": el, "Z": c.Atomic_number, "uhf": UHF[el],
                    "chrg": int(d.total_charge),   # complex charge; not +3 for the 8 charged ligands
                    "SMILES": c.SMILES, "family": d.family, "group": d.group, "splits": c.splits,
                    "n_rows": c.n_rows, "denticity": int(d.denticity), "n_water": int(d.n_water),
                    "n_heavy": int(d.n_heavy), "n_atoms": int(d.n_atoms_complex),
                    "donors": d.donor_atoms, "gfn2_source": b["source"], "gfn2_seed": b["seed"],
                    "gfn2_E_min_Eh": b["E_min_Eh"], "gfn2_n_converged": b["n_converged"],
                    "gfn2_xyz": str(b["xyz"].relative_to(HERE.parent))})
    lst = pd.DataFrame(out).sort_values(["idx", "element"])
    lst.to_csv(HERE / "list.csv", index=False)
    # one task per line for run_gxtb.sh: name uhf chrg
    (HERE / "tasks.txt").write_text("".join(f"{r.name} {r.uhf} {r.chrg}\n" for r in lst.itertuples()))
    pd.DataFrame(bad).to_csv(HERE / "unmatched.csv", index=False)

    print(f"combinations with a start structure: {len(lst)}  (ligands {lst.idx.nunique()})")
    print("metal breakdown:")
    print(lst.groupby(["element", "uhf"]).size().rename("combinations").reset_index().to_string(index=False))
    print(f"unmatched combinations: {len(bad)}")
    if bad:
        print(pd.DataFrame(bad).groupby("reason").agg(combinations=("SMILES", "size"),
                                                     ligands=("SMILES", "nunique")).to_string())
        print(pd.DataFrame(bad).head(10).to_string())
    print(f"total atoms to optimise: {lst.n_atoms.sum()}, largest {lst.n_atoms.max()}")
    print("complex charge:", lst.chrg.value_counts().to_dict())


if __name__ == "__main__":
    main()
