"""Steps 1-2: list the Oxidation_state == 4 combinations and build ETKDG conformers for each.

conda activate env-gxtb && python xtb-complex-gxtb-ox4/build_list.py [--seeds 1 20]

Complex building is reused from xtb-complex-1/src/build.py (same donor table, same ETKDGv3 settings,
same water placement), with the metal swapped for the +4 element: build_complex() hardcodes Eu, so
BUILD.EU_Z is patched per metal and the metal's formal charge is set to +4 afterwards.
xtb-complex-1/ is only read (no bytecode written there).

Outputs: list.csv (one row per combination), tasks.txt (name uhf chrg per seed), starts/<idx>_<El>_s<seed>.xyz,
         build_log.csv (per combination: seeds embedded, method used).
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "orig-data").is_dir())   # project root
XTB = HERE.parent / "xtb-complex-1"
ORIG = ROOT / "orig-data" / "0.1-1"
sys.dont_write_bytecode = True
sys.path.insert(0, str(XTB / "src"))
import build as BUILD  # noqa: E402

SYMBOL = {89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm"}
# unpaired f electrons of the +4 ion
UHF4 = {"Th": 0, "Pa": 1, "U": 2, "Np": 3, "Pu": 4, "Am": 5, "Cm": 6}
Z_OF = {v: k for k, v in SYMBOL.items()}


def combinations():
    d = pd.concat([pd.read_csv(ORIG / f"{n}_dataset.csv",
                               usecols=["SMILES", "Atomic_number", "Oxidation_state"]).assign(split=n)
                   for n in ("trainVal", "test")], ignore_index=True)
    d = d[d.Oxidation_state == 4]
    return (d.groupby(["SMILES", "Atomic_number"], as_index=False)
             .agg(n_rows=("split", "size"), splits=("split", lambda s: "+".join(sorted(set(s))))))


def build_seeds(row, seeds):
    """(mol, method) per seed; ETKDG first, geometric water shell as fallback (as in xtb-complex-1)."""
    BUILD.EU_Z = Z_OF[row["element"]]
    donors = [int(x) for x in row["donors"].split()]
    out = {}
    for seed in seeds:
        mol, m_idx, _, _ = BUILD.build_complex(row["build_smiles"], donors, row["n_water"])
        if BUILD.embed(mol, seed) == 0:
            mol.GetAtomWithIdx(m_idx).SetFormalCharge(4)
            out[seed] = (mol, "etkdg")
            continue
        got = BUILD.build_complex_geom(row["build_smiles"], donors, row["n_water"], seed)
        if got is not None:
            mol = got[0]
            mol.GetAtomWithIdx(got[1]).SetFormalCharge(4)
            out[seed] = (mol, "geom")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", nargs=2, type=int, default=[1, 20], metavar=("FIRST", "LAST"))
    ap.add_argument("--only", nargs="*", default=None, help="restrict to these combination names")
    args = ap.parse_args()
    seeds = range(args.seeds[0], args.seeds[1] + 1)
    (HERE / "starts").mkdir(parents=True, exist_ok=True)

    combos = combinations()
    dt = pd.read_csv(XTB / "results" / "donor_table.csv").set_index("SMILES")
    rows, log = [], []
    for c in combos.itertuples():
        el = SYMBOL[c.Atomic_number]
        d = dt.loc[c.SMILES]
        rows.append({"name": f"{d.idx}_{el}", "idx": int(d.idx), "element": el, "Z": c.Atomic_number,
                     "uhf": UHF4[el], "chrg": int(d.ligand_charge) + 4, "SMILES": c.SMILES,
                     "build_smiles": d.build_smiles, "family": d.family, "group": d.group,
                     "splits": c.splits, "n_rows": c.n_rows, "denticity": int(d.denticity),
                     "n_water": int(d.n_water), "n_heavy": int(d.n_heavy),
                     "n_atoms": int(d.n_atoms_complex), "donors": d.donor_atoms,
                     "ligand_charge": int(d.ligand_charge)})
    lst = pd.DataFrame(rows).sort_values(["idx", "element"])
    if args.only:
        lst = lst[lst.name.isin(args.only)]
    print(f"ox4 combinations {len(combos)} (ligands {combos.SMILES.nunique()}), building {len(lst)}")

    tasks = []
    for r in lst.to_dict("records"):
        built = build_seeds(r, seeds)
        for seed, (mol, method) in built.items():
            name = f"{r['name']}_s{seed:02d}"
            BUILD.write_xyz(mol, HERE / "starts" / f"{name}.xyz",
                            f"{r['name']} seed {seed} {method} chrg={r['chrg']} uhf={r['uhf']}")
            tasks.append(f"{name} {r['uhf']} {r['chrg']}")
        log.append({"name": r["name"], "element": r["element"], "n_atoms": r["n_atoms"],
                    "seeds_requested": len(list(seeds)), "seeds_built": len(built),
                    "n_geom": sum(1 for _, m in built.values() if m == "geom")})
        print(f"  {r['name']:>12s} {r['n_atoms']:>4d} atoms  seeds {len(built)}/{len(list(seeds))}"
              + ("  (geometric water shell used)" if log[-1]["n_geom"] else ""))

    if args.only or args.seeds != [1, 20]:  # extension run: append
        lst_old = pd.read_csv(HERE / "list.csv")
        lst = pd.concat([lst_old, lst[~lst.name.isin(lst_old.name)]], ignore_index=True)
        with open(HERE / "tasks.txt", "a") as f:
            f.write("\n".join(tasks) + "\n")
        pd.concat([pd.read_csv(HERE / "build_log.csv"), pd.DataFrame(log)]).to_csv(HERE / "build_log.csv", index=False)
    else:
        (HERE / "tasks.txt").write_text("\n".join(tasks) + "\n")
        pd.DataFrame(log).to_csv(HERE / "build_log.csv", index=False)
    lst.to_csv(HERE / "list.csv", index=False)

    b = pd.DataFrame(log)
    print(f"\ncombinations {len(b)}; start structures {len(tasks)}; "
          f"seeds that needed the geometric water shell: {int(b.n_geom.sum())}")
    print("metals:", lst.groupby(["element", "uhf"]).size().to_dict())
    print("complex charge:", lst.chrg.value_counts().to_dict())


if __name__ == "__main__":
    main()
