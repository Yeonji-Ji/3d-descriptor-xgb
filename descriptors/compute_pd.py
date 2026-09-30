"""Step 1-2: persistence diagrams for every structure in complex-dataset/, for any metal.

conda activate env-pd
python complex-dataset-pd/compute_pd.py --check     # Eu equivalence check only
python complex-dataset-pd/compute_pd.py             # check, then all structures

The maths and the numeric settings come from xtb-complex-1/pd/compute_pd.py (imported, not copied):
Vietoris-Rips with ripser maxdim=1 on the pairwise distance matrix, H0 deaths cross-checked against
the scipy MST edges. Only the subset definition is generalised: the metal is detected per structure
instead of being hard-coded to Eu, so
    heavy  all non-H atoms          (unchanged)
    MNO    metal + N + O            (was EuNO)
    MNOS   metal + N + O + S        (was EuNOS)
Everything else - cutoffs, bins, maxdim - is untouched.

Outputs: pd/<structure_id>.npz (one per structure), pd_log.csv, eu_check.csv
"""
import argparse
import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "orig-data").is_dir())   # project root
DATASET = HERE.parent / "complex-dataset"
XTB = ROOT / "xtb-prod" / "xtb-complex-1"
sys.dont_write_bytecode = True
sys.path.insert(0, str(XTB / "pd"))
import compute_pd as cpd  # noqa: E402  (diagram(), BINS, NO_CN, S_EDGE)

ORGANIC = {"H", "C", "N", "O", "S", "P", "Cl", "F", "Br", "I", "Si", "B"}
SUBSETS = {"heavy": None, "MNO": ("N", "O"), "MNOS": ("N", "O", "S")}   # metal added per structure
# the three Eu structures used for the equivalence check against the original code
EU_CHECK = [("G3", 40, 11), ("G4", 187, 15), ("G7", 34, 40)]
OLD_NAME = {"heavy": "heavy", "MNO": "EuNO", "MNOS": "EuNOS"}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    el = np.array([l.split()[0] for l in lines[2:2 + n]])
    X = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
    assert len(el) == n, f"{path}: header says {n} atoms, found {len(el)}"
    return el, X


def find_metal(el):
    """The one element that is not part of the organic ligand / water."""
    m = sorted(set(el) - ORGANIC)
    assert len(m) == 1, f"expected exactly one metal, found {m}"
    return m[0]


def diagrams(el, X):
    """PDs per subset, plus the MST pair tables the PI step needs."""
    metal = find_metal(el)
    keep = np.where(el != "H")[0]
    out = {"metal": np.str_(metal), "metal_index": int(np.where(el == metal)[0][0]),
           "n_atoms": len(el)}
    for s, allowed in SUBSETS.items():
        idx = keep if allowed is None else keep[np.isin(el[keep], (metal,) + allowed)]
        H0, H1, edges, lens = cpd.diagram(X[idx])
        assert np.allclose(H0[:, 0], 0.0), "H0 births are not all zero"
        out.update({f"{s}_atom_index": idx, f"{s}_elements": el[idx], f"{s}_coords": X[idx],
                    f"{s}_H0": H0, f"{s}_H1": H1,
                    f"{s}_H0_pair_local": edges, f"{s}_H0_pair_atom_index": idx[edges],
                    f"{s}_H0_pair_elements": el[idx][edges], f"{s}_H0_pair_death": lens})
    return out


def one(path, out_dir):
    sid = Path(path).stem
    try:
        el, X = read_xyz(path)
        z = diagrams(el, X)
        np.savez_compressed(out_dir / f"{sid}.npz", **z)
        return {"structure_id": sid, "metal": str(z["metal"]), "n_atoms": z["n_atoms"],
                **{f"{s}_n": len(z[f"{s}_atom_index"]) for s in SUBSETS},
                **{f"{s}_H1_n": len(z[f"{s}_H1"]) for s in SUBSETS}, "error": ""}
    except Exception as e:
        return {"structure_id": sid, "error": f"{type(e).__name__}: {e}"}


def eu_check():
    """Same PDs as the original Eu-only code? Compare against the stored pd/results npz."""
    rows = []
    for group, idx, seed in EU_CHECK:
        xyz = XTB / "results" / group / str(idx) / f"seed_{seed:02d}_opt.xyz"
        old = np.load(XTB / "pd" / "results" / "raw" / group / f"idx{idx}_seed{seed:02d}.npz",
                      allow_pickle=True)
        new = diagrams(*read_xyz(xyz))
        r = {"idx": idx, "group": group, "seed": seed, "metal": str(new["metal"])}
        for s in SUBSETS:
            o = OLD_NAME[s]
            r[f"{s}_atoms_equal"] = bool(np.array_equal(new[f"{s}_atom_index"], old[f"{o}_atom_index"]))
            for h in ("H0", "H1"):
                a, b = new[f"{s}_{h}"], old[f"{o}_{h}"]
                r[f"{s}_{h}_shape_equal"] = a.shape == b.shape
                if a.shape != b.shape:
                    r[f"{s}_{h}_max_abs_diff"] = np.nan
                elif a.size == 0:                       # no loops in this subset
                    r[f"{s}_{h}_max_abs_diff"] = 0.0
                else:                                   # the infinite H0 bar compares as 0 vs 0
                    r[f"{s}_{h}_max_abs_diff"] = float(np.max(np.abs(np.where(np.isinf(a), 0.0, a)
                                                                    - np.where(np.isinf(b), 0.0, b))))
            r[f"{s}_pairs_equal"] = bool(np.array_equal(new[f"{s}_H0_pair_elements"],
                                                        old[f"{o}_H0_pair_elements"]))
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "eu_check.csv", index=False)
    cols = [c for c in df.columns if c.endswith("_equal")]
    diffs = [c for c in df.columns if c.endswith("_max_abs_diff")]
    ok = bool(df[cols].all().all()) and float(df[diffs].max().max()) < 1e-9
    print("Eu equivalence check (generalised code vs the original Eu-only PDs):")
    print(df[["idx", "group", "seed", "metal"] + cols].to_string(index=False))
    print(f"largest absolute difference in any diagram: {df[diffs].max().max():.3e}")
    print("identical\n" if ok else "DIFFERENT\n")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="run only the Eu equivalence check")
    ap.add_argument("--nproc", type=int, default=9)
    args = ap.parse_args()

    if not eu_check():
        sys.exit("generalised code does not reproduce the original Eu PDs - stopping")
    if args.check:
        return

    out_dir = HERE / "pd"
    out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted((DATASET / "structures").glob("*.xyz"))
    print(f"{len(files)} structures")
    with Pool(args.nproc) as p:
        rows = p.map(partial(one, out_dir=out_dir), files, chunksize=8)
    log = pd.DataFrame(rows)
    log.to_csv(HERE / "pd_log.csv", index=False)
    bad = log[log.error != ""]
    print(f"done: {len(log) - len(bad)} ok, {len(bad)} failed")
    if len(bad):
        print(bad[["structure_id", "error"]].to_string(index=False))
    print(log.groupby("metal").size().rename("structures").to_string())
    print(log[[f"{s}_n" for s in SUBSETS] + [f"{s}_H1_n" for s in SUBSETS]].describe().round(1).to_string())


if __name__ == "__main__":
    main()
