"""Step 3-4: persistence images for every structure, then attach them to the experimental rows.

conda activate env-pd && python complex-dataset-pd/build_pi.py

Kernel settings are taken from xtb-complex-1/pi/config.json and left unchanged:
  H0  sigma = 0.10 * (1 + 0.5 * |EN_a - EN_b|) per MST edge, weight 1, infinite bar skipped
  H1  sigma = 0.10, weight = min(persistence / 0.5, 1)
  pixel-integrated Gaussians truncated at 3 sigma, no per-image normalisation
The electronegativity table is extended from the organic elements to all 23 metals (Pauling).

Grids are fitted ONCE over the PDs of all structures together, per subset and homology dimension
(never per metal or per split), and then applied to every structure:
  H0  every birth is exactly 0 (Vietoris-Rips), so the image is 1-D over persistence, 40 bins
  H1  2-D birth x persistence, 20 x 20

Outputs: pi_config.json, pi_features.csv, rows_with_pi.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "orig-data").is_dir())   # project root
DATASET = HERE.parent / "complex-dataset"
ORIG = ROOT / "orig-data" / "0.1-1"
CFG_SRC = json.load(open(ROOT / "xtb-prod" / "xtb-complex-1" / "pi" / "config.json"))

SUBSETS = ["heavy", "MNO", "MNOS"]
RES_H0, RES_H1 = 40, 20
TRUNC = CFG_SRC["kernel_truncation_sigma"]
H0C, H1C = CFG_SRC["H0"], CFG_SRC["H1"]
EN = dict(CFG_SRC["electronegativity_pauling"])          # H, C, N, O, P, S, Cl, Eu
EN.update({"La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14, "Pm": 1.13, "Sm": 1.17, "Gd": 1.20,
           "Tb": 1.10, "Dy": 1.22, "Ho": 1.23, "Er": 1.24, "Tm": 1.25, "Yb": 1.10, "Lu": 1.27,
           "Ac": 1.10, "Th": 1.30, "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28, "Am": 1.13,
           "Cm": 1.28, "Bk": 1.30, "Cf": 1.30})


def points(z, s):
    """H0 (finite bars: birth 0, persistence = MST edge, sigma from the EN difference) and H1."""
    death = z[f"{s}_H0_pair_death"].astype(float)
    el = z[f"{s}_H0_pair_elements"]
    H0 = z[f"{s}_H0"]
    assert np.allclose(np.sort(death), np.sort(H0[np.isfinite(H0[:, 1]), 1]), atol=1e-4)
    den = np.array([abs(EN[str(a)] - EN[str(b)]) for a, b in el]) if len(el) else np.zeros(0)
    h0 = {"b": np.zeros(len(death)), "p": death,
          "sigma": H0C["sigma_base"] * (1.0 + H0C["sigma_en_factor"] * den),
          "w": np.full(len(death), float(H0C["weight"]))}
    H1 = z[f"{s}_H1"]
    pers = (H1[:, 1] - H1[:, 0]) if len(H1) else np.zeros(0)
    h1 = {"b": H1[:, 0] if len(H1) else np.zeros(0), "p": pers,
          "sigma": np.full(len(pers), float(H1C["sigma"])),
          "w": np.minimum(pers / H1C["weight_ramp_end"], 1.0)}
    return {"H0": h0, "H1": h1}


def axis_mass(edges, centre, sigma):
    z = np.clip((edges[None, :] - centre[:, None]) / sigma[:, None], -TRUNC, TRUNC)
    return np.diff(ndtr(z), axis=1)


def render(pts, py, bx=None):
    """1-D (persistence only) if bx is None, else 2-D persistence x birth, flattened row-major."""
    n_p = len(py) - 1
    if bx is None:
        if len(pts["p"]) == 0:
            return np.zeros(n_p)
        return (pts["w"][:, None] * axis_mass(py, pts["p"], pts["sigma"])).sum(0)
    n_b = len(bx) - 1
    if len(pts["p"]) == 0:
        return np.zeros(n_p * n_b)
    Y = axis_mass(py, pts["p"], pts["sigma"])
    X = axis_mass(bx, pts["b"], pts["sigma"])
    return np.einsum("n,ni,nj->ij", pts["w"], Y, X).ravel()


def main():
    ids = sorted(p.stem for p in (HERE / "pd").glob("*.npz"))
    print(f"{len(ids)} diagrams")
    allpts = {}
    for sid in ids:
        z = np.load(HERE / "pd" / f"{sid}.npz", allow_pickle=True)
        allpts[sid] = {s: points(z, s) for s in SUBSETS}

    # ---- one grid per (subset, dim), fitted on every structure at once
    grids, cfg = {}, {}
    for s in SUBSETS:
        for h in ("H0", "H1"):
            p = np.concatenate([allpts[i][s][h]["p"] for i in ids])
            if p.size == 0:                       # subset with no bars at all
                p = np.zeros(1)
            b = np.concatenate([allpts[i][s][h]["b"] for i in ids])
            sig = np.concatenate([allpts[i][s][h]["sigma"] for i in ids])
            pad = TRUNC * (sig.max() if len(sig) else H1C["sigma"])
            p_range = (max(0.0, float(p.min()) - pad), float(p.max()) + pad)
            if h == "H0":
                py = np.linspace(*p_range, RES_H0 + 1)
                grids[(s, h)] = (py, None)
                cfg[f"{s}_{h}"] = {"persistence_range": list(p_range), "resolution": RES_H0,
                                   "birth": "all births are exactly 0, so the image is 1-D",
                                   "n_points": int(len(p))}
            else:
                b_range = (float(b.min()) - pad, float(b.max()) + pad) if len(b) else (0.0, 1.0)
                py, bx = np.linspace(*p_range, RES_H1 + 1), np.linspace(*b_range, RES_H1 + 1)
                grids[(s, h)] = (py, bx)
                cfg[f"{s}_{h}"] = {"persistence_range": list(p_range), "birth_range": list(b_range),
                                   "resolution": RES_H1, "n_points": int(len(b))}
    json.dump({"source_config": "xtb-complex-1/pi/config.json", "kernel": {"H0": H0C, "H1": H1C,
               "truncation_sigma": TRUNC, "normalize_per_image": False},
               "electronegativity_pauling": EN, "grids": cfg},
              open(HERE / "pi_config.json", "w"), indent=1)

    # ---- render
    cols, mats = [], []
    for s in SUBSETS:
        for h in ("H0", "H1"):
            py, bx = grids[(s, h)]
            M = np.stack([render(allpts[i][s][h], py, bx) for i in ids])
            n = M.shape[1]
            cols += [f"PI_{s}_{h}_p{k:03d}" for k in range(n)]
            mats.append(M)
            print(f"{s} {h}: {n} pixels, {100 * (M == 0).all(0).mean():.1f}% always zero, "
                  f"mass captured {M.sum(1).mean():.3f} per structure")
    values = np.hstack(mats)
    # 6 significant digits: keeps the CSV small without touching any original column later on
    mag = np.where(values == 0, 0.0, np.floor(np.log10(np.abs(np.where(values == 0, 1.0, values)))))
    factor = 10.0 ** (5 - mag)
    values = np.round(values * factor) / factor
    pi = pd.DataFrame(values, columns=cols)
    pi.insert(0, "structure_id", ids)
    keep = [c for c in cols if (pi[c] != 0).any()]
    dropped = len(cols) - len(keep)
    pi = pi[["structure_id"] + keep]
    pi.to_csv(HERE / "pi_features.csv", index=False)
    print(f"\npi_features.csv: {len(pi)} structures x {len(keep)} PI columns "
          f"({dropped} all-zero pixels dropped of {len(cols)})")

    # ---- attach to the experimental rows, rebuilt from source so every original column is present
    index = pd.read_csv(DATASET / "index.csv")
    key = index[["structure_id", "SMILES", "Z", "ox"]].rename(
        columns={"Z": "Atomic_number", "ox": "Oxidation_state"})
    parts = []
    for split in ("trainVal", "test"):
        part = pd.read_csv(ORIG / f"{split}_dataset.csv")
        part["split"] = split
        parts.append(part.copy())
    raw = pd.concat(parts, ignore_index=True)
    raw = raw[raw.Oxidation_state.isin([3, 4, 5, 6])]
    merged = raw.merge(key, on=["SMILES", "Atomic_number", "Oxidation_state"], how="inner")

    ref = pd.read_csv(DATASET / "rows.csv", usecols=["structure_id", "split", "SMILES", "Class_index"])
    assert len(merged) == len(ref), f"row count {len(merged)} != rows.csv {len(ref)}"
    assert merged.structure_id.value_counts().equals(ref.structure_id.value_counts()), \
        "structure_id distribution differs from rows.csv"
    assert merged.split.value_counts().equals(ref.split.value_counts()), "split distribution differs"

    out = merged.merge(pi, on="structure_id", how="left")
    assert out[keep].notna().all().all(), "some rows have no PI vector"
    front = ["structure_id", "split"]
    out = out[front + [c for c in out.columns if c not in front]]
    out.to_csv(HERE / "rows_with_pi.csv", index=False)   # no float_format: original columns stay exact
    n_fp = sum(c.startswith("CircularFP_") for c in out.columns)
    print(f"rows_with_pi.csv: {len(out)} rows x {out.shape[1]} columns "
          f"({n_fp} CircularFP, {len(keep)} PI, {out.shape[1] - n_fp - len(keep)} other)")
    print("split:", out.split.value_counts().to_dict(), "| structures used:", out.structure_id.nunique())


if __name__ == "__main__":
    main()
