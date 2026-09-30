"""Build the six feature-set CSVs for the "do 3D descriptors help?" experiment.

conda activate env-pd && python dataset/ml-featuresets-3d-vs-2d/build_featuresets.py

Source (read-only):
  dataset/complex-dataset-pd/rows_with_pi.csv   7053 g-xTB rows, original columns + PI
  dataset/complex-dataset-pd/pd/*.npz           persistence diagrams (PD summary is computed here)
  dataset/complex-dataset/index.csv             denticity per structure

Blocks
  2D        CircularFP_* fingerprint columns
  rdkit     RDKit descriptors (MaxAbsEStateIndex ... fr_unbrch_alkane)
  cond      experimental conditions (extractant/acid/metal concentration, temperature, diluent A/B)
  metal     metal properties + Oxidation_state (3 or 4 here, so it varies and is a real feature)
  dent      denticity
  PDsummary topology summary features, computed with xtb-complex-1/pd/compute_pd.py subset_features()
            for the heavy / MNO / MNOS subsets (metal columns renamed Eu -> M)
  PI        persistence-image pixels

rdkit + cond + metal are in every set, so the sets differ only by 2D / dent / PDsummary / PI.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "orig-data").is_dir())
PDDIR = HERE.parent / "complex-dataset-pd"
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "xtb-prod" / "xtb-complex-1" / "pd"))
import compute_pd as cpd  # noqa: E402  (subset_features, BINS)

SUBSETS = {"heavy": None, "MNO": "EuNO", "MNOS": "EuNOS"}   # new name -> name used by subset_features
META = ["structure_id", "split", "SMILES", "SMILES_class", "element", "Oxidation_state_meta",
        "idx", "family", "Class_index", "Distribution_ratio", "Log_D"]
SETS = {
    "fs1_2D_rdkit_cond":            ["2D"],
    "fs2_2D_rdkit_cond_denticity":  ["2D", "dent"],
    "fs3_PDsummary_rdkit_cond":     ["PDsummary"],
    "fs4_PI_rdkit_cond":            ["PI"],
    "fs5_2D_PDsummary_rdkit_cond":  ["2D", "PDsummary"],
    "fs6_2D_PI_rdkit_cond":         ["2D", "PI"],
}
BASE_BLOCKS = ["rdkit", "cond", "metal"]


def pd_summary():
    """One row per structure: subset_features() for each subset, metal columns renamed Eu -> M."""
    rows = []
    for f in sorted((PDDIR / "pd").glob("*.npz")):
        z = np.load(f, allow_pickle=True)
        r = {"structure_id": f.stem}
        for s, old in SUBSETS.items():
            local = (None if old is None
                     else int(np.where(z[f"{s}_atom_index"] == z["metal_index"])[0][0]))
            feat = cpd.subset_features(old or s, z[f"{s}_elements"], z[f"{s}_H0_pair_death"],
                                       z[f"{s}_H1"], z[f"{s}_H0_pair_local"], local)
            for k, v in feat.items():
                k = k.replace("EuNOS", "MNOS").replace("EuNO", "MNO")
                k = k.replace("_Eu_mst", "_M_mst").replace("H0_coord_Eu", "H0_coord_M")
                r[k] = v
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    df = pd.read_csv(PDDIR / "rows_with_pi.csv")
    print(f"rows {len(df)}, columns {df.shape[1]}")

    cols = list(df.columns)
    i = cols.index
    blocks = {
        "2D": [c for c in cols if c.startswith("CircularFP_")],
        "rdkit": cols[i("MaxAbsEStateIndex"): i("fr_unbrch_alkane") + 1],
        "cond": cols[i("Extractant_conc_M"): i("Log_metal_conc_mM") + 1],
        "metal": cols[i("Atomic_number"): i("Ionic_radius_nm") + 1] + ["Oxidation_state"],
        "PI": [c for c in cols if c.startswith("PI_")],
    }

    # denticity, and the element / idx used only as metadata
    index = pd.read_csv(HERE.parent / "complex-dataset" / "index.csv",
                        usecols=["structure_id", "idx", "element", "denticity"])
    df = df.merge(index, on="structure_id", how="left", validate="m:1")
    assert df.denticity.notna().all(), "some rows have no denticity"
    blocks["dent"] = ["denticity"]

    # PD summary
    summ = pd_summary()
    df = df.merge(summ, on="structure_id", how="left", validate="m:1")
    blocks["PDsummary"] = [c for c in summ.columns if c != "structure_id"]
    assert df[blocks["PDsummary"]].notna().all().all(), "some rows have no PD summary"

    extra_meta = {"Oxidation_state_meta": df.Oxidation_state}   # readable copy; the feature stays in `metal`
    if "family" not in df:
        extra_meta["family"] = ""
    df = pd.concat([df, pd.DataFrame(extra_meta, index=df.index)], axis=1).copy()

    json.dump({"rows": len(df), "target": "Class_index",
               "note": "rdkit + cond + metal are included in every set; "
                       "Distribution_ratio and Log_D are metadata, never features",
               "metadata_columns": META,
               "blocks": {k: v for k, v in blocks.items()},
               "block_sizes": {k: len(v) for k, v in blocks.items()},
               "sets": {name: BASE_BLOCKS + extra for name, extra in SETS.items()}},
              open(HERE / "feature_blocks.json", "w"), indent=1)

    summary = []
    for name, extra in SETS.items():
        feats = [c for b in BASE_BLOCKS + extra for c in blocks[b]]
        assert len(feats) == len(set(feats))
        out = df[META + feats]
        assert not out[feats].isna().any().any(), f"{name}: NaN in features"
        path = HERE / f"{name}.csv"
        out.to_csv(path, index=False)
        summary.append({"file": path.name, "blocks": "+".join(BASE_BLOCKS + extra),
                        "n_features": len(feats), "n_rows": len(out),
                        "size_MB": round(path.stat().st_size / 1048576, 1)})
        print(f"{path.name:32s} {len(feats):5d} features  {summary[-1]['size_MB']:6.1f} MB")

    s = pd.DataFrame(summary)
    s.to_csv(HERE / "featureset_summary.csv", index=False)

    # unique ligand list for train/val splitting, cut down to the ligands present here
    uniq = pd.read_csv(ROOT / "orig-data" / "unique_smiles.csv")
    per = df.groupby("SMILES").agg(n_rows=("split", "size"), n_structures=("structure_id", "nunique"),
                                   n_metals=("Atomic_number", "nunique"),
                                   split_in_data=("split", "first")).reset_index()
    assert df.groupby("SMILES").split.nunique().eq(1).all(), "a ligand appears in both splits"
    out = uniq[uniq.SMILES.isin(per.SMILES)].merge(per, on="SMILES", validate="1:1")
    assert len(out) == df.SMILES.nunique(), "some ligand is missing from orig-data/unique_smiles.csv"
    assert (out.split == out.split_in_data).all(), "split label disagrees with the rows"
    out = out.drop(columns=["split_in_data"])
    out.to_csv(HERE / "unique_smiles.csv", index=False)
    print(f"\nunique_smiles.csv: {len(out)} ligands "
          f"({out.split.value_counts().to_dict()}), SMILES_class {out.SMILES_class.nunique()} kinds, "
          f"smallest class {out.SMILES_class.value_counts().min()} ligands")
    print("\nblock sizes:", {k: len(v) for k, v in blocks.items()})
    print("split:", df.split.value_counts().to_dict())
    print("class:", df.Class_index.value_counts().sort_index().to_dict())
    print("oxidation state:", df.Oxidation_state.value_counts().sort_index().to_dict())
    print("ligands:", df.SMILES.nunique(), "| structures:", df.structure_id.nunique())


if __name__ == "__main__":
    main()
