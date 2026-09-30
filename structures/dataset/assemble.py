"""Assemble the usable g-xTB structures into one dataset and pair them with the experimental rows.

conda activate env-pd && python complex-dataset/assemble.py

Sources (read-only): xtb-complex-gxtb (+3), xtb-complex-gxtb-ox4 (+4),
xtb-complex-gxtb-actinyl (+5/+6, dropped entirely), xtb-complex-1 (donor table).

Selection
  +3  clean combinations, plus flagged ones whose only problem is a lost water:
      CN fell and every design donor is still within the donor cutoff (n_donor_bound >= denticity).
      Dropped: a donor left, SCF failed, the optimisation did not converge, or CN rose
      (not a lost water, so outside the inclusion rule - listed separately as cn_gain).
  +4  all Th; Np and Pu only with metal spin excess <= 0.5 (the rest behave as reduced +3).
  +5/+6  dropped.

Outputs: structures/<idx>_<El><ox>.xyz, index.csv, rows.csv, dropped.csv
"""
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "orig-data").is_dir())   # project root
BASE = ROOT / "xtb-prod"                     # the g-xTB / GFN2 run folders
ORIG = ROOT / "orig-data" / "0.1-1"
SRC3, SRC4, SRC56 = "xtb-complex-gxtb", "xtb-complex-gxtb-ox4", "xtb-complex-gxtb-actinyl"
SPIN_MAX = 0.5


def select_ox3():
    d = pd.read_csv(BASE / SRC3 / "structures.csv")
    d["ox"], d["chrg"] = 3, 3
    f = d["flags"].fillna("")
    run_bad = f.str.contains("missing|scf_fail|not_converged")
    donors_ok = d.n_donor_bound >= d.denticity
    cn_down, cn_up = d.CN < d.start_CN, d.CN > d.start_CN

    d["reason"] = np.select(
        [run_bad,
         f.eq(""),
         cn_down & donors_ok,
         cn_up & ~donors_ok,
         cn_up & donors_ok,
         ~donors_ok],
        ["run failed (SCF / not converged / no output)",
         "clean",
         "water lost, all design donors still bound",
         "CN rose and a design donor left",
         "CN rose (not a lost water); outside the inclusion rule",
         "a design donor left"],
        default="other flag")
    d["keep"] = d.reason.isin(["clean", "water lost, all design donors still bound"])
    return d


def select_ox4():
    d = pd.read_csv(BASE / SRC4 / "structures.csv")
    d["ox"] = 4
    d["reason"] = np.where(d.element == "Th", "Th, kept by rule",
                  np.where(d.spin_excess <= SPIN_MAX, f"spin excess <= {SPIN_MAX}",
                           f"spin excess > {SPIN_MAX} (behaves as a reduced +3 complex)"))
    d["keep"] = (d.element == "Th") | (d.spin_excess <= SPIN_MAX)
    return d


def select_ox56():
    d = pd.read_csv(BASE / SRC56 / "structures.csv")
    d["keep"] = False
    d["reason"] = "oxidation state 5/6 dropped (bent O=M=O, short M=O; see that folder's README)"
    return d


def main():
    (HERE / "structures").mkdir(parents=True, exist_ok=True)
    d3, d4, d56 = select_ox3(), select_ox4(), select_ox56()

    # ---- structures + index
    rows = []
    for d, src, best in ((d3, SRC3, None), (d4, SRC4, "best")):
        for r in d[d.keep].to_dict("records"):
            sid = f"{r['idx']}_{r['element']}{r['ox']}"
            src_xyz = (BASE / src / "best" / f"{r['name']}.xyz" if best
                       else BASE / src / "runs" / r["name"] / "xtbopt.xyz")
            shutil.copy(src_xyz, HERE / "structures" / f"{sid}.xyz")
            rows.append({"structure_id": sid, "idx": r["idx"], "SMILES": r["SMILES"],
                         "element": r["element"], "Z": r["Z"], "ox": r["ox"],
                         "chrg": r.get("chrg", r["ox"]), "uhf": r["uhf"], "source": src,
                         "family": r["family"], "denticity": r["denticity"],
                         "n_water": r.get("n_water"), "CN": r["CN"],
                         "M_O_mean": r["M_O_mean"], "M_N_mean": r["M_N_mean"],
                         "M_q": r["M_q"], "M_spin": r["M_spin"],
                         "spin_excess": r.get("spin_excess", r["M_spin"] - r["uhf"]),
                         "n_donor_bound": r["n_donor_bound"], "n_water_bound": r["n_water_bound"],
                         "energy_Eh": r["energy_Eh"], "kept_reason": r["reason"],
                         "source_flags": r["flags"] if isinstance(r["flags"], str) else ""})
    index = pd.DataFrame(rows).sort_values(["idx", "element", "ox"])
    index.to_csv(HERE / "index.csv", index=False)

    # ---- experimental rows
    parts = []
    for n in ("trainVal", "test"):
        part = pd.read_csv(ORIG / f"{n}_dataset.csv")
        part["split"] = n
        parts.append(part.copy())            # copy(): these frames are ~1865 columns wide
    raw = pd.concat(parts, ignore_index=True)
    raw = raw[raw.Oxidation_state.isin([3, 4, 5, 6])]
    key = index[["structure_id", "SMILES", "Z", "ox"]].rename(columns={"Z": "Atomic_number", "ox": "Oxidation_state"})
    merged = raw.merge(key, on=["SMILES", "Atomic_number", "Oxidation_state"], how="left")
    keep_cols = [c for c in merged.columns if not c.startswith("CircularFP_")]
    covered = merged[merged.structure_id.notna()]
    covered[["structure_id", "split"] + [c for c in keep_cols if c not in ("structure_id", "split")]] \
        .to_csv(HERE / "rows.csv", index=False)

    # ---- dropped: combinations and rows
    drop_comb = pd.concat([
        d[~d.keep].assign(source=src, level="combination")[
            ["level", "source", "name", "idx", "SMILES", "element", "ox", "reason"]]
        for d, src in ((d3, SRC3), (d4, SRC4), (d56, SRC56))], ignore_index=True)
    lost = merged[merged.structure_id.isna()]
    lookup = pd.concat([d.assign(source=s) for d, s in ((d3, SRC3), (d4, SRC4), (d56, SRC56))])
    lookup = lookup.set_index(["SMILES", "Z", "ox"]).reason.to_dict()
    drop_rows = lost.groupby(["SMILES", "Atomic_number", "Oxidation_state", "split"], as_index=False).size()
    drop_rows["reason"] = [lookup.get((s, z, o), "no structure was ever built for this ligand")
                           for s, z, o in zip(drop_rows.SMILES, drop_rows.Atomic_number, drop_rows.Oxidation_state)]
    drop_rows = drop_rows.rename(columns={"Atomic_number": "Z", "Oxidation_state": "ox", "size": "n_rows"})
    drop_rows["level"], drop_rows["source"] = "rows", ""
    pd.concat([drop_comb, drop_rows], ignore_index=True).to_csv(HERE / "dropped.csv", index=False)

    # ---- report
    print("STRUCTURES")
    print(f"  +3: kept {int(d3.keep.sum())} of {len(d3)}")
    print(d3.groupby("reason").agg(n=("name", "size"), kept=("keep", "sum")).to_string())
    print(f"\n  +4: kept {int(d4.keep.sum())} of {len(d4)}")
    print(d4.groupby(["element", "reason"]).agg(n=("name", "size"), kept=("keep", "sum")).to_string())
    print(f"\n  +5/+6: dropped all {len(d56)}")
    print(f"\n  total structures written: {len(index)}")
    print(index.groupby(["ox", "element"]).size().rename("structures").to_string())

    print("\nROWS")
    tot = merged.groupby("Oxidation_state").size().rename("rows")
    cov = covered.groupby("Oxidation_state").size().rename("covered")
    t = pd.concat([tot, cov], axis=1).fillna(0).astype(int)
    t["lost"] = t.rows - t.covered
    t["covered_%"] = (100 * t.covered / t.rows).round(1)
    print(t.to_string())
    print("\nby metal (oxidation state 3 and 4 only):")
    m = merged[merged.Oxidation_state.isin([3, 4])].copy()
    m["el"] = m.Atomic_number.map(dict(zip(index.Z, index.element))
                                  | {61: "Pm", 98: "Cf", 96: "Cm", 93: "Np", 94: "Pu"})
    mm = pd.concat([m.groupby(["Oxidation_state", "el"]).size().rename("rows"),
                    m[m.structure_id.notna()].groupby(["Oxidation_state", "el"]).size().rename("covered")],
                   axis=1).fillna(0).astype(int)
    mm["lost"] = mm.rows - mm.covered
    print(mm[mm.lost > 0].to_string())
    print(f"\ntrain/test split of covered rows: {covered.split.value_counts().to_dict()}")
    print(f"covered rows {len(covered)} / {len(merged)}; ligands {covered.SMILES.nunique()}")


if __name__ == "__main__":
    main()
