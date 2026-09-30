#!/usr/bin/env python
"""Persistence diagrams + summary topology features for the lowest-energy xTB conformer
of each [Eu(L)(H2O)n] complex.

Run end to end (conda env env-pd), from anywhere:
    python pd/compute_pd.py

Inputs (read-only):  results/, results-waterfix/, results/donor_table.csv
Outputs:
    pd/results/   complexes whose best seed comes from results/
    pd/waterfix/  the 17 complexes recovered in results-waterfix/
  each containing features.csv, skipped.csv, sanity_check.txt,
  raw/<group>/idx<idx>_seed<XX>.npz, plots/idx<idx>.png (heavy | EuNO)
  and plots/idx<idx>_EuNO_vs_EuNOS.png for sulfur-binding examples
"""
import os, re, argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree
from ripser import ripser

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PD = os.path.join(BASE, "pd")
SOURCES = {"results": "results", "waterfix": "results-waterfix"}
BINS = {"cov": (0.0, 1.8), "coord": (1.8, 3.0), "outer": (3.0, 8.0)}
SUBSETS = {"heavy": None, "EuNO": ("Eu", "N", "O"), "EuNOS": ("Eu", "N", "O", "S")}
PLOT_IDX = {"results": [174, 12, 51, 107, 169], "waterfix": [262, 81, 122]}
# sulfur-binding examples, EuNO vs EuNOS side by side (all S-donor ligands are in results/)
PLOT_S = {"results": [148, 150, 216], "waterfix": []}
FEW_CONVERGED = 15
# Eu-S bonds run 2.79-3.27 A (next S donor 3.90 A, non-donor thiophene S of idx 142 at 3.38 A), so
# Eu-S edges and S atoms use their own 3.3 A edge in the *_elem columns (decision B1). N/O keep 3.0/3.2.
S_EDGE = 3.3
NO_CN = 3.2


# ---------------------------------------------------------------- io
def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    m = re.search(r"energy:\s*(-?\d+\.\d+)", lines[1])
    rows = [ln.split() for ln in lines[2:2 + n]]
    assert len(rows) == n, f"{path}: header says {n} atoms, found {len(rows)}"
    elems = np.array([r[0] for r in rows])
    xyz = np.array([[float(v) for v in r[1:4]] for r in rows])
    return elems, xyz, (float(m.group(1)) if m else None)


# ---------------------------------------------------------------- topology
def diagram(coords):
    """Vietoris-Rips PD (maxdim=1, Angstrom) plus the MST that realises the finite H0 bars."""
    D = squareform(pdist(coords))
    dg = ripser(D, maxdim=1, distance_matrix=True)["dgms"]
    H0, H1 = dg[0], dg[1]
    mst = minimum_spanning_tree(D).tocoo()
    order = np.argsort(mst.data, kind="stable")
    edges = np.stack([mst.row[order], mst.col[order]], axis=1).astype(int)
    lengths = mst.data[order]
    deaths = np.sort(H0[np.isfinite(H0[:, 1]), 1])
    assert len(deaths) == len(lengths) == len(coords) - 1, "H0 bar count != MST edge count"
    assert np.allclose(deaths, lengths, atol=1e-4), "ripser H0 deaths != scipy MST edge lengths"
    return H0, H1, edges, lengths


def count_bin(v, lo, hi):
    return int(((v >= lo) & (v < hi)).sum())


def subset_features(s, elems, H0lens, H1, edges, eu_local):
    f = {f"{s}_n_atoms": len(elems)}
    life1 = H1[:, 1] - H1[:, 0] if len(H1) else np.zeros(0)
    birth1 = H1[:, 0] if len(H1) else np.zeros(0)
    for b, (lo, hi) in BINS.items():
        f[f"{s}_H0_{b}"] = count_bin(H0lens, lo, hi)
        f[f"{s}_H1_{b}"] = count_bin(birth1, lo, hi)
    f[f"{s}_H0_beyond8"] = int((H0lens >= 8.0).sum())
    f[f"{s}_H1_beyond8"] = int((birth1 >= 8.0).sum())
    f[f"{s}_H0_totlife"] = float(H0lens.sum())
    f[f"{s}_H0_meanlife"] = float(H0lens.mean()) if len(H0lens) else 0.0
    f[f"{s}_H0_maxlife"] = float(H0lens.max()) if len(H0lens) else 0.0
    f[f"{s}_H1_n"] = int(len(H1))
    f[f"{s}_H1_totlife"] = float(life1.sum())
    f[f"{s}_H1_meanlife"] = float(life1.mean()) if len(life1) else 0.0
    f[f"{s}_H1_maxlife"] = float(life1.max()) if len(life1) else 0.0
    f[f"{s}_H1_per_atom"] = len(H1) / len(elems)
    for b in BINS:
        f[f"{s}_H1_frac_{b}"] = f[f"{s}_H1_{b}"] / len(H1) if len(H1) else 0.0
    assert sum(f[f"{s}_H0_{b}"] for b in BINS) + f[f"{s}_H0_beyond8"] == len(elems) - 1
    if eu_local is not None:
        touch = (edges[:, 0] == eu_local) | (edges[:, 1] == eu_local)
        lo, hi = BINS["coord"]
        f[f"{s}_Eu_mst_degree"] = int(touch.sum())
        f[f"{s}_H0_coord_Eu"] = int((touch & (H0lens >= lo) & (H0lens < hi)).sum())
        f[f"{s}_Eu_mst_minlen"] = float(H0lens[touch].min()) if touch.any() else np.nan
        if s == "EuNOS":
            # B1: an Eu-S MST edge is "coord" up to S_EDGE; every other edge keeps the global bins
            eu_s = touch & (elems[edges] == "S").any(axis=1)
            top = np.where(eu_s, S_EDGE, hi)
            coord = (H0lens >= lo) & (H0lens < top)
            outer = (H0lens >= top) & (H0lens < BINS["outer"][1])
            f["EuNOS_H0_coord_elem"] = int(coord.sum())
            f["EuNOS_H0_outer_elem"] = int(outer.sum())
            f["EuNOS_H0_coord_Eu_elem"] = int((touch & coord).sum())
            assert f["EuNOS_H0_cov"] + f["EuNOS_H0_coord_elem"] + f["EuNOS_H0_outer_elem"] \
                + f["EuNOS_H0_beyond8"] == len(elems) - 1
    return f


# ---------------------------------------------------------------- one complex
def process(task):
    try:
        elems, xyz, e_xyz = read_xyz(task["xyz"])
        L, nw = task["n_heavy"], task["n_water"]
        # atom order written by the build: ligand heavy atoms, Eu, water O, then all H
        assert (elems == "Eu").sum() == 1 and elems[L] == "Eu", f"Eu not at index {L}"
        wat = np.arange(L + 1, L + 1 + nw)
        assert np.all(elems[wat] == "O"), "water O indices do not hold O"
        don = np.array(task["donors"], dtype=int)
        assert list(elems[don]) == task["donor_elems"], \
            f"donor elements {list(elems[don])} != donor table {task['donor_elems']}"
        if e_xyz is not None:
            assert abs(e_xyz - task["energy_Eh"]) < 1e-6, "xyz energy != summary E_min_Eh"

        keep = np.where(elems != "H")[0]
        row = {k: task[k] for k in ("source", "group", "idx", "family", "SMILES", "split", "seed",
                                    "energy_Eh", "denticity", "n_water", "total_charge", "n_converged")}
        npz = dict(idx=task["idx"], seed=task["seed"], energy_Eh=task["energy_Eh"],
                   group=np.str_(task["group"]), source=np.str_(task["source"]),
                   eu_index=L, water_O_index=wat, donor_index=don)
        for s, allowed in SUBSETS.items():
            idx = keep if allowed is None else keep[np.isin(elems[keep], allowed)]
            H0, H1, edges, lens = diagram(xyz[idx])
            eu_local = int(np.where(idx == L)[0][0]) if s != "heavy" else None
            row.update(subset_features(s, elems[idx], lens, H1, edges, eu_local))
            npz.update({f"{s}_atom_index": idx, f"{s}_elements": elems[idx], f"{s}_coords": xyz[idx],
                        f"{s}_H0": H0, f"{s}_H1": H1,
                        f"{s}_H0_pair_local": edges, f"{s}_H0_pair_atom_index": idx[edges],
                        f"{s}_H0_pair_elements": elems[idx][edges], f"{s}_H0_pair_death": lens})

        # direct (non-topological) coordination geometry
        d = np.linalg.norm(xyz - xyz[L], axis=1)
        isNO, isS = np.isin(elems, ("N", "O")), elems == "S"
        row["CN_NO_3p0"] = int((isNO & (d < 3.0)).sum())
        row["CN_NO_3p2"] = int((isNO & (d < 3.2)).sum())
        row["CN_NOS_3p0"] = int(((isNO | isS) & (d < 3.0)).sum())
        row["CN_NOS_3p2"] = int(((isNO | isS) & (d < 3.2)).sum())
        row["CN_NOS_elem"] = int((isNO & (d < NO_CN)).sum() + (isS & (d < S_EDGE)).sum())
        row["CN_design"] = task["denticity"] + nw
        row["n_donor_within_3p2"] = int((d[don] < 3.2).sum())
        row["n_donor_within_elem"] = int((d[don] < np.where(elems[don] == "S", S_EDGE, NO_CN)).sum())
        row["n_water_within_3p2"] = int((d[wat] < 3.2).sum()) if nw else 0
        for name, sel in (("donor", don), ("water", wat)):
            row[f"d_Eu_{name}_min"] = float(d[sel].min()) if len(sel) else np.nan
            row[f"d_Eu_{name}_mean"] = float(d[sel].mean()) if len(sel) else np.nan
            row[f"d_Eu_{name}_max"] = float(d[sel].max()) if len(sel) else np.nan
        row["n_S_donor"] = task["donor_elems"].count("S")
        row["flag_S_donor"] = int(row["n_S_donor"] > 0)
        row["flag_few_converged"] = int(task["n_converged"] < FEW_CONVERGED)
        row["flag_donor_detached"] = bool(row["n_donor_within_elem"] < task["denticity"])
        row["flag_donor_detached_3p2"] = bool(row["n_donor_within_3p2"] < task["denticity"])

        os.makedirs(os.path.dirname(task["npz"]), exist_ok=True)
        np.savez_compressed(task["npz"], **npz)
        return row, None
    except Exception as e:
        return None, (task["idx"], task["group"], f"{type(e).__name__}: {e}")


# ---------------------------------------------------------------- report
def sanity(df, out_dir):
    out = []

    def agree(a, b, label, frame=None):
        f = df if frame is None else frame
        diff = f[a] - f[b]
        r = ("n/a (a column is constant)" if f[a].nunique() < 2 or f[b].nunique() < 2
             else f"{f[a].corr(f[b]):.3f}")
        out.append(f"{a} vs {b} ({label}, n={len(f)}): exact {100*(diff==0).mean():.1f}%, within +-1 "
                   f"{100*(diff.abs()<=1).mean():.1f}%, mean({a}-{b}) {diff.mean():+.2f}, "
                   f"MAE {diff.abs().mean():.2f}, Pearson r {r}")
        out.append("    difference distribution: " + ", ".join(
            f"{int(k):+d}:{v}" for k, v in diff.value_counts().sort_index().items()))

    def npz(r):
        return np.load(os.path.join(out_dir, "raw", r.group, f"idx{int(r.idx)}_seed{int(r.seed):02d}.npz"))

    S, noS = df[df.flag_S_donor == 1], df[df.flag_S_donor == 0]
    out.append(f"n complexes: {len(df)} ({len(S)} with S donors)")
    out.append(f"global coord bin {BINS['coord']} A; CN cutoffs 3.0 / 3.2 A; "
               f"*_elem columns: Eu-S edges and S atoms up to {S_EDGE} A, N/O unchanged")

    out.append("\n-- EuNO (reference) --")
    agree("EuNO_H0_coord", "CN_NO_3p2", "requested comparison")
    agree("EuNO_H0_coord", "CN_NO_3p0", "same 3.0 A edge as the coord bin")
    agree("EuNO_H0_coord_Eu", "CN_NO_3p0", "only MST edges touching Eu")
    agree("CN_NO_3p2", "CN_design", "direct CN vs designed CN 9")

    out.append("\n-- EuNOS, global bins --")
    agree("EuNOS_H0_coord", "CN_NOS_3p2", "all complexes")
    agree("EuNOS_H0_coord_Eu", "CN_NOS_3p2", "all complexes")
    agree("EuNOS_H0_coord_Eu", "CN_NOS_3p0", "same 3.0 A edge as the coord bin")
    agree("CN_NOS_3p2", "CN_design", "direct CN vs designed CN 9")
    out.append("\n-- EuNOS, element-specific S edge (B1, main) --")
    agree("EuNOS_H0_coord_Eu_elem", "CN_NOS_elem", "all complexes")
    agree("CN_NOS_elem", "CN_design", "direct CN vs designed CN 9")
    if len(S):
        out.append("  sulfur-binding complexes only:")
        agree("EuNO_H0_coord_Eu", "CN_NO_3p0", "S-donor, EuNO reference", S)
        agree("EuNOS_H0_coord_Eu", "CN_NOS_3p2", "S-donor, global bins", S)
        agree("EuNOS_H0_coord_Eu", "CN_NOS_3p0", "S-donor, global bins", S)
        agree("EuNOS_H0_coord_Eu_elem", "CN_NOS_elem", "S-donor, B1", S)
        agree("CN_NOS_3p2", "CN_design", "S-donor, global", S)
        agree("CN_NOS_elem", "CN_design", "S-donor, B1", S)
    if len(noS):
        out.append(f"  non-S complexes: EuNOS_H0_coord_Eu == EuNO_H0_coord_Eu in "
                   f"{int((noS.EuNOS_H0_coord_Eu == noS.EuNO_H0_coord_Eu).sum())}/{len(noS)}; "
                   f"*_elem == global in {int((noS.EuNOS_H0_coord_Eu_elem == noS.EuNOS_H0_coord_Eu).sum())}/{len(noS)} "
                   f"(coord_Eu) and {int((noS.CN_NOS_elem == noS.CN_NOS_3p2).sum())}/{len(noS)} (CN)")

    out.append("\nmean(H0_coord - CN) by group:   EuNO vs CN_NO_3p2 | EuNOS vs CN_NOS_3p2 | "
               "EuNOS_Eu vs CN_NOS_3p2 | EuNOS_Eu_elem vs CN_NOS_elem")
    for g, sub in df.groupby("group"):
        out.append(f"    {g:<4} n={len(sub):>3}   {(sub.EuNO_H0_coord - sub.CN_NO_3p2).mean():+6.2f}"
                   f"            {(sub.EuNOS_H0_coord - sub.CN_NOS_3p2).mean():+6.2f}"
                   f"              {(sub.EuNOS_H0_coord_Eu - sub.CN_NOS_3p2).mean():+6.2f}"
                   f"                 {(sub.EuNOS_H0_coord_Eu_elem - sub.CN_NOS_elem).mean():+6.2f}")

    lo, hi = BINS["coord"]
    for name, targets in (("EuNO", ("N", "O")), ("EuNOS", ("N", "O", "S"))):
        pairs, missed, missed_by, eus_elem = {}, 0, {}, 0
        for _, r in df.iterrows():
            z = npz(r)
            el, dth, loc = z[f"{name}_H0_pair_elements"], z[f"{name}_H0_pair_death"], z[f"{name}_H0_pair_local"]
            for (x, y), v in zip(el, dth):
                k = "-".join(sorted((x, y)))
                if lo <= v < hi:
                    pairs[k] = pairs.get(k, 0) + 1
                if k == "Eu-S" and lo <= v < S_EDGE:
                    eus_elem += 1
            els, c = z[f"{name}_elements"], z[f"{name}_coords"]
            eu = int(np.where(els == "Eu")[0][0])
            d = np.linalg.norm(c - c[eu], axis=1)
            adj = set(loc[(loc[:, 0] == eu) | (loc[:, 1] == eu)].ravel()) - {eu}
            for i in np.where((d < hi) & np.isin(els, targets))[0]:
                if i not in adj:
                    missed += 1
                    inc = np.where((loc[:, 0] == i) | (loc[:, 1] == i))[0]
                    k = "-".join(sorted(el[inc[np.argmin(dth[inc])]]))
                    missed_by[k] = missed_by.get(k, 0) + 1
        tot = sum(pairs.values())
        out.append(f"\ncoord-bin {name} MST edges by element pair ({tot} edges):")
        for k, v in sorted(pairs.items(), key=lambda kv: -kv[1]):
            out.append(f"    {k:<6} {v:>5}  ({100*v/tot:.1f}%)")
        if name == "EuNOS":
            out.append(f"    Eu-S edges in [{lo}, {S_EDGE}) under B1: {eus_elem}")
        out.append(f"{'/'.join(targets)} atoms within {hi} A of Eu NOT MST-adjacent to Eu: {missed} (attached via: "
                   + ", ".join(f"{k} {v}" for k, v in sorted(missed_by.items(), key=lambda kv: -kv[1])) + ")")

    dS, extra = [], []
    for _, r in df.iterrows():
        z = npz(r)
        ai, els, c = z["EuNOS_atom_index"], z["EuNOS_elements"], z["EuNOS_coords"]
        eu = int(np.where(els == "Eu")[0][0])
        d = np.linalg.norm(c - c[eu], axis=1)
        des = set(map(int, z["donor_index"]))
        for k in np.where(els == "S")[0]:
            if int(ai[k]) in des:
                dS.append(float(d[k]))
            elif d[k] < 3.6:
                extra.append((int(r.idx), r.family, int(ai[k]), float(d[k])))
    if dS:
        dS = np.sort(np.array(dS))
        edges = [0, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.6, np.inf]
        out.append(f"\ndesigned Eu-S donor distances (n={len(dS)}):")
        for x, y in zip(edges[:-1], edges[1:]):
            out.append(f"    [{x:.1f}, {y:.1f}) : {int(((dS >= x) & (dS < y)).sum())}")
        bnd = dS[dS < 3.6]
        out.append(f"    below 3.6 A (bonded): n={len(bnd)}, min {bnd.min():.2f}, median {np.median(bnd):.2f}, "
                   f"max {bnd.max():.2f}; >=3.0 A: {int((bnd >= 3.0).sum())}; >=3.2 A: {int((bnd >= 3.2).sum())}; "
                   f"within B1 {S_EDGE} A: {int((bnd < S_EDGE).sum())}")
        out.append("    sorted: " + ", ".join(f"{v:.2f}" for v in dS))
        out.append("non-designated S within 3.6 A of Eu: " + ("; ".join(
            f"idx {i} {fam} S(xyz {a}) {v:.2f} A{' (counted in CN_NOS_elem)' if v < S_EDGE else ''}"
            for i, fam, a, v in extra) or "none"))

    for col in ("CN_NO_3p2", "CN_NOS_3p2", "CN_NOS_elem"):
        off = df[df[col] != df.CN_design]
        out.append(f"\ncomplexes with {col} != CN_design: {len(off)}")
        for _, r in off.sort_values(["group", "idx"]).iterrows():
            out.append(f"    idx {int(r.idx):>3} {r.group:<4} {r.family:<32} CN_NO_3p2 {int(r.CN_NO_3p2)} "
                       f"CN_NOS_3p2 {int(r.CN_NOS_3p2)} CN_NOS_elem {int(r.CN_NOS_elem)} "
                       f"donors(elem) {int(r.n_donor_within_elem)}/{int(r.denticity)} "
                       f"waters<3.2 {int(r.n_water_within_3p2)}/{int(r.n_water)} S_donor {int(r.flag_S_donor)}")

    for flag, ncol in (("flag_donor_detached", "n_donor_within_elem"), ("flag_donor_detached_3p2", "n_donor_within_3p2")):
        det = df[df[flag]]
        out.append(f"\n{flag} ({ncol} < denticity): {len(det)}")
        for _, r in det.sort_values(["group", "idx"]).iterrows():
            out.append(f"    idx {int(r.idx):>3} {r.group:<4} {r.family:<32} donors {int(r[ncol])}/{int(r.denticity)}  "
                       f"d_Eu_donor_max {r.d_Eu_donor_max:.2f}  S_donor {int(r.flag_S_donor)}")

    for sub in SUBSETS:
        out.append(f"bars at/after 8 A ({sub}): H0 deaths {int(df[f'{sub}_H0_beyond8'].sum())} "
                   f"in {int((df[f'{sub}_H0_beyond8'] > 0).sum())} complexes; "
                   f"H1 births {int(df[f'{sub}_H1_beyond8'].sum())}")
    return "\n".join(out)


def plot(npz_path, png, title, subsets=("heavy", "EuNO")):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from persim import plot_diagrams
    z = np.load(npz_path)
    fig, axes = plt.subplots(1, len(subsets), figsize=(5.5 * len(subsets), 5))
    for ax, s in zip(axes, subsets):
        plot_diagrams([z[f"{s}_H0"], z[f"{s}_H1"]], ax=ax, show=False, labels=["H0", "H1"],
                      title=f"{s} ({len(z[f'{s}_elements'])} atoms)")
        for edge in (1.8, 3.0):
            ax.axhline(edge, ls=":", lw=0.8, c="grey")
            ax.axvline(edge, ls=":", lw=0.8, c="grey")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", nargs="+", default=list(SOURCES), choices=list(SOURCES))
    ap.add_argument("--nproc", type=int, default=9)
    a = ap.parse_args()

    table = pd.read_csv(os.path.join(BASE, "results", "donor_table.csv")).set_index("idx")
    waterfix_idx = set(pd.concat(pd.read_csv(os.path.join(BASE, "results-waterfix", g, "summary.csv"))
                                 for g in os.listdir(os.path.join(BASE, "results-waterfix"))
                                 if g.startswith("G")).idx)

    for key in a.sources:
        tree = os.path.join(BASE, SOURCES[key])
        out = os.path.join(PD, key)
        os.makedirs(out, exist_ok=True)
        summ = pd.concat(pd.read_csv(os.path.join(tree, g, "summary.csv"))
                         for g in sorted(os.listdir(tree)) if g.startswith("G"))
        tasks, skipped = [], []
        for _, r in summ.sort_values("idx").iterrows():
            idx = int(r.idx)
            if pd.isna(r.best_seed):
                why = ("no converged structure in results/; recovered in results-waterfix (see pd/waterfix)"
                       if key == "results" and idx in waterfix_idx else "no converged seed")
                skipped.append((idx, r.group, why)); continue
            seed = int(r.best_seed)
            t = table.loc[idx]
            donors = [int(x) for x in str(t.donor_atoms).split()]
            delems = [re.match(r"(O|N|S)", x).group(1) for x in t.donor_rule.split("|")[2].split(".")]
            tasks.append(dict(
                source=SOURCES[key], group=r.group, idx=idx, family=r.family, SMILES=r.SMILES,
                split=r.split, seed=seed, energy_Eh=float(r.E_min_Eh), denticity=int(r.denticity),
                n_water=int(r.n_water), total_charge=int(r.total_charge),
                n_converged=int(r.n_converged), n_heavy=int(t.n_heavy), donors=donors, donor_elems=delems,
                xyz=os.path.join(tree, r.group, str(idx), f"seed_{seed:02d}_opt.xyz"),
                npz=os.path.join(out, "raw", r.group, f"idx{idx}_seed{seed:02d}.npz")))

        with Pool(a.nproc) as p:
            res = p.map(process, tasks, chunksize=1)
        rows = [r for r, _ in res if r is not None]
        skipped += [s for _, s in res if s is not None]
        df = pd.DataFrame(rows).sort_values("idx").reset_index(drop=True)
        df.to_csv(os.path.join(out, "features.csv"), index=False)
        pd.DataFrame(skipped, columns=["idx", "group", "reason"]).sort_values("idx") \
          .to_csv(os.path.join(out, "skipped.csv"), index=False)

        os.makedirs(os.path.join(out, "plots"), exist_ok=True)
        for ids, subsets, suffix in ((PLOT_IDX[key], ("heavy", "EuNO"), ""),
                                     (PLOT_S[key], ("EuNO", "EuNOS"), "_EuNO_vs_EuNOS")):
            for idx in ids:
                row = df[df.idx == idx]
                if len(row):
                    r = row.iloc[0]
                    plot(os.path.join(out, "raw", r.group, f"idx{idx}_seed{int(r.seed):02d}.npz"),
                         os.path.join(out, "plots", f"idx{idx}{suffix}.png"),
                         f"idx {idx} | {r.group} {r.family} | d={r.denticity}, {r.n_water} H2O | seed {r.seed}",
                         subsets)

        rep = sanity(df, out)
        open(os.path.join(out, "sanity_check.txt"), "w").write(rep + "\n")
        print(f"\n===== {key}: {len(df)} complexes written, {len(skipped)} skipped -> {out}")
        for s in skipped:
            print(f"  skipped idx {s[0]} ({s[1]}): {s[2]}")
        print(rep)


if __name__ == "__main__":
    main()
