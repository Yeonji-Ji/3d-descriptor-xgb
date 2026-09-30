"""Steps 1-2: list the Oxidation_state 5/6 combinations and build actinyl [MO2(L)(H2O)n] conformers.

conda activate env-gxtb && python xtb-complex-gxtb-actinyl/build_list.py [--seeds 1 20] [--only NAME ...]

The +3/+4 builder has no oxo atoms, so the geometry is assembled here (same idea as the geometric
water shell in xtb-complex-1/src/build.py):
  1. ligand + metal only, embedded with the SAME ETKDGv3 settings/seed as everywhere else;
  2. equatorial plane = best-fit plane through the donor directions (its normal is the O=M=O axis);
  3. two oxo O placed at +-axis * 1.78 A (linear O=M=O by construction);
  4. waters placed in the equatorial plane by angular repulsion, M-O 2.45 A, H-O-H 104.5 deg.

Atom order: ligand heavy, metal, oxo O x2, water O, then all H.
Waters: n_water = 5 - equatorial denticity (min 0). Charge = actinyl charge + ligand charge.
Other folders are read-only (no bytecode written there).

Outputs: list.csv, tasks.txt, build_log.csv, starts/<idx>_<El><ox>_s<seed>.xyz
"""
import argparse
import math
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Geometry import Point3D

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "orig-data").is_dir())   # project root
XTB = HERE.parent / "xtb-complex-1"
ORIG = ROOT / "orig-data" / "0.1-1"
sys.dont_write_bytecode = True
sys.path.insert(0, str(XTB / "src"))
import build as BUILD  # noqa: E402  (embed, write_xyz, ETKDG settings)

SYMBOL = {91: "Pa", 92: "U", 93: "Np", 94: "Pu"}
Z_OF = {v: k for k, v in SYMBOL.items()}
# uhf = unpaired f electrons of the actinyl ion
UHF = {("U", 6): 0, ("Np", 6): 1, ("Pu", 6): 2, ("Np", 5): 2, ("Pa", 5): 0}
ACTINYL_CHARGE = {6: 2, 5: 1}
D_OXO, D_OW, EQ_WATER_MAX = 1.78, 2.45, 5


def combinations():
    d = pd.concat([pd.read_csv(ORIG / f"{n}_dataset.csv",
                               usecols=["SMILES", "Atomic_number", "Oxidation_state"]).assign(split=n)
                   for n in ("trainVal", "test")], ignore_index=True)
    d = d[d.Oxidation_state.isin([5, 6])]
    return (d.groupby(["SMILES", "Atomic_number", "Oxidation_state"], as_index=False)
             .agg(n_rows=("split", "size"), splits=("split", lambda s: "+".join(sorted(set(s))))))


def equatorial_frame(donor_dirs, rng):
    """Axis (O=M=O direction) = normal of the best-fit plane through the donor directions."""
    if len(donor_dirs) >= 2:
        u, s, vt = np.linalg.svd(np.asarray(donor_dirs) - 0.0)
        axis = vt[-1]                                   # smallest singular direction
    else:                                               # one donor: any perpendicular
        d = np.asarray(donor_dirs[0])
        t = np.array([1.0, 0.0, 0.0]) if abs(d[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(d, t)
    axis = axis / np.linalg.norm(axis)
    e1 = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = e1 - axis * np.dot(axis, e1)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    return axis, e1, e2


def water_angles(occupied, n, rng):
    """n angles on the equatorial circle, kept away from the occupied (donor) angles."""
    if n == 0:
        return []
    ang = []
    for _ in range(n):                                  # greedy max-min start
        best, bestgap = 0.0, -1.0
        for k in range(720):
            a = 2 * math.pi * k / 720
            gap = min(abs((a - b + math.pi) % (2 * math.pi) - math.pi) for b in list(occupied) + ang) \
                if (len(occupied) + len(ang)) else math.pi
            if gap > bestgap:
                best, bestgap = a, gap
        ang.append(best)
    for _ in range(200):                                # short 1-D repulsion relaxation
        new = []
        for i, a in enumerate(ang):
            f = 0.0
            for b in list(occupied) + [x for j, x in enumerate(ang) if j != i]:
                d = (a - b + math.pi) % (2 * math.pi) - math.pi
                f += d / max(abs(d), 0.15) ** 2 * 0.02
            new.append((a + f) % (2 * math.pi))
        ang = new
    return ang


def build_one(row, seed):
    """Return (mol with conformer, metal index, donor indices, oxo indices, water O indices) or None."""
    donors = [int(x) for x in row["donors"].split()]
    n_water = row["n_water_eq"]
    lig = Chem.MolFromSmiles(row["build_smiles"])
    L = lig.GetNumAtoms()

    em = Chem.RWMol(lig)
    m = Chem.Atom(Z_OF[row["element"]])
    m.SetFormalCharge(ACTINYL_CHARGE[row["ox"]])
    m.SetNoImplicit(True)
    m_idx = em.AddAtom(m)
    for d in donors:
        em.AddBond(d, m_idx, Chem.BondType.DATIVE)
    oxo = []
    for _ in range(2):
        o = Chem.Atom(8)
        o.SetNoImplicit(True)          # bare oxo, no H
        oi = em.AddAtom(o)
        em.AddBond(oi, m_idx, Chem.BondType.DATIVE)
        oxo.append(oi)
    wat = []
    for _ in range(n_water):
        o = Chem.Atom(8)               # AddHs turns each into H2O
        oi = em.AddAtom(o)
        em.AddBond(oi, m_idx, Chem.BondType.DATIVE)
        wat.append(oi)
    full = em.GetMol()
    Chem.SanitizeMol(full)
    full = Chem.AddHs(full)

    # ligand + metal only, embedded with the shared ETKDG settings
    em2 = Chem.RWMol(lig)
    m2 = Chem.Atom(Z_OF[row["element"]])
    m2.SetFormalCharge(ACTINYL_CHARGE[row["ox"]])
    m2.SetNoImplicit(True)
    c_idx = em2.AddAtom(m2)
    for d in donors:
        em2.AddBond(d, c_idx, Chem.BondType.DATIVE)
    core = em2.GetMol()
    Chem.SanitizeMol(core)
    core = Chem.AddHs(core)
    if BUILD.embed(core, seed) != 0:
        return None
    cc = core.GetConformer()

    conf = Chem.Conformer(full.GetNumAtoms())
    for i in range(L + 1):                                   # ligand heavy + metal share indices
        p = cc.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x, p.y, p.z))
    hcore = [a.GetIdx() for a in core.GetAtoms() if a.GetAtomicNum() == 1]
    hfull = [a.GetIdx() for a in full.GetAtoms()
             if a.GetAtomicNum() == 1 and a.GetNeighbors()[0].GetIdx() < L]
    assert len(hcore) == len(hfull), "ligand H mapping mismatch"
    for hc, hf in zip(hcore, hfull):
        p = cc.GetAtomPosition(hc)
        conf.SetAtomPosition(hf, Point3D(p.x, p.y, p.z))

    M = np.array(cc.GetAtomPosition(c_idx))
    dd = [np.array(cc.GetAtomPosition(d)) - M for d in donors]
    dd = [v / np.linalg.norm(v) for v in dd]
    rng = random.Random(seed)
    axis, e1, e2 = equatorial_frame(dd, rng)
    for oi, sgn in zip(oxo, (1, -1)):
        conf.SetAtomPosition(oi, Point3D(*(M + sgn * D_OXO * axis)))

    occupied = [math.atan2(np.dot(v, e2), np.dot(v, e1)) for v in dd]
    lig_pos = [np.array(cc.GetAtomPosition(i)) for i in range(core.GetNumAtoms()) if i != c_idx]
    placed = lig_pos + [M + s * D_OXO * axis for s in (1, -1)]
    half = math.radians(104.5 / 2)
    for oi, a in zip(wat, water_angles(occupied, n_water, rng)):
        u = math.cos(a) * e1 + math.sin(a) * e2           # equatorial direction
        O = M + D_OW * u
        best, bestd = None, -9e9
        for k in range(18):                                # rotate H-O-H about M->O, keep H's clear
            t = 2 * math.pi * k / 18
            v = math.cos(t) * axis + math.sin(t) * np.cross(axis, u)
            hs = [O + 0.97 * (math.cos(half) * u + sgn * math.sin(half) * v) for sgn in (1, -1)]
            md = min(np.linalg.norm(h - q) for h in hs for q in placed)
            if md > bestd:
                best, bestd = hs, md
        hs_idx = [x.GetIdx() for x in full.GetAtomWithIdx(oi).GetNeighbors() if x.GetAtomicNum() == 1]
        assert len(hs_idx) == 2, f"water {oi} has {len(hs_idx)} H"
        conf.SetAtomPosition(oi, Point3D(*O))
        for hi, h in zip(hs_idx, best):
            conf.SetAtomPosition(hi, Point3D(*h))
        placed.extend([O] + list(best))

    full.RemoveAllConformers()
    full.AddConformer(conf, assignId=True)
    return full, m_idx, donors, oxo, wat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", nargs=2, type=int, default=[1, 20], metavar=("FIRST", "LAST"))
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    seeds = range(args.seeds[0], args.seeds[1] + 1)
    (HERE / "starts").mkdir(parents=True, exist_ok=True)

    combos = combinations()
    dt = pd.read_csv(XTB / "results" / "donor_table.csv").set_index("SMILES")
    rows = []
    for c in combos.itertuples():
        el, ox = SYMBOL[c.Atomic_number], int(c.Oxidation_state)
        d = dt.loc[c.SMILES]
        n_water_eq = max(0, EQ_WATER_MAX - int(d.denticity))
        rows.append({"name": f"{d.idx}_{el}{ox}", "idx": int(d.idx), "element": el, "ox": ox,
                     "Z": c.Atomic_number, "uhf": UHF[(el, ox)],
                     "chrg": ACTINYL_CHARGE[ox] + int(d.ligand_charge), "SMILES": c.SMILES,
                     "build_smiles": d.build_smiles, "family": d.family, "group": d.group,
                     "splits": c.splits, "n_rows": c.n_rows, "denticity": int(d.denticity),
                     "n_water_eq": n_water_eq, "n_heavy": int(d.n_heavy),
                     "donors": d.donor_atoms, "ligand_charge": int(d.ligand_charge),
                     "eq_donors_gt_5": int(d.denticity) > EQ_WATER_MAX,
                     "flag_build": "Pa(V) may not form a true PaO2+ ion" if el == "Pa" else
                                   ("equatorial denticity > 5" if int(d.denticity) > EQ_WATER_MAX else "")})
    lst = pd.DataFrame(rows).sort_values(["idx", "element", "ox"])
    if args.only:
        lst = lst[lst.name.isin(args.only)]
    print(f"combinations {len(combos)} (ligands {combos.SMILES.nunique()}), building {len(lst)}")

    tasks, log = [], []
    for r in lst.to_dict("records"):
        built = 0
        for seed in seeds:
            got = build_one(r, seed)
            if got is None:
                continue
            mol = got[0]
            name = f"{r['name']}_s{seed:02d}"
            BUILD.write_xyz(mol, HERE / "starts" / f"{name}.xyz",
                            f"{r['name']} seed {seed} chrg={r['chrg']} uhf={r['uhf']} actinyl")
            tasks.append(f"{name} {r['uhf']} {r['chrg']}")
            built += 1
        n_atoms = mol.GetNumAtoms() if built else 0
        log.append({"name": r["name"], "element": r["element"], "ox": r["ox"], "n_atoms": n_atoms,
                    "seeds_built": built, "seeds_requested": len(list(seeds))})
        print(f"  {r['name']:>10s} {n_atoms:>4d} atoms  seeds {built}/{len(list(seeds))}"
              + (f"  [{r['flag_build']}]" if r["flag_build"] else ""))

    lst = lst.assign(n_atoms=lst.name.map({x["name"]: x["n_atoms"] for x in log}))
    if (args.only or args.seeds != [1, 20]) and (HERE / "list.csv").exists():
        old = pd.read_csv(HERE / "list.csv")
        lst = pd.concat([old, lst[~lst.name.isin(old.name)]], ignore_index=True)
        with open(HERE / "tasks.txt", "a") as f:
            f.write("\n".join(tasks) + "\n")
        pd.concat([pd.read_csv(HERE / "build_log.csv"), pd.DataFrame(log)]).to_csv(HERE / "build_log.csv", index=False)
    else:
        (HERE / "tasks.txt").write_text("\n".join(tasks) + "\n")
        pd.DataFrame(log).to_csv(HERE / "build_log.csv", index=False)
    lst.to_csv(HERE / "list.csv", index=False)
    print(f"\nstart structures {len(tasks)}; metals {lst.groupby(['element','ox','uhf']).size().to_dict()}")
    print("charges:", lst.chrg.value_counts().to_dict(),
          "| equatorial waters:", lst.n_water_eq.value_counts().to_dict())
    print("build flags:", lst[lst.flag_build != ""].name.tolist())


if __name__ == "__main__":
    main()
