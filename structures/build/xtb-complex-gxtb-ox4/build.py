"""Build [Eu(L)(H2O)n] and run GFN2-xTB. Fixed settings per task spec."""
import os, re, subprocess, time
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

EU_Z = 63
CUTOFF = 3.0                      # A, "bound" criterion for donors and waters
XTB_BASE = ["--gfn", "2", "--opt", "--uhf", "0", "--alpb", "water"]


def build_complex(build_smiles, donor_atoms, n_water):
    """ligand + Eu(+3) + n_water O, DATIVE bonds donor->Eu and O_water->Eu."""
    lig = Chem.MolFromSmiles(build_smiles)
    assert lig is not None, "ligand SMILES unparsable"
    em = Chem.RWMol(lig)
    eu = Chem.Atom(EU_Z)
    eu.SetFormalCharge(3)
    eu.SetNoImplicit(True)
    eu_idx = em.AddAtom(eu)
    for d in donor_atoms:
        em.AddBond(int(d), eu_idx, Chem.BondType.DATIVE)
    wat = []
    for _ in range(n_water):
        o = Chem.Atom(8)          # bare O; AddHs turns each into H2O
        oi = em.AddAtom(o)
        em.AddBond(oi, eu_idx, Chem.BondType.DATIVE)
        wat.append(oi)
    m = em.GetMol()
    Chem.SanitizeMol(m)
    m = Chem.AddHs(m)
    return m, eu_idx, list(donor_atoms), wat


def embed(mol, seed):
    ps = AllChem.ETKDGv3()
    ps.randomSeed = int(seed)
    ps.useRandomCoords = True
    ps.maxIterations = 2000
    return AllChem.EmbedMolecule(mol, ps)


def write_xyz(mol, path, comment=""):
    conf = mol.GetConformer()
    with open(path, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n{comment}\n")
        for a in mol.GetAtoms():
            p = conf.GetAtomPosition(a.GetIdx())
            f.write(f"{a.GetSymbol():<3s} {p.x:>14.8f} {p.y:>14.8f} {p.z:>14.8f}\n")


def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0]); xyz.append(tuple(float(v) for v in p[1:4]))
    return syms, xyz


def run_xtb(xyz_path, charge, workdir, log_path, timeout=None):
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_STACKSIZE="4G")
    cmd = ["xtb", os.path.basename(xyz_path)] + XTB_BASE + ["--chrg", str(int(charge))]
    t0 = time.time()
    with open(log_path, "w") as lg:
        try:
            p = subprocess.run(cmd, cwd=workdir, stdout=lg, stderr=subprocess.STDOUT,
                               env=env, timeout=timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            return "timeout", None, time.time() - t0
    wall = time.time() - t0
    log = open(log_path, errors="replace").read()
    e = None
    mm = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", log)
    if mm:
        e = float(mm[-1])
    if "GEOMETRY OPTIMIZATION CONVERGED" in log and "normal termination" in log:
        return "converged", e, wall
    if "FAILED TO CONVERGE GEOMETRY" in log:
        return "not_converged", e, wall
    return f"xtb_error(rc={rc})", e, wall


def xtb_error_message(log_path):
    """Pull xtb's real reason out of the log, skipping its '####' banner lines."""
    if not os.path.exists(log_path):
        return ""
    txt = open(log_path, errors="replace").read()
    msgs = re.findall(r"^\s*-\d+-\s*(.+?)\s*$", txt, re.M)
    if msgs:
        return "; ".join(dict.fromkeys(m for m in msgs))[:240]
    for pat in (r"\[ERROR\]\s*(.+)", r"\[WARNING\]\s*(.+)"):
        m = re.findall(pat, txt)
        if m:
            return m[-1].strip()[:240]
    lines = [l.strip() for l in txt.splitlines() if l.strip() and set(l.strip()) != {"#"}]
    return (lines[-1][:240] if lines else "")


def dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def coord_counts(xyz, eu_idx, donor_atoms, water_O):
    eu = xyz[eu_idx]
    dd = [dist(eu, xyz[i]) for i in donor_atoms]
    dw = [dist(eu, xyz[i]) for i in water_O]
    nd = sum(1 for d in dd if d <= CUTOFF)
    nw = sum(1 for d in dw if d <= CUTOFF)
    return nd, nw, (max(dd) if dd else 0.0), len(dw) - nw


# ---------------------------------------------------------------- geometric water shell
# Recovery path for complexes ETKDG cannot embed when 5-8 waters are dative-bonded to Eu
# at once (diagnosed: all such failures embed fine at <=2 waters). Here the ligand+Eu is
# embedded by the SAME ETKDGv3 call/seed as everywhere else, and the waters are then placed
# on the vacant coordination directions instead of being handed to distance geometry.
import math, random


def _unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v) if n else (0.0, 0.0, 1.0)


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _scale(a, k):
    return tuple(x * k for x in a)


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _fibonacci_sphere(n=512):
    pts, ga = [], math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        y = 1 - 2 * i / float(n - 1)
        r = math.sqrt(max(0.0, 1 - y * y))
        th = ga * i
        pts.append((math.cos(th) * r, y, math.sin(th) * r))
    return pts


def place_waters(eu_pos, occupied_dirs, n_water, lig_positions, d_ow=2.45, seed=0):
    """Choose n_water coordination directions by repulsion relaxation on the unit sphere
    (donor directions fixed, waters mobile, ligand atoms as weak repulsors), then build
    H2O on each: O at d_ow, H-O-H 104.5 deg, bisector pointing away from Eu.

    Greedy max-min selection was tried first and packs points far too tightly
    (38 deg between neighbours where 9-coordination wants ~70), so relaxation it is.
    """
    rng = random.Random(seed)
    fixed = [_unit(d) for d in occupied_dirs]
    # weak repulsors: ligand atoms near the metal, projected onto the sphere
    near = []
    for q in lig_positions:
        v = _sub(q, eu_pos)
        r2 = _dot(v, v)
        if r2 < 36.0 and r2 > 1e-6:
            near.append((_unit(v), 0.35))
    mob = []
    for _ in range(n_water):
        v = _unit((rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1)))
        mob.append(v)

    for step in range(400):
        step_size = 0.25 * (1.0 - step / 400.0) + 0.02
        new_mob = []
        for i, a in enumerate(mob):
            fx = [0.0, 0.0, 0.0]
            others = ([(d, 1.0) for d in fixed]
                      + [(d, 1.0) for j, d in enumerate(mob) if j != i]
                      + near)
            for b, wgt in others:
                d = _sub(a, b)
                r2 = _dot(d, d)
                if r2 < 1e-9:
                    d = (rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1))
                    r2 = _dot(d, d)
                f = wgt / max(r2, 0.05)
                fx = [fx[k] + d[k] * f for k in range(3)]
            # keep the move tangential, then renormalise back onto the sphere
            rad = _dot(fx, a)
            tan = tuple(fx[k] - rad * a[k] for k in range(3))
            new_mob.append(_unit(_add(a, _scale(tan, step_size))))
        mob = new_mob

    out = mob
    waters = []
    half = math.radians(104.5 / 2.0)
    placed = list(lig_positions)
    for d in out:
        o = _add(eu_pos, _scale(d, d_ow))
        u = _unit(d)                                   # Eu -> O, H's go on the far side
        ref = (0.0, 0.0, 1.0) if abs(u[2]) < 0.9 else (1.0, 0.0, 0.0)
        v0 = _unit(_cross(u, ref))
        v1 = _unit(_cross(u, v0))
        bestH, bestd = None, -9e9
        for k in range(18):
            a = 2 * math.pi * k / 18.0
            v = _unit(_add(_scale(v0, math.cos(a)), _scale(v1, math.sin(a))))
            hs = [_add(o, _scale(_add(_scale(u, math.cos(half)), _scale(v, sgn * math.sin(half))), 0.97))
                  for sgn in (1, -1)]
            md = min((sum((x - y) ** 2 for x, y in zip(h, q)) for h in hs for q in placed),
                     default=99.0)
            if md > bestd:
                bestH, bestd = hs, md
        waters.append((o, bestH[0], bestH[1]))
        placed.extend([o, bestH[0], bestH[1]])
    return waters


def build_complex_geom(build_smiles, donor_atoms, n_water, seed):
    """Same complex, same atom ordering and same ETKDGv3 seed as build_complex/embed, but the
    water shell is placed geometrically. Returns (mol_with_conformer, eu, donors, water_O) or None."""
    full, eu, dn, wat = build_complex(build_smiles, donor_atoms, n_water)

    # ligand + Eu only (no waters) -> identical heavy-atom order, embed with the same seed
    lig = Chem.MolFromSmiles(build_smiles)
    em = Chem.RWMol(lig)
    a = Chem.Atom(EU_Z); a.SetFormalCharge(3); a.SetNoImplicit(True)
    eu2 = em.AddAtom(a)
    for d in donor_atoms:
        em.AddBond(int(d), eu2, Chem.BondType.DATIVE)
    core = em.GetMol()
    Chem.SanitizeMol(core)
    core = Chem.AddHs(core)
    if embed(core, seed) != 0:
        return None
    cc = core.GetConformer()

    conf = Chem.Conformer(full.GetNumAtoms())
    nheavy_lig = lig.GetNumAtoms()
    # heavy atoms: ligand 0..L-1 then Eu -> same indices in both molecules
    for i in range(nheavy_lig + 1):
        p = cc.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x, p.y, p.z))
    # hydrogens on the ligand: AddHs appends them grouped by parent, in parent order, in both
    hcore = [x.GetIdx() for x in core.GetAtoms() if x.GetAtomicNum() == 1]
    hfull = [x.GetIdx() for x in full.GetAtoms()
             if x.GetAtomicNum() == 1 and x.GetNeighbors()[0].GetIdx() < nheavy_lig]
    assert len(hcore) == len(hfull), f"H mapping mismatch {len(hcore)} vs {len(hfull)}"
    for hc, hf in zip(hcore, hfull):
        p = cc.GetAtomPosition(hc)
        conf.SetAtomPosition(hf, Point3D(p.x, p.y, p.z))

    eu_pos = tuple(cc.GetAtomPosition(eu2))
    occ = [_unit(_sub(tuple(cc.GetAtomPosition(int(d))), eu_pos)) for d in donor_atoms]
    ligpos = [tuple(cc.GetAtomPosition(i)) for i in range(core.GetNumAtoms()) if i != eu2]
    shells = place_waters(eu_pos, occ, n_water, ligpos, seed=seed)
    for oi, (o, h1, h2) in zip(wat, shells):
        hs = [x.GetIdx() for x in full.GetAtomWithIdx(oi).GetNeighbors() if x.GetAtomicNum() == 1]
        assert len(hs) == 2, f"water {oi} has {len(hs)} H"
        conf.SetAtomPosition(oi, Point3D(*o))
        conf.SetAtomPosition(hs[0], Point3D(*h1))
        conf.SetAtomPosition(hs[1], Point3D(*h2))
    full.RemoveAllConformers()
    full.AddConformer(conf, assignId=True)
    return full, eu, dn, wat
