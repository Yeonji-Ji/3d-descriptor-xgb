# Does constraining the actinyl axis fix the bent O=M=O geometries?

Small test: 5 complexes × 4 seeds, re-run with the O=M=O angle constrained to 180°, compared with the existing unconstrained runs of the same seeds. The parent folder was only read.

## Step 1: the 5 complexes

Values below are the parent set's unconstrained result (`../structures.csv`).

| Complex | Chosen as | Atoms | Donors + waters | Angle | M=O (Å) | Spin excess |
|---|---|---|---|---|---|---|
| 17_U6 | U(VI), severely bent | 122 | 3 + 2 | 97.2° | 1.670 / 1.682 | 0.00 |
| 259_U6 | U(VI), middle range | 75 | 1 + 4 | 162.0° | 1.711 / 1.711 | 0.00 |
| 280_Np5 | Np(V) | 74 | 3 + 2 | 109.2° | 1.716 / 1.713 | +0.44 |
| 12_Pu6 | Pu(VI), high spin excess | 122 | 3 + 2 | 94.2° | 1.648 / 1.631 | +0.78 |
| 3_U6 | 6 equatorial donors (build-flagged) | 101 | 6 + 0 | 101.5° | 1.690 / 1.701 | 0.00 |

On Pu(VI): the spin excesses closest to +1.0 (0.89–1.68) all belong to the 234–240-atom tripodal complexes, so 12_Pu6 (+0.78) was taken instead — it is also one of the worst bent cases.

## Step 2: the runs

Same start structures as the parent run (`../starts/`), seeds 1–4, same charge and uhf, gas phase:

```
xtb in.xyz --gxtb --opt --chrg <q> --uhf <n> --input constrain.inp
```

with, per complex, the verified 1-based indices (checked in every file by element and M=O distance, not assumed — e.g. 17_U6 gives metal 42, oxo 43 and 44):

```
$constrain
  force constant=1.0
  angle: <O1>,<M>,<O2>,180.0
$end
```

**All 20 runs converged at 180.0°, so the 178° fallback was never needed.** Run time 33–212 s per structure.

## Step 3: comparison

Best converged run of each mode per complex (full per-seed table in `results.csv`):

| Complex | Mode | Angle | M=O (Å) | eq CN | Donor off-plane | Metal charge | Spin excess | Final grad | ΔE (kcal/mol) |
|---|---|---|---|---|---|---|---|---|---|
| 17_U6 | unconstrained | 97.2° | 1.670 / 1.682 | 5 | 40.7° | +2.60 | 0.00 | 0.0004 | |
| | constrained | 179.6° | 1.692 / 1.722 | 5 | 4.6° | +2.45 | 0.00 | 0.0006 | +4.98 |
| 259_U6 | unconstrained | 159.7° | 1.717 / 1.701 | 5 | 8.3° | +2.38 | 0.00 | 0.0007 | |
| | constrained | 179.4° | 1.713 / 1.696 | 5 | 2.1° | +2.35 | 0.00 | 0.0006 | +2.49 |
| 280_Np5 | unconstrained | 174.5° | 1.724 / 1.734 | 5 | 10.2° | +2.09 | +0.46 | 0.0006 | |
| | constrained | 179.6° | 1.737 / 1.739 | 5 | 7.3° | +2.23 | +0.44 | 0.0004 | +4.44 |
| 12_Pu6 | unconstrained | 94.2° | 1.648 / 1.631 | 5 | 34.9° | +2.80 | +0.78 | 0.0007 | |
| | constrained | 179.5° | 1.663 / 1.659 | 5 | 12.9° | +2.55 | +0.75 | 0.0005 | +14.32 |
| 3_U6 | unconstrained | 101.5° | 1.690 / 1.701 | 6 | 48.4° | +2.61 | 0.00 | 0.0007 | |
| | constrained | 179.5° | 1.764 / 1.683 | 5 | 29.5° | +2.72 | 0.00 | 0.0004 | +13.69 |

ΔE compares the best constrained structure with the best unconstrained one; under a constraint this is not a true minimum, so it measures how hard the complex resists a linear axis. Per-seed penalties (same seed, both modes) scatter much more — mean 2.8 kcal/mol for 259_U6 up to 23 kcal/mol for 12_Pu6 — because each seed relaxes into a different local minimum.

## Step 4: the four answers

**1. Does the constraint hold the axis linear without distorting the rest? Yes to the axis; partly to the rest.**
Every run lands at 178.7–179.8°. The equatorial plane usually improves with it: the worst donor deviation drops from 40.7° to 4.6° (17_U6), 8.3° to 2.1° (259_U6) and 34.9° to 12.9° (12_Pu6). But the belt pays for it: **equatorial CN drops from 5 to 4 in 8 of the 20 constrained runs** (1 of 20 unconstrained), so a water or donor is pushed out. Averaged over seeds, donor deviation still exceeds 17° for 17_U6, 280_Np5 and 12_Pu6.

**2. Is the residual gradient small, or does the structure fight the constraint? Gradient small, but the bent cases do fight.**
All 20 converged with a final gradient of 0.0007 Eh/bohr or less. The resistance shows up as energy instead: +2.5 kcal/mol for the mildly bent 259_U6, but +13.7 and +14.3 kcal/mol for 3_U6 and 12_Pu6, which started near 95–101°.

**3. Does M=O land in 1.75–1.80 Å? No.**
Only 2 of 40 constrained M=O bonds fall in that window. The typical range is 1.65–1.77: Pu(VI) 1.65–1.71, U(VI) 1.69–1.77, Np(V) 1.71–1.77. The constraint straightens the axis but leaves the M=O bonds as short as before — for Pu(VI) about 0.1 Å shorter than a real actinyl.

**4. Does the 6-donor case work? No — it is hopeless with this model.**
3_U6 reaches a linear axis, but its donors stay ~39° out of the equatorial plane on average (one seed 73°), equatorial CN falls from 6 to 5 (4 in one seed), and it pays ~14 kcal/mol. Six donors do not fit an actinyl equatorial belt, and forcing the axis straight just moves the strain into the ligand.

**Is the full +5/+6 set worth rerunning this way? Not as it stands.**
The constraint does the one thing it promises — a linear axis at about 100–200 s per structure — and it also tidies the equatorial plane in the mildly bent cases. But it does not fix the two other problems flagged in the parent set: M=O stays outside 1.70–1.90 Å just as often, and the spin excess is unchanged (Pu +0.78 → +0.75, Np +0.46 → +0.44). It also introduces a new defect, the equatorial CN dropping to 4 in 40% of runs, and the resulting structures are no longer minima, so their energies cannot be used to rank conformers.

A rerun makes sense only if a linear axis is specifically required for the descriptors, and then only for complexes with 3 or fewer equatorial donors that are already above ~150°. The better next step is to check one constrained structure against DFT — as was done for Am(TEDGA)3 in `xtb-gxtb-test/` — since the short M=O and the spin excess suggest the actinyl electronic structure itself is off, not only the angle.

## Files

- `run_constrained.py` (env-gxtb): verifies the atom indices per file, writes `constrain.inp`, runs the 20 optimisations; `--angle 178` and `--out runs_178` for the fallback, `--only` to restrict.
- `compare.py` (env-pd): reuses `../analyze.py` (and through it `xtb-complex-gxtb/qc.py`) so every quantity is defined as in the parent QC; writes `results.csv`.
- `runs/<name>_s<seed>/`: `in.xyz`, `constrain.inp`, `xtbopt.xyz`, `xtb.log`.
- `results.csv`: one row per run (both modes) — angle, both M=O, equatorial CN, donor and water off-plane deviation, metal charge, spin excess, final gradient, energy, convergence.
- `run_log_runs.csv`, `run.log`: run bookkeeping.
