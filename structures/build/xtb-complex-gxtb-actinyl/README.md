# g-xTB build of the +5 / +6 actinyl complexes

[MO2(L)(H2O)n] complexes for the Oxidation_state 5 and 6 rows, built here (the +3/+4 builder has no oxo atoms) and optimised with g-xTB in gas phase. Other folders were only read.

## Step 0 (side check on the +4 set): yes, the high-spin-excess complexes have longer M–O

In `xtb-complex-gxtb-ox4/structures.csv`, the Pu and Np complexes split into two groups by metal spin excess, and the high group has clearly longer M–O and a lower metal charge:

| | Pu, excess ≈ +1.0 | Pu, excess ≈ +0.25 | Np, ≈ +1.0 | Np, ≈ +0.25 |
|---|---|---|---|---|
| Combinations | 13 | 7 | 2 | 5 |
| Mean M–O (Å) | 2.473 | 2.370 | 2.517 | 2.392 |
| Mean metal charge | +2.31 | +2.68 | +2.29 | +2.62 |

The two groups contain different ligands, so the comparison was repeated inside one family: for DGA complexes with 6 waters, the high-excess ones (12, 14, 15, 17, 253_Pu) sit at 2.434–2.468 Å and the low-excess ones (16, 32, 255, 280_Pu) at 2.364–2.371 Å. Same ligand type, same water count, still about 0.09 Å longer with a lower charge. **This is consistent with the metal being effectively reduced towards +3 in those calculations.**

## Headline for the actinyl set

- **53 combinations** (U(VI) 30, Np(V) 11, Np(VI) 4, Pu(VI) 7, Pa(V) 1); **all 53 have a converged structure** in `best/`.
- **1190 seed runs, 1105 converged (92.9%)**, about 65 CPU-h.
- **Only 1 combination (32_U6) passes every check.** The dominant problem is the actinyl geometry itself: **44 of 53 have O=M=O below 170°** and 31 have an M=O outside 1.70–1.90 Å.
- The bending is a real g-xTB result, not a build artifact — see "Is the bending real?".

## Step 1: the list

790 rows with Oxidation_state 5 or 6, 30 ligands, 53 unique (SMILES, atomic number, oxidation state) combinations. All ligands were already in the `xtb-complex-1` donor table and all are neutral, so the complex charge is the actinyl charge: +2 for VI (41 combinations), +1 for V (12).

## Step 2: the builder

The geometry is assembled directly (same idea as the geometric water shell in `xtb-complex-1/src/build.py`):

1. ligand + metal only, embedded with the shared ETKDGv3 settings and seed;
2. equatorial plane = best-fit plane through the donor directions, its normal being the O=M=O axis;
3. two oxo O at ±1.78 Å along that axis, so the start is exactly linear;
4. waters in the equatorial plane by angular repulsion, M–O 2.45 Å, H–O–H 104.5°.

Atom order: ligand heavy, metal, oxo O ×2, water O, then H. Waters: n_water = 5 − denticity (min 0). 20 seeds per combination; all 1060 embedded.

**Check before running** (one U, one Np, 2 seeds each): the built structures had M=O 1.78 Å, O=M=O 180.0°, waters equatorial at 2.45 Å, equatorial CN 5. After optimisation the donors and waters settled into the equatorial plane (mostly ≤6° out). No builder change was needed.

**Build flags** (9 combinations):
- `Pa(V) may not form a true PaO2+ ion`: 174_Pa5, built like the others as instructed.
- `equatorial denticity > 5`: 3_U6, 4_U6 (bisDGA, 6 donors) and the six 9-dentate tripodal trisDGA combinations (33, 34, 281 × U6/Pu6). Five or more donors cannot sit in an actinyl equatorial plane, and these are also the worst geometries below.

## Step 3: runs

`xtb in.xyz --gxtb --opt --chrg <q> --uhf <n>`, gas phase, 1 thread each, 9 in parallel, resumable. uhf: U(VI) 0, Np(VI) 1, Pu(VI) 2, Np(V) 2, Pa(V) 0.

**63 seeds hit "SCF not converged"**; the rest of the 85 failures are unconverged optimisations. No combination lost all its seeds.

## Step 4: seed check

- **13 combinations had their best seed in 15–20** and got seeds 21–30 (one extra round, as specified).
- **5 improved**; 8 still peak in 15–20 and keep the `late_best_seed` flag: 12_U6, 15_Np5, 33_U6, 44_Np5, 82_Np5, 173_Pu6, 184_U6, 258_U6.
- **Gap to the second-best seed:** median 1.57 kcal/mol, mean 3.09, range 0.03–18.7.

## Step 5: QC

| Flag | Combinations |
|---|---|
| oxo_angle (O=M=O < 170°) | 44 |
| oxo_dist (M=O outside 1.70–1.90 Å) | 31 |
| spin_mismatch | 22 |
| eq_out_of_plane (donor > 25° off-plane) | 16 |
| cn_change | 13 |
| late_best_seed | 8 |
| short_bond (< 2.0 Å) | 7 |
| eq_CN_range (outside 4–6) | 4 |
| water_left | 2 |
| missing / scf_fail (whole combination) | 0 |

**O=M=O angle** — the main finding:

| Metal | Combinations | Median angle | Min | ≥ 170° |
|---|---|---|---|---|
| Np(V) | 11 | 165.9° | 99.0° | 4 |
| Np(VI) | 4 | 165.3° | 162.7° | 0 |
| Pa(V) | 1 | 155.0° | 155.0° | 0 |
| Pu(VI) | 7 | 152.8° | 94.2° | 3 |
| U(VI) | 30 | 162.5° | 97.2° | 2 |

Across the set: 9 combinations are at 170° or above, 32 lie between 150° and 170°, and 12 fall below 150° — including 9 below 120°. The worst are the build-flagged ones: 3_U6 and 4_U6 (bisDGA, 6 equatorial donors) end at ~101°, essentially cis-dioxo.

**M=O distance:** 29 combinations have at least one M=O shorter than 1.70 Å and 4 longer than 1.90 Å. Mean by metal: Np(VI) 1.679, Pu(VI) 1.844, U(VI) 1.703, Np(V) 1.727, Pa(V) 1.805 Å.

**Equatorial CN** is 5 in 42 combinations and 4 or 6 in 7 more. The 4 outside the 4–6 window (7–8 donors) are all 9-dentate tripodal trisDGA, as expected from the build flag.

**Metal spin excess** (metal Mulliken spin − uhf), the same pattern as the +4 set:

| Metal | uhf | Spin excess (mean / range) |
|---|---|---|
| U(VI) | 0 | 0.00 |
| Pa(V) | 0 | 0.00 |
| Np(VI) | 1 | +0.43 (0.42–0.43) |
| Np(V) | 2 | +0.45 (0.43–0.46) |
| Pu(VI) | 2 | +1.02 (0.75–1.68) |

Closed-shell U(VI) and Pa(V) are exact; every open-shell actinyl carries extra spin density on the metal, compensated by the ligands. Pu(VI) is near +1, the same size of excess that, in the +4 set, came with a longer M–O and a lower metal charge.

## Is the bending real?

Three checks say the bent geometries are genuine g-xTB minima for these complexes, not artifacts of the build:

1. **The method handles bare actinyl correctly.** Starting linear, g-xTB keeps UO2(2+) at 180.0° (U=O 1.778 Å) and [UO2(H2O)5](2+) at 180.0° (1.708 Å).
2. **Straightening does not help.** Taking two bent results, forcing O=M=O back to 180° and re-optimising returns them to exactly the same angle and energy: 260_U6 163.0° → 163.0° (−0.00 kcal/mol), 3_U6 101.5° → 101.5° (−0.03 kcal/mol).
3. **It is not oxo protonation.** Only 2 of 53 complexes have an H within 1.2 Å of an oxo; the other 51 keep bare oxo at 1.71–1.72 Å while still bending.

So with an organic ligand bound, g-xTB prefers a bent actinyl, while real actinyl ions are close to linear. **These structures should not be used as actinyl geometries without checking them against DFT**, in the way the Am(TEDGA)3 comparison was done in `xtb-gxtb-test/`.

## Files

- `build_list.py` (env-gxtb): builds `list.csv`, `tasks.txt`, `build_log.csv`, `starts/<idx>_<El><ox>_s<seed>.xyz`. `--seeds 21 30 --only <names>` appends the extra round.
- `run_gxtb.sh` (env-gxtb): runs g-xTB; resumable; `runs/status.txt` logs uhf, charge, rc and wall time.
- `analyze.py` (env-pd): writes `seeds.csv`, `structures.csv` and `best/`; reuses parsing from `xtb-complex-gxtb/qc.py` and adds the actinyl checks.
- `runs/<idx>_<El><ox>_s<seed>/`, `best/<idx>_<El><ox>.xyz`.
- `structures.csv`: per combination — best seed, energy, gap to second best, oxo distances and angle, equatorial CN and out-of-plane deviations, M–O and M–N, water range, metal charge, spin and spin excess, build flag, QC flags.
- `seeds.csv`: per seed — convergence, energy, gradient, metal charge and spin, wall time.
