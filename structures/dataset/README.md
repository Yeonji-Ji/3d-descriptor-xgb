# Complex dataset: usable g-xTB structures paired with the experimental rows

One folder holding every g-xTB structure that passed the inclusion rules, plus the experimental rows they cover. Source folders were only read.

## What is in it

- **1471 structures** for 282 ligands: 1431 at +3 and 40 at +4.
- **7053 experimental rows** (trainVal 6604, test 449), 87.3% of the 8075 rows with oxidation state 3–6. The original train/test split is untouched.
- Every structure is used by at least one row.

## Step 1: structure selection

### +3 (source `xtb-complex-gxtb`): 1431 of 1459 kept

| Outcome | Combinations | Kept |
|---|---|---|
| clean | 1381 | yes |
| water lost, all design donors still bound | 50 | yes |
| a design donor left | 7 | no |
| CN rose and a design donor left | 4 | no |
| CN rose (not a lost water) | 6 | no |
| run failed (SCF / not converged / no output) | 11 | no |

**The CN judgement call.** 67 combinations changed coordination number with a successful run. Splitting them by direction and by whether every design donor is still inside the donor cutoff (N/O 3.2 Å, S 3.3 Å):

| | Donors intact | A donor left |
|---|---|---|
| **CN fell** | **50 → kept (water loss)** | 7 → dropped |
| **CN rose** | 6 → dropped | 4 → dropped |

- **Water loss (50):** these lost one or more waters while the ligand stayed fully bound, which is the case the inclusion rule allows. 18 of the 67 also carried the explicit `water_left` flag.
- **Donor loss (11 in total):** dropped, since the ligand no longer binds as designed — e.g. 57_Gd ends with 0 of 1 donors bound, 193_Eu with 6 of 8.
- **CN rose with donors intact (6):** 58_Tm, 67_Tm, 101_Tm, 120_Am, 120_Eu, 290_Am. Something extra entered the cutoff rather than anything being lost. They fail no drop condition, but they are not "a lost water" either, so the inclusion rule leaves them out. This is the one borderline group — they are marked in `dropped.csv` with reason "CN rose (not a lost water); outside the inclusion rule" and can be added back with a one-line change if you want them.

### +4 (source `xtb-complex-gxtb-ox4`): 40 of 55 kept

| Metal | Combinations | Kept | Dropped |
|---|---|---|---|
| Th | 28 | 28 (all, by rule) | 0 |
| Np | 7 | 5 (spin excess ≤ 0.5) | 2 |
| Pu | 20 | 7 (spin excess ≤ 0.5) | 13 |

The dropped ones all sit at a spin excess of 1.04–1.06, the group that behaves as a reduced +3 complex (longer M–O, lower metal charge — see that folder's README). The kept Np/Pu are all at 0.15–0.34, a clean gap.

### +5/+6 (source `xtb-complex-gxtb-actinyl`): all 53 dropped

Bent O=M=O (44 of 53 below 170°) and M=O outside 1.70–1.90 Å in 31; the constrained-axis test did not fix the distances or the spin.

## Step 2: row coverage

| Oxidation state | Rows | Covered | Lost | Covered |
|---|---|---|---|---|
| +3 | 6557 | 6461 | 96 | 98.5% |
| +4 | 728 | 592 | 136 | 81.3% |
| +5 | 188 | 0 | 188 | 0% |
| +6 | 602 | 0 | 602 | 0% |
| **total** | **8075** | **7053** | **1022** | **87.3%** |

**Lost rows by reason:**

| Reason | Rows |
|---|---|
| oxidation state 5/6 dropped | 790 |
| spin excess > 0.5 (reduced Pu/Np at +4) | 136 |
| a design donor left | 42 |
| run failed | 29 |
| CN rose (not a lost water) | 13 |
| CN rose and a donor left | 10 |
| no structure was ever built (ligand 219) | 2 |

**Metals that lose rows** (only +3 and +4 shown; every other metal is fully covered):

| Oxidation state | Metal | Rows | Covered | Lost |
|---|---|---|---|---|
| +4 | Pu | 175 | 47 | 128 |
| +4 | Np | 100 | 92 | 8 |
| +3 | Am | 1406 | 1357 | 49 |
| +3 | Eu | 1164 | 1149 | 15 |
| +3 | La | 342 | 331 | 11 |
| +3 | Sm | 327 | 321 | 6 |
| +3 | Ce | 390 | 385 | 5 |
| +3 | Nd | 356 | 351 | 5 |
| +3 | Tm | 200 | 197 | 3 |
| +3 | Er, Gd | 299, 324 | 298, 323 | 1 each |

Pu(IV) is the one badly hit case: 73% of its rows have no structure, because 13 of its 20 combinations were the reduced ones.

## What is missing

1. **All +5/+6 chemistry (790 rows, 78% of the losses).** Uranyl, neptunyl and plutonyl are absent, so the dataset covers no actinyl chemistry at all — that is 602 U(VI) rows plus Np(V)/Np(VI)/Pu(VI). g-xTB bends O=M=O in these complexes and constraining the axis does not repair the M=O distances or the spin, so nothing here was trusted.
2. **Reduced Pu(IV) and Np(IV) (136 rows).** Only 47 of 175 Pu(IV) rows survive. Any Pu(IV) model built on this dataset rests on 7 ligands.
3. **Ligand 219 (2 rows), `CCCCNC(=O)CC(=O)NCCCC`, a malonamide.** It never embedded in the original GFN2 work, so it has no structure at any oxidation state.
4. **Structures with an imperfect but usable shell.** 50 of the kept +3 structures lost one or more waters, so their CN is below the CN 9 design. They are marked in `index.csv` (`kept_reason` = "water lost, all design donors still bound") and can be filtered out if the descriptors need a complete water shell.
5. **Spin caveat carried over.** Kept structures include open-shell metals whose metal spin exceeds the requested uhf by up to 0.5 (Np/Pu at +4, and the late lanthanides at +3 are exact). `index.csv` has `spin_excess` for filtering.

## Files

- `structures/<idx>_<El><ox>.xyz` — 1471 structures, copied unchanged from the source folders (+3 from `runs/<name>/xtbopt.xyz`, which is that set's single optimised structure; +4 from `best/`).
- `index.csv` — one row per structure: ligand idx, SMILES, element, atomic number, oxidation state, charge, uhf, source folder, family, denticity, CN, mean M–O and M–N, metal charge and spin, spin excess, donors and waters bound, energy, `kept_reason`, and the original QC flags.
- `rows.csv` — the 7053 covered experimental rows with `structure_id` and the original `split`. All original columns are kept except the 1671 `CircularFP_*` fingerprint columns, which are recomputable from SMILES.
- `dropped.csv` — every excluded combination (96: 53 actinyl, 28 at +3, 15 at +4) and every lost row group (98), each with its reason.
- `assemble.py` — rebuilds all of the above; rerunning it gives byte-identical files.
