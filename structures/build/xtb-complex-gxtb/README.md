# g-xTB rebuild of the +3 metal complexes

All trivalent (SMILES, metal) combinations from `orig-data/0.1-1/` re-optimised with g-xTB in gas phase, starting from the existing GFN2 conformers. No new conformer search. `xtb-complex-1/` was only read.

## Headline

- **1459 combinations** (278 ligands × 20 metals), of which **1381 are clean (94.7%)** and 78 are flagged.
- Two settings problems in the first pass were found and fixed (see "Problems found and fixed"); 438 combinations were rerun.
- **Cost:** 79.6 CPU-h in total (61.6 h first pass, 18.0 h rerun), about 8 h wall time on 9 workers.

## Step 1: the list

- **Rows:** 6557 trivalent rows (Oxidation_state == 3) from trainVal + test. Oxidation states 4, 5 and 6 were skipped.
- **Combinations:** 1461 unique (SMILES, atomic number) pairs; 279 unique ligands.
- **Matched:** 1459 combinations / 278 ligands. SMILES matched the `xtb-complex-1` donor table exactly (string match), so no canonicalisation was needed.
- **Unmatched (in `unmatched.csv`):** `CCCCNC(=O)CC(=O)NCCCC` (idx 219, malonamide) as Eu and as Am. That ligand produced no converged GFN2 structure, so there is nothing to start from.
- **Start structures:** the lowest-energy seed of that ligand in `xtb-complex-1` (`best_seed` in the group's `summary.csv`), copied unchanged except that the Eu label became the target element. Coordinated waters were kept. Written to `starts/<idx>_<El>.xyz`, listed in `list.csv`.

| Metal | Combinations | | Metal | Combinations | | Metal | Combinations |
|---|---|---|---|---|---|---|---|
| Eu | 247 | | Sm | 93 | | Yb | 63 |
| Am | 167 | | Ce | 90 | | Ho | 62 |
| La | 103 | | Dy | 81 | | Tm | 58 |
| Nd | 97 | | Gd | 79 | | Cm | 7 |
| Pr | 94 | | Er | 77 | | Pu | 3 |
| | | | Tb | 68 | | Cf | 2 |
| | | | Lu | 66 | | Np, Pm | 1 each |

## Step 2: the runs

`xtb in.xyz --gxtb --opt --chrg <q> --uhf <n>`, gas phase, 1 thread each, 9 in parallel (`run_gxtb.sh`, env `env-gxtb`). Reruns skip any combination whose log already says the geometry optimisation converged. Output in `runs/<idx>_<El>/`.

- **Charge `q`:** the complex charge from `xtb-complex-1/results/donor_table.csv` — +3 for 1423 combinations, and +2 / 0 / −1 for the 36 whose ligand is charged.
- **Spin `n`:** unpaired electrons of the +3 ion. f^n has n unpaired up to the half-filled shell and 14 − n beyond it: La 0, Ce 1, Pr 2, Nd 3, Pm 4, Sm 5, Eu 6, Gd 7, Tb 6, Dy 5, Ho 4, Er 3, Tm 2, Yb 1, Lu 0, Ac 0, Th 1, Pa 2, U 3, Np 4, Pu 5, Am 6, Cm 7, Cf 5.

## Step 3: QC

Rules: `scf_fail`, `not_converged`, `cn_change` (CN differs from the start structure), `short_bond` (metal–donor < 2.0 Å), `water_left` (water O > 3.2 Å from the metal), `missing` (no output). CN counts N/O within 3.2 Å and S within 3.3 Å, as in `xtb-complex-1/pd`. Per-combination numbers are in `structures.csv`.

| Result | Combinations |
|---|---|
| clean | 1381 |
| CN change only | 49 |
| CN change + water left | 18 |
| SCF failure (no output) | 6 |
| not converged only | 3 |
| not converged + CN change | 2 |

- **No `short_bond` cases remain** (there were 9 in the first pass, all Tm).
- **SCF failures (6):** 73_Nd, 79_Nd, 283_Nd, 75_Ce, 204_Am, 232_Am. All say "SCF not converged in 250 cycles". They are spread over 3 metals and 6 different ligands, so there is no pattern.
- **Clean rate by metal:** 0.85 (Tb) to 1.00, most above 0.93. The lowest are Tb 0.85, Ce 0.91, Am 0.93, Tm 0.93, Dy 0.93.
- **Water loss:** at most 12% of complexes for Tb, below 5% for every other metal.

**Spin check.** The Mulliken spin on the metal now matches the requested value everywhere, to within 0.05 (for example Tb 6.02 for uhf 6, Yb 0.99 for uhf 1, Eu 6.02, Gd 7.02, Am 6.01, Cm 6.90).

## Problems found and fixed

### 1. The first pass used f-electron counts as uhf

For f^n ions with n > 7 the number of unpaired electrons is 14 − n, not n. The first pass used n (Tb 8, Dy 9, Ho 10, Er 11, Tm 12, Yb 13, Cf 9), which forced more unpaired electrons than the f shell can hold, so the surplus went onto the ligands. The metal spin showed it directly: Tb came out 7.08 instead of 8, Yb 1.45 instead of 13.

Consequences, and what the corrected rerun gave:

| | First pass | After the fix |
|---|---|---|
| clean rate, the 7 affected metals | 0.35–0.67 | 0.934 |
| complexes losing water (Er, Tm, Ho, Yb, Tb) | 40–55% | ≤ 12% |
| short bonds | 9 (all Tm) | 0 |

### 2. `--chrg 3` was wrong for the 8 charged ligands

36 combinations use ligands that carry a charge, so their complex is +2, 0 or −1. Running them at +3 adds or removes electrons. Their clean rate was 12–25%; with the correct charge it is 0.78.

### Rerun

438 combinations were rerun with both fixes (411 late-f, 36 charged ligands, 9 in both). The superseded runs were kept in `runs_superseded/`.

## Files

- `build_list.py` (env-pd): builds `list.csv`, `tasks.txt`, `unmatched.csv` and `starts/`.
- `run_gxtb.sh` (env-gxtb): runs g-xTB; resumable; `runs/status.txt` logs uhf, charge, rc and wall time.
- `qc.py` (env-pd): writes `structures.csv` — per combination: ligand idx, SMILES, element, uhf, convergence, wall time, CN (start and final), mean M–O and M–N, donor and water distance ranges, metal charge and spin, and QC flags.
- `runs/<idx>_<El>/`: `in.xyz`, `xtbopt.xyz`, `xtb.log` and the other xtb outputs.
- `runs_superseded/<idx>_<El>/`: the first-pass runs of the 438 rerun combinations (wrong uhf or charge).
