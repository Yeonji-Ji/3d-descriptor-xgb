# g-xTB build of the +4 metal complexes

New conformers for the Oxidation_state == 4 rows, built and optimised entirely here — there are no GFN2 structures for these metals. `xtb-complex-1/` was only read (code reused, nothing written there).

## Headline

- **55 combinations** (37 ligands × Th/Np/Pu), **all 55 have a converged structure** in `best/`.
- **Geometry QC is clean everywhere:** CN 9 in every complex, no water left the metal, no metal–donor bond below 2.0 Å, no CN change from the start structure.
- **24 of 55 carry no flag at all; 33 pass every geometry and spin check** (the other 9 differ only by a late best seed). The remaining flags are metal spin (22) and late best seed (13); both are discussed below.
- **Cost:** 1260 seed runs. Summed per-run wall time 80.9 h (median 109 s, mean 231 s); xtb's own timers give about 67 CPU-h, the difference being periods when the machine was idle. Roughly 9 h of wall time on 9 workers.

## Step 1: the list

- **Rows:** 728 rows with Oxidation_state == 4 (Th 453, Pu 175, Np 100), 37 ligands.
- **Combinations:** 55 unique (SMILES, atomic number) pairs — Th 28, Pu 20, Np 7. Every ligand was already in the `xtb-complex-1` donor table, so all donor assignments were reused and nothing had to be matched by hand.
- **Ligands:** all neutral, so every complex is **+4**. 9 families, water counts 0–8, complexes 62–238 atoms.

## Step 2: conformers

`build_list.py` reuses `xtb-complex-1/src/build.py` — same donor table, same ETKDGv3 settings (`useRandomCoords`, `maxIterations` 2000, no force field), same water placement. The builder hardcodes Eu, so the metal is patched per element and the metal's formal charge set to +4. Waters follow the +3 set: n_water = 9 − denticity.

- **20 seeds per combination** (seeds 1–20): 1100 start structures, all embedded by ETKDG.
- **The geometric water shell fallback was never needed** (0 seeds), unlike the +3 set.

## Step 3: runs

`xtb in.xyz --gxtb --opt --chrg 4 --uhf <n>`, gas phase, 1 thread each, 9 in parallel (`run_gxtb.sh`). uhf = unpaired f electrons of the +4 ion: Th 0, Np 3, Pu 4.

- **1260 seed runs** (1100 + 160 extra), **1141 converged (90.6%)**.
- **Failures:** 119 seeds, of which 106 are "SCF not converged in 250 cycles". By metal: Pu 79, Th 31, Np 9. They are individual seeds, not whole combinations — the worst case, 12_Pu, still has 8 converged seeds of 30.

## Step 4: seed check

For every combination the lowest converged energy was taken (`seeds.csv` has every seed's energy).

- **16 combinations had their best seed in 15–20** and got 10 extra seeds (21–30).
- **3 of those improved:** 174_Np (best seed now 25), 174_Th (25) and 249_Th (26).
- **13 still peak at seed 15–20.** The rule allows only one extra round, so they were kept as they are; they carry the flag `late_best_seed`. Their conformer search may not be converged.
- **Energy gap to the second-best seed:** median 1.75 kcal/mol, mean 2.52, range 0.01–13.5. In 21 of 55 the gap is below 1 kcal/mol, so the best and second-best conformers are nearly degenerate.

## Step 5: QC

Same rules as `xtb-complex-gxtb/qc.py` (imported directly), plus a spin check.

| Flag | Combinations |
|---|---|
| none | 24 |
| spin_mismatch | 18 |
| late_best_seed | 9 |
| spin_mismatch + late_best_seed | 4 |
| cn_change, short_bond, water_left, missing | 0 |

Counting overlaps: 22 spin_mismatch, 13 late_best_seed. 24 combinations have no flag at all, and 33 have nothing beyond `late_best_seed`.

**Geometry** (means over combinations):

| Metal | Combinations | CN | Mean M–O | Mean M–N | Metal charge |
|---|---|---|---|---|---|
| Th | 28 | 9.0 | 2.419 | 2.502 | +2.37 |
| Np | 7 | 9.0 | 2.428 | – | +2.52 |
| Pu | 20 | 9.0 | 2.437 | 2.565 | +2.44 |

**Spin.** The flag is `|metal spin − uhf| > 0.2`:

- **Th (uhf 0):** exactly 0.000 everywhere.
- **Np (uhf 3):** metal spin 3.15–4.04, so 0.15–1.04 above the requested value.
- **Pu (uhf 4):** metal spin 4.22–5.06, i.e. 0.22–1.06 above, in two clusters (13 combinations near +1.0, 7 near +0.25).

The *total* unpaired count is fixed by `--uhf`, so a metal spin above it means the ligands carry compensating negative spin density — spin polarisation, not a different total spin state. It is not obviously an error, but the near-integer excess of +1 for most Pu (and 2 Np) complexes is worth checking before these structures are used, since it suggests one electron's worth of spin density shifted onto the ligands.

**Metals that failed systematically:** none. Th, Np and Pu all produced a structure for every combination, with identical CN and similar M–O. Pu has the most failed individual seeds (79 of its 600), but still 18 converged seeds per combination on average.

## Files

- `build_list.py` (env-gxtb): builds `list.csv`, `tasks.txt`, `build_log.csv` and `starts/<idx>_<El>_s<seed>.xyz`. `--seeds 21 30 --only <names>` appends an extra round.
- `run_gxtb.sh` (env-gxtb): runs g-xTB; resumable (skips converged runs); `runs/status.txt` logs uhf, charge, rc and wall time.
- `analyze.py` (env-pd): writes `seeds.csv`, `structures.csv` and `best/`; reuses the QC code from `xtb-complex-gxtb/qc.py`.
- `runs/<idx>_<El>_s<seed>/`: `in.xyz`, `xtbopt.xyz`, `xtb.log`.
- `best/<idx>_<El>.xyz`: the lowest-energy converged structure per combination.
- `seeds.csv`: per seed — convergence, energy, gradient, metal charge and spin, wall time.
- `structures.csv`: per combination — best seed, energy, gap to second best, CN before and after, M–O and M–N, water range, metal charge and spin, spin excess, flags.
