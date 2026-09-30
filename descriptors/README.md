# PD and PI descriptors for the complex dataset

Persistence diagrams and persistence images for all 1471 structures in `complex-dataset/`, attached to the 7053 experimental rows. `complex-dataset/` and `xtb-complex-1/` were only read.

## Headline

- **1471 of 1471 structures** produced diagrams and images; none failed.
- **829 PI columns** after dropping pixels that are zero for every structure (1320 before).
- **`rows_with_pi.csv`: 7053 rows × 2696 columns** — every original column (196 experimental + 1671 CircularFP) plus the 829 PI columns, with the original split kept (trainVal 6604, test 449).
- Runtime: 3 s for the diagrams, 6 s for the images.

## Step 1: generalising the metal

`compute_pd.py` imports `xtb-complex-1/pd/compute_pd.py` rather than copying it, so the maths and every numeric setting are unchanged: Vietoris–Rips with ripser `maxdim=1` on the pairwise distance matrix, H0 deaths cross-checked against the scipy MST edges.

Only the subset definition changed. The metal is detected per structure as the single element outside the organic set (H, C, N, O, S, P, Cl, F, Br, I, Si, B), and the subsets are renamed:

| New | Old | Contents |
|---|---|---|
| heavy | heavy | all non-H atoms (unchanged) |
| MNO | EuNO | metal + N + O |
| MNOS | EuNOS | metal + N + O + S |

This covers all 23 metals in the dataset. The element-specific cutoffs (N/O 3.2 Å, S 3.3 Å) live in the imported module and were not touched; they affect the summary features, not the diagrams.

**Equivalence check** (`eu_check.csv`): the generalised code was run on three Eu structures from the original work — idx 40 (G3, seed 11), idx 187 (G4, seed 15) and idx 34 (G7, seed 40) — and compared with the diagrams stored in `xtb-complex-1/pd/results/raw/`.

- Atom selection identical for all three subsets, in all three structures.
- H0 and H1 identical in shape, and the **largest absolute difference in any diagram is 0.000e+00**.
- The MST pair element tables match as well.

The check runs automatically before the full pass (`--check` runs it alone), and the script stops if it fails.

## Step 2: diagrams

`pd/<structure_id>.npz`, one per structure, holding per subset: atom indices, elements, coordinates, H0, H1, and the MST pair table (local indices, atom indices, element pairs, edge lengths) that the PI step needs.

| | heavy | MNO | MNOS |
|---|---|---|---|
| atoms per structure (median, max) | 41, 202 | 12, 34 | 12, 34 |
| H1 loops per structure (median, max) | 6, 54 | 0, 8 | 0, 8 |

Most MNO/MNOS subsets have no loops at all, so their H1 images are empty for many structures.

## Step 3: images

Kernel settings taken unchanged from `xtb-complex-1/pi/config.json`:

- **H0:** σ = 0.10 × (1 + 0.5 · |ΔEN|) per MST edge, weight 1, the infinite bar skipped.
- **H1:** σ = 0.10, weight = min(persistence / 0.5 Å, 1).
- Pixel-integrated Gaussians truncated at 3σ, no per-image normalisation.
- The Pauling electronegativity table was extended from the organic elements to all 23 metals (La 1.10 … Cf 1.30; Eu stays 1.20 as before).

**Grids were fitted once over the diagrams of all 1471 structures together**, per subset and homology dimension, then applied to every structure — never per metal, per oxidation state or per split. Ranges are padded by 3σ so no mass falls off the edge.

| Subset / dim | Persistence range (Å) | Birth range (Å) | Resolution | Points used | Pixels kept |
|---|---|---|---|---|---|
| heavy H0 | 0.544 – 3.678 | all 0 (1-D) | 40 | 62053 | 35 |
| heavy H1 | 0.0 – 4.102 | 1.070 – 7.235 | 20 × 20 | 10143 | 169 |
| MNO H0 | 0.642 – 6.178 | all 0 (1-D) | 40 | 16645 | 36 |
| MNO H1 | 0.0 – 2.096 | 2.118 – 6.654 | 20 × 20 | 592 | 260 |
| MNOS H0 | 0.642 – 6.178 | all 0 (1-D) | 40 | 16755 | 36 |
| MNOS H1 | 0.0 – 2.096 | 2.118 – 6.654 | 20 × 20 | 650 | 293 |

Two deliberate choices, both different from the single combined image used in `xtb-complex-1/pi`:

1. **H0 images are one-dimensional.** Every H0 bar of a Vietoris–Rips filtration is born at 0; this is asserted for all 1471 structures. A birth axis would therefore carry no information, so H0 uses persistence only, at 40 bins.
2. **H1 uses 20 × 20 instead of 40 × 40.** Splitting by dimension turns one image per subset into two, so six images in total. At 40 × 40 that would be 9600 columns; 20 × 20 keeps the total at 1320 before pruning (829 after), close to the 1600 of the earlier single-image design.

Per-image totals: heavy H0 carries ~42 units of mass per structure (one per MST edge), heavy H1 ~4.2, and the MNO/MNOS H1 images only ~0.1–0.14 because most of those subsets have no loops.

## Step 4: attaching to the rows

`rows_with_pi.csv` is rebuilt from `trainVal_dataset.csv` and `test_dataset.csv` (not from `rows.csv`), so every original column is present, and then joined to the PI table on `structure_id`.

Checks that ran as part of the build:

- Row count, `structure_id` distribution and split distribution all match `complex-dataset/rows.csv`.
- Every row received a PI vector (no missing values).
- The original columns are byte-identical to `complex-dataset/rows.csv` — PI values are rounded to 6 significant digits, but no original column is reformatted.
- `CircularFP_*` values match the source datasets for all 282 ligands.

## Failures

None. All 1471 structures produced diagrams and images, and all 7053 rows were matched.

## Files

- `compute_pd.py` — Step 1–2; `--check` runs only the Eu equivalence check. Writes `pd/`, `pd_log.csv`, `eu_check.csv`.
- `build_pi.py` — Step 3–4. Writes `pi_config.json`, `pi_features.csv`, `rows_with_pi.csv`.
- `pd/<structure_id>.npz` — 1471 diagrams (17 MB total).
- `pd_log.csv` — per structure: metal, atom counts and H1 counts per subset, any error.
- `eu_check.csv` — the three-structure equivalence check.
- `pi_config.json` — kernel settings, the electronegativity table, and the fitted grids.
- `pi_features.csv` — 1471 × 829 PI columns, named `PI_<subset>_<dim>_p<pixel>` (5.4 MB).
- `rows_with_pi.csv` — the model input, 7053 × 2696 (61 MB).
