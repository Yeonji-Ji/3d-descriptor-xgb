In the JACS paper, newly generated molecules were scored with an ML model trained on ~8,000 experimental rows covering 295 unique ligands. Each row had 1,860 features: 2D ECFP fingerprints, RDKit descriptors, and experimental conditions and metal properties. These 2D features are identical for a given ligand regardless of the metal, and they carry no information about how the ligand binds. 

Here, this experiment tests whether 3D descriptors add more information: 
- g-xTB structures of each metal–ligand complex
- summarized with persistent homology (PD summary features and persistence images). 
- Scope: 1,471 complexes, 7,053 of the 8,075 rows (+3 metals and Th(IV); +5/+6 actinyls and reduced Pu/Np(IV) excluded). 

**Result:** the 3D descriptors did not improve on the 2D model. Their information largely overlaps with the 2D features.
# Data Preparation
## xTB
- Structures were first built with GFN2-xTB (ALPB water, 50 seeds, lowest energy kept). 
	- GFN2-xTB covers the lanthanides but not the actinides (Z > 86), so every complex was built as Eu. 
- To give each metal its own structure, all complexes were rebuilt with g-xTB, gas phase (its implicit solvent models failed for these complexes). g-xTB reproduced the DFT structure of Am(TEDGA)3 (core RMSD 0.11–0.25 Å). 
- By oxidation state: 
	- +3 (1459 combinations): the lowest-energy GFN2 conformer, with the metal swapped, re-optimised. 1431 kept. 
		- Dropped: 11 failed runs, 11 where a design donor left the metal, 6 where the coordination number rose. 
	- +4 (55): 20 new seeds per combination; 10 more if the best seed was in 15–20. 
		- Th was clean, but most Pu and some Np came out effectively reduced to +3 (metal spin ~+1 too high, M–O 0.09 Å longer), so 15 were dropped. 40 kept. 
	- +5/+6 (53, actinyl): O=M=O bent below 170° in 44 of 53 and M=O was too short.
		- Constraining the angle to 180° did not fix the M=O distance or the spin and pushed out equatorial ligands. All dropped. 
- Final dataset: 1471 structures, 7053 of 8075 rows (87%).

## Persistent Homology
### Persistent Diagram:
-  Vietoris-Rips with ripser maxdim=1 on the pairwise distance matrix.
- Eu-specific names were generalised for 20 metals: `EuNO` → `MNO`, `EuNOS` → `MNOS`, `_Eu_` → `_M_`.
#### Used PD Summary Features:
- Computed with `subset_features()` from `xtb-prod/xtb-complex-1/pd/compute_pd.py`, unchanged (same bin edges). 
- Per subset (`heavy`, `MNO`, `MNOS`): 
	- Bins: cov < 1.8 Å, coord 1.8–3.0 Å, outer 3.0–8.0 Å, beyond8 > 8 Å.

| Group                  | Columns                                                          |
| ---------------------- | ---------------------------------------------------------------- |
| Counts by distance bin | `*_H0_cov/coord/outer/beyond8`, `*_H1_cov/coord/outer`           |
| Lifetimes              | `*_H0_totlife/meanlife/maxlife`, `*_H1_totlife/meanlife/maxlife` |
| Ratios                 | `*_H1_per_atom`, `*_H1_frac_cov/coord/outer`                     |
| Metal (MNO, MNOS)      | `*_M_mst_degree`, `*_H0_coord_M`, `*_M_mst_minlen`               |
| Per element (MNOS)     | `*_H0_coord_elem`, `*_H0_outer_elem`, `*_H0_coord_M_elem`        |

### PI features
- Persistence images flattened pixel by pixel, one column per pixel. Six images per structure (`heavy`, `MNO`, `MNOS` × H0, H1): H0 as 40 bins (persistence only), H1 as 20 × 20. Pixels that are zero in every structure were dropped: 1320 → 829 columns, named `PI_<subset>_<dim>_p<pixel>`.

# XGBoost
## Datasets: six combination comparisons

| File                              | Combinations                           | Number of features |
| --------------------------------- | -------------------------------------- | ------------------ |
| `fs1_2D_rdkit_cond.csv`           | 2D + rdkit + cond + metal              | 1860               |
| `fs2_2D_rdkit_cond_denticity.csv` | fs1 + denticity                        | 1861               |
| `fs3_PDsummary_rdkit_cond.csv`    | PD summary + rdkit + cond + metal      | 258                |
| `fs4_PI_rdkit_cond.csv`           | PI + rdkit + cond + metal              | 1018               |
| `fs5_2D_PDsummary_rdkit_cond.csv` | 2D + PD summary + rdkit + cond + metal | 1929               |
| `fs6_2D_PI_rdkit_cond.csv`        | 2D + PI + rdkit + cond + metal         | 2689               |
- 2D = ECFP fingerprint (`CircularFP_*`, 1671 bits).
- fs2: Denticity was added as a control to separate the contribution of PD topology from that of denticity alone.
## Training
### Target
- `Class_index` (3 classes: aqueous / unselective / organic). 
- Excluded from features: `Log_D`, `Distribution_ratio` (derived from the target), `SMILES`, `SMILES_class`. 
- No feature scaling (tree model).

### Feature sets
- All sets share RDKit descriptors (161), experimental conditions (17) and metal properties (11) = 189 columns.
- Six datasets, see above.

### Data split

|          | Ligands | Rows |
| -------- | ------- | ---- |
| trainVal | 267     | 6604 |
| test     | 15      | 449  |

- CV: 3 ligand-level *StratifiedShuffleSplit* folds on trainVal (test_size 0.15, stratified by SMILES_class), ~41 validation ligands per fold. 
	- Validation rows per fold: 946, 835, 567. Folds can share validation ligands.

### Model and Tuning
- XGBClassifier (`objective='multi:softmax'`, `tree_method='hist'`).
- RandomizedSearchCV, 60 candidates per set, scoring = accuracy, same CV folds and `random_state` for every set.
- Search space: `n_estimators` 200/300/500, `learning_rate` 0.1/0.2/0.3, `max_depth` 6/9/12/15, `min_child_weight` 1/3/5/10, `colsample_bytree` 0.4/0.6/0.8, `alpha` 0.01/0.1/0.5, `lambda` 5/10/20.
- The best setting was refit on all trainVal and scored once on the test set.

## Results

All sets also include RDKit descriptors, experimental conditions and metal properties.

| Feature set          | n_features | CV acc (mean ± sd) | CV train acc | Test acc | Test F1 (macro) | Importance of added block |
| -------------------- | ---------- | ------------------ | ------------ | -------- | --------------- | ------------------------- |
| fs1: 2D              | 1860       | 0.660 ± 0.030      | 0.970        | 0.711    | 0.691           | –                         |
| fs2: 2D + denticity  | 1861       | 0.665 ± 0.026      | 0.953        | 0.728    | 0.719           | 0.001                     |
| fs3: PD summary      | 258        | 0.646 ± 0.031      | 0.959        | 0.684    | 0.686           | 0.141                     |
| fs4: PI              | 1018       | 0.650 ± 0.040      | 0.973        | 0.682    | 0.648           | 0.366                     |
| fs5: 2D + PD summary | 1929       | 0.642 ± 0.051      | 0.973        | 0.695    | 0.673           | 0.029                     |
| fs6: 2D + PI         | 2689       | 0.645 ± 0.055      | 0.970        | 0.659    | 0.643           | 0.165                     |
### What can we read from the results:
1. CV accuracy was 0.642–0.665 for all six sets, with fold-to-fold sd of 0.03–0.055, which can be said that no set is distinguishable from 2D.
2. The 2D + denticity (fs2) scored highest (acc 0.728, F1 0.719) for the test acc, but denticity had near-zero importance (0.001) and the test set holds only 15 ligands, so this is not attributed to denticity.
3. PD summary alone (fs3, 258 features) reached ~0.68 test accuracy and F1, close to 2D (1860 features) but still lower, suggesting it carries much of the same information.
4. Adding PD summary or PI to 2D did not improve on 2D in CV or on the test set.
5. From the importance, PD summary and PI were used heavily on their own (0.14 and 0.37) but when 2D was present, the values decreased (0.03 and 0.17): their information overlaps with 2D.
6. All sets overfit strongly (train accuracy 0.95–0.97 vs CV ~0.65).

Caveats: 3 shuffle-split folds (validation ligands can overlap), and a 15-ligand test set.

# Conclusion
- Adding 3D descriptors from g-xTB complexes (PD summary, PI) did not improve on the 2D model. Their information largely overlaps with the 2D ECFP fingerprint.
- All six feature sets overfit to the training ligands and stalled at a similar level of CV accuracy, which suggests that the number of ligands (267 in trainVal) might be the main limit or bottleneck, not the choice of descriptor.