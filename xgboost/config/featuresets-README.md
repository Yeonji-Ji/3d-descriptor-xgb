# 3D descriptor 실험용 feature set (2D vs PD/PI)

g-xTB 구조에서 뽑은 3D descriptor(PD summary, PI)가 기존 2D feature 대비 도움이 되는지 확인하려고, 같은 7,053개 row를 6가지 feature 조합으로 나눠 놓은 CSV들입니다.

## 데이터

- **Row 7,053개** — `dataset/complex-dataset-pd/rows_with_pi.csv` (g-xTB 구조가 붙은 실험 row 전부)
- **구조 1,471개 / ligand 282개**, 산화수 +3 6,461 row + +4 592 row
- **Split**: trainVal 6,604 / test 449 (원본 그대로, 재분할 안 함)
- **Target**: `Class_index` (0 aqueous 1,689 / 1 unselective 1,521 / 2 organic 3,843)
- 6개 파일 모두 **행 순서와 메타 컬럼이 동일**하므로 fold를 한 번 만들어 모든 세트에 그대로 쓸 수 있습니다.

## 6개 feature set

`rdkit + cond + metal`은 모든 세트에 공통으로 들어갑니다. 세트 간 차이는 2D / denticity / PD summary / PI 뿐입니다.

| 파일 | 조합 | feature 수 | 크기 |
|---|---|---|---|
| `fs1_2D_rdkit_cond.csv` | 2D + rdkit + cond + metal | 1860 | 35 MB |
| `fs2_2D_rdkit_cond_denticity.csv` | fs1 + denticity | 1861 | 35 MB |
| `fs3_PDsummary_rdkit_cond.csv` | PD summary + rdkit + cond + metal | 258 | 16 MB |
| `fs4_PI_rdkit_cond.csv` | PI + rdkit + cond + metal | 1018 | 39 MB |
| `fs5_2D_PDsummary_rdkit_cond.csv` | 2D + PD summary + rdkit + cond + metal | 1929 | 39 MB |
| `fs6_2D_PI_rdkit_cond.csv` | 2D + PI + rdkit + cond + metal | 2689 | 62 MB |

## Block 구성

| Block | 컬럼 수 | 내용 |
|---|---|---|
| `2D` | 1671 | `CircularFP_*` |
| `rdkit` | 161 | RDKit descriptor (`MaxAbsEStateIndex` ~ `fr_unbrch_alkane`) |
| `cond` | 17 | 실험 조건 — extractant/acid/metal 농도, 온도, diluent A/B 물성 |
| `metal` | 11 | 금속 물성 (atomic number, radii, EN, IE 등) + **`Oxidation_state`** |
| `dent` | 1 | `denticity` |
| `PDsummary` | 69 | heavy / MNO / MNOS subset의 topology 요약 |
| `PI` | 829 | `PI_<subset>_<dim>_p<pixel>` |

- **`Oxidation_state`를 feature로 넣었습니다.** 이전 +3 전용 실험과 달리 여기는 +3과 +4가 섞여 있어 상수가 아닙니다.
- **metal block이 모든 세트에 들어갑니다.** 금속이 20종이라 금속 물성 없이는 비교가 무의미합니다.

## PD summary 컬럼

`xtb-prod/xtb-complex-1/pd/compute_pd.py`의 `subset_features()`를 그대로 호출해 계산했습니다 (bin 경계 등 수치 설정 동일). subset별로:

- 개수: `*_H0_cov/coord/outer/beyond8`, `*_H1_*` (bin: cov < 1.8 Å, coord 1.8–3.0 Å, outer 3.0–8.0 Å)
- 수명: `*_H0_totlife/meanlife/maxlife`, `*_H1_totlife/meanlife/maxlife`
- 비율: `*_H1_per_atom`, `*_H1_frac_cov/coord/outer`
- 금속 subset 추가: `*_M_mst_degree`, `*_H0_coord_M`, `*_M_mst_minlen`, MNOS는 `*_H0_coord_elem`, `*_H0_outer_elem`, `*_H0_coord_M_elem`

원래 Eu 전용 이름(`EuNO`, `EuNOS`, `_Eu_`)은 금속이 20종이므로 `MNO`, `MNOS`, `_M_`로 바꿨습니다.

## 메타 컬럼 (feature 아님)

모든 파일 앞쪽에 동일하게 들어갑니다:

`structure_id`, `split`, `SMILES`, `SMILES_class`, `element`, `Oxidation_state_meta`, `idx`, `family`, `Class_index`, `Distribution_ratio`, `Log_D`

- **`Distribution_ratio`와 `Log_D`는 target이므로 feature로 쓰면 안 됩니다.** 회귀 실험용으로 남겨둔 것뿐입니다.
- `Oxidation_state_meta`는 읽기용 복사본이고, 실제 feature는 `metal` block 안의 `Oxidation_state`입니다.

## `unique_smiles.csv` — train/val 분할용 ligand 목록

`orig-data/unique_smiles.csv`(295개)를 지금 데이터에 실제로 들어있는 **282개 ligand만** 남겨 자른 파일입니다. 컬럼 구성은 원본과 동일하고 참고용 3개를 덧붙였습니다.

| 컬럼 | 내용 |
|---|---|
| `SMILES`, `SMILES_class`, `split` | 원본 그대로 |
| `n_rows` | 이 ligand가 들어간 실험 row 수 |
| `n_structures` | 이 ligand의 g-xTB 구조 수 (금속·산화수 조합 수) |
| `n_metals` | 이 ligand와 짝지어진 금속 종류 수 |

확인한 내용:

- 282개 ligand 전부 원본 `unique_smiles.csv`에 있었습니다 (빠진 것 없음). 원본에서 제외된 13개는 지금 데이터에 구조가 없는 ligand입니다.
- **split 라벨이 row 데이터와 100% 일치**하고, 한 ligand가 두 split에 걸친 경우는 없습니다. trainVal 267 / test 15.
- `SMILES_class`는 14종이고 가장 작은 class도 ligand 4개라, 이 컬럼으로 stratify 해서 train/val을 나눌 수 있습니다 (예전 notebook의 `Class_by_SMARTS_combineRare` 자리에 그대로 쓰면 됩니다).

trainVal만 필요하면 `split == "trainVal"`로 걸러 쓰면 되고, `orig-data/unique_smiles_EuAm_label.csv`(Eu/Am 라벨 포함)도 같은 방식으로 자를 수 있습니다 — 지금은 금속이 20종이라 그 라벨의 의미가 줄어들어 넣지 않았습니다.

## 학습 시 주의

1. **Split은 SMILES 기준입니다.** 같은 ligand가 여러 row에 나오므로, CV를 돌릴 때도 `StratifiedGroupKFold(groups=SMILES)`처럼 ligand 단위로 묶어야 합니다. row 단위로 섞으면 성능이 크게 부풀려집니다.
2. **PI는 항상 0인 pixel을 이미 제거했습니다** (1320 → 829). fold 안에서 추가로 거의 안 쓰이는 pixel을 잘라내면 더 줄일 수 있습니다.
3. **데이터의 한계** — +5/+6 (actinyl) 구조는 전부 빠졌고, Pu(IV)는 row의 27%만 남아 있습니다. 자세한 내용은 `dataset/complex-dataset/README.md` 참고.

## 파일

- `fs1` ~ `fs6` CSV 6개
- `unique_smiles.csv` — train/val 분할용 ligand 282개
- `feature_blocks.json` — block별 컬럼 목록과 세트 정의
- `featureset_summary.csv` — 위 표와 동일한 내용
- `build_featuresets.py` — 재생성 스크립트 (약 17초)
