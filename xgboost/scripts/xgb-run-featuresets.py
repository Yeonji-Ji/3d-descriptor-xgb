# ============================================================
# XGBoost on the six feature sets (dataset/ml-featuresets-3d-vs-2d)
#
#   python xgb-run-featuresets.py                 # run fs1..fs6, then compare
#   python xgb-run-featuresets.py --fs fs3 fs4    # run only some sets, then compare
#   python xgb-run-featuresets.py --compare-only  # rebuild comparison from saved results
#
# Output: <RUN>/results/<fs name>/  run.log, cv_results_top100.csv, top5_models.csv,
#                             best_params.json, confusion_matrix.csv,
#                             classification_report.csv, feature_importance.csv,
#                             best_model.json, summary.json
#         <RUN>/results/comparison.csv   one row per feature set
# ============================================================

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, classification_report,
                             confusion_matrix, f1_score)
from sklearn.model_selection import StratifiedShuffleSplit, RandomizedSearchCV
import xgboost as xgb

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[1] / "dataset" / "ml-featuresets-3d-vs-2d"
RUN = "run2"                                  # results go to <RUN>/results; change per run
OUT = HERE / RUN / "results"
BLOCKS = json.load(open(DATA / "feature_blocks.json"))
FEATURESETS = sorted(BLOCKS["sets"])          # fs1_... fs6_...
META_COLS = BLOCKS["metadata_columns"]        # SMILES, split, Class_index, Log_D, ... never features
TARGET = "Class_index"                        # 0 aqueous / 1 unselective / 2 organic

param_dist = {
    'n_estimators': [200, 300, 500],
    'learning_rate': [0.1, 0.2, 0.3],
    'max_depth': [6, 9, 12, 15],
    'min_child_weight': [1, 3, 5, 10],
    'colsample_bytree': [0.4, 0.6, 0.8],
    'alpha': [0.01, 0.1, 0.5],
    'lambda': [5, 10, 20],
}


def make_log(fs_dir):
    log = logging.getLogger(fs_dir.name)
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s  %(message)s", "%H:%M:%S")
    for h in (logging.FileHandler(fs_dir / "run.log", mode="w"), logging.StreamHandler()):
        h.setFormatter(fmt)
        log.addHandler(h)
    return log


def ligand_folds(train_val_df):
    """Stratified ligand-level splits (by SMILES_class), as positional row indices.
    Row order is identical in all six CSVs, so every feature set gets the same folds."""
    smiles = pd.read_csv(DATA / "unique_smiles.csv")
    smiles = smiles[smiles["split"] == "trainVal"].reset_index(drop=True)
    splitter = StratifiedShuffleSplit(n_splits=3, test_size=0.15, random_state=14)
    folds = []
    for tr, va in splitter.split(smiles, smiles["SMILES_class"]):
        in_val = train_val_df["SMILES"].isin(smiles.loc[va, "SMILES"]).to_numpy()
        folds.append((np.where(~in_val)[0], np.where(in_val)[0]))
    return folds


def run(fs):
    fs_dir = OUT / fs
    fs_dir.mkdir(parents=True, exist_ok=True)
    log = make_log(fs_dir)
    t0 = time.time()

    df = pd.read_csv(DATA / f"{fs}.csv")
    train_val_df = df[df["split"] == "trainVal"].reset_index(drop=True)
    test_df = df[df["split"] == "test"].reset_index(drop=True)
    feature_cols = [c for c in df.columns if c not in META_COLS]

    X_trainVal, y_trainVal = train_val_df[feature_cols], train_val_df[TARGET]
    X_test, y_test = test_df[feature_cols], test_df[TARGET]

    folds = ligand_folds(train_val_df)
    log.info(f"{fs}: {len(feature_cols)} features | trainVal {X_trainVal.shape} | test {X_test.shape}")
    log.info(f"class distribution (trainVal): {y_trainVal.value_counts().sort_index().to_dict()}")
    for i, (tr, va) in enumerate(folds):
        log.info(f"fold {i}: train rows {len(tr)}, val rows {len(va)}, "
                 f"val ligands {train_val_df.loc[va, 'SMILES'].nunique()}")

    xgb_clf = xgb.XGBClassifier(
        objective='multi:softmax', num_class=3, eval_metric='merror',
        tree_method='hist',
        n_jobs=1,
    )
    search = RandomizedSearchCV(
        estimator=xgb_clf,
        param_distributions=param_dist,
        n_iter=60,
        scoring='accuracy',
        cv=folds,
        refit=True,
        return_train_score=True,
        random_state=0,
        n_jobs=-1
    )
    search.fit(X_trainVal, y_trainVal)

    # CV results
    cv = pd.DataFrame(search.cv_results_).sort_values("rank_test_score")
    cv.head(100).to_csv(fs_dir / "cv_results_top100.csv", index=False)
    top5 = cv.head(5)[["rank_test_score", "params", "mean_train_score",
                       "mean_test_score", "std_test_score", "mean_fit_time"]]
    top5.to_csv(fs_dir / "top5_models.csv", index=False)
    json.dump(search.best_params_, open(fs_dir / "best_params.json", "w"), indent=1)
    log.info(f"best params: {search.best_params_}")
    log.info(f"best CV val accuracy: {search.best_score_:.4f}")

    # Test evaluation
    best_xgb = search.best_estimator_
    best_xgb.save_model(fs_dir / "best_model.json")
    y_pred = best_xgb.predict(X_test)
    labels = [0, 1, 2]
    pd.DataFrame(confusion_matrix(y_test, y_pred, labels=labels),
                 index=[f"Actual {c}" for c in labels],
                 columns=[f"Predicted {c}" for c in labels]).to_csv(fs_dir / "confusion_matrix.csv")
    pd.DataFrame(classification_report(y_test, y_pred, labels=labels, output_dict=True, zero_division=0)
                 ).transpose().to_csv(fs_dir / "classification_report.csv")

    # Feature importance (gain), tagged with its block so 2D vs PD/PI contribution is visible
    block_of = {c: b for b, cols in BLOCKS["blocks"].items() for c in cols}
    imp = pd.DataFrame({"feature": feature_cols, "importance": best_xgb.feature_importances_})
    imp["block"] = imp["feature"].map(block_of)
    imp.sort_values("importance", ascending=False).to_csv(fs_dir / "feature_importance.csv", index=False)
    block_share = imp.groupby("block")["importance"].sum().round(4).to_dict()

    best = cv.iloc[0]
    summary = {
        "featureset": fs,
        "blocks": "+".join(BLOCKS["sets"][fs]),
        "n_features": len(feature_cols),
        "cv_val_acc_mean": round(best["mean_test_score"], 4),
        "cv_val_acc_std": round(best["std_test_score"], 4),
        "cv_train_acc_mean": round(best["mean_train_score"], 4),
        "test_acc": round(accuracy_score(y_test, y_pred), 4),
        "test_balanced_acc": round(balanced_accuracy_score(y_test, y_pred), 4),
        "test_f1_macro": round(f1_score(y_test, y_pred, average="macro"), 4),
        **{f"test_f1_class{c}": round(v, 4)
           for c, v in zip(labels, f1_score(y_test, y_pred, labels=labels, average=None))},
        **{f"importance_{b}": v for b, v in block_share.items()},
        "best_params": json.dumps(search.best_params_),
        "runtime_min": round((time.time() - t0) / 60, 1),
    }
    json.dump(summary, open(fs_dir / "summary.json", "w"), indent=1)
    log.info(f"test acc {summary['test_acc']} | balanced {summary['test_balanced_acc']} "
             f"| F1 macro {summary['test_f1_macro']} | {summary['runtime_min']} min")


def compare():
    rows = [json.load(open(p)) for p in sorted(OUT.glob("*/summary.json"))]
    if not rows:
        print("no results yet")
        return
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT / "comparison.csv", index=False)
    print(comp[["featureset", "n_features", "cv_val_acc_mean", "cv_val_acc_std",
                "test_acc", "test_balanced_acc", "test_f1_macro"]].to_string(index=False))
    print(f"\nsaved {OUT / 'comparison.csv'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fs", nargs="+", default=["fs1", "fs2", "fs3", "fs4", "fs5", "fs6"],
                    help="feature set prefixes, e.g. fs1 fs3")
    ap.add_argument("--compare-only", action="store_true")
    args = ap.parse_args()

    if not args.compare_only:
        for prefix in args.fs:
            run(next(f for f in FEATURESETS if f.startswith(prefix + "_")))
    compare()
