# ============================================================
# Fold-by-fold comparison of the six feature sets against fs1 (2D), from one finished run.
# No hyperparameter search: fold accuracy comes from the saved cv_results; macro F1 needs
# predictions, so each set's best_params are refit on the same folds (3 fits per set).
#
#   python fold-compare.py            # uses run2
#   python fold-compare.py --run run1
#
# Output: <run>/fold_compare/  fold_scores.csv, paired_vs_2D.csv, fold_accuracy.png, README.md
# ============================================================

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import accuracy_score, f1_score
import xgboost as xgb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
# reuse data paths, meta columns and the exact fold construction of the training script
_spec = importlib.util.spec_from_file_location("xgbfs", HERE / "xgb-run-featuresets.py")
xgbfs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(xgbfs)

BASE = xgbfs.FEATURESETS[0]   # fs1_2D_rdkit_cond


def short(fs):
    return fs.split("_")[0]


def fold_scores(fs, res_dir):
    """Per-fold val accuracy from cv_results at the best index, plus refit accuracy and macro F1."""
    cv = pd.read_csv(res_dir / fs / "cv_results_top100.csv")
    best = cv.loc[cv["rank_test_score"].idxmin()]
    best_params = json.load(open(res_dir / fs / "best_params.json"))
    assert best_params == eval(best["params"]), f"{fs}: best_params.json != rank-1 row"

    df = pd.read_csv(xgbfs.DATA / f"{fs}.csv")
    tv = df[df["split"] == "trainVal"].reset_index(drop=True)
    X = tv[[c for c in df.columns if c not in xgbfs.META_COLS]]
    y = tv[xgbfs.TARGET]

    rows = []
    for k, (tr, va) in enumerate(xgbfs.ligand_folds(tv)):
        clf = xgb.XGBClassifier(objective='multi:softmax', num_class=3, eval_metric='merror',
                                tree_method='hist', n_jobs=-1, **best_params)
        clf.fit(X.iloc[tr], y.iloc[tr])
        pred = clf.predict(X.iloc[va])
        rows.append({"featureset": fs, "fold": k,
                     "accuracy": best[f"split{k}_test_score"],
                     "accuracy_refit": accuracy_score(y.iloc[va], pred),
                     "f1_macro": f1_score(y.iloc[va], pred, average="macro")})
        print(f"{fs} fold {k}: acc {rows[-1]['accuracy']:.4f} "
              f"(refit {rows[-1]['accuracy_refit']:.4f}), F1 {rows[-1]['f1_macro']:.4f}")
    return rows


def paired(scores, metric):
    """Each set minus fs1, fold by fold."""
    wide = scores.pivot(index="fold", columns="featureset", values=metric)
    n = len(wide)
    out = []
    for fs in wide.columns.drop(BASE):
        d = (wide[fs] - wide[BASE]).to_numpy()
        mean, se = d.mean(), d.std(ddof=1) / np.sqrt(n)
        half = stats.t.ppf(0.975, n - 1) * se
        tied = np.isclose(d, 0)
        out.append({
            "featureset": fs, "metric": metric, "n_folds": n,
            "mean_diff": mean, "ci95_low": mean - half, "ci95_high": mean + half,
            "better": int(((d > 0) & ~tied).sum()),
            "worse": int(((d < 0) & ~tied).sum()),
            "tied": int(tied.sum()),
            "t_test_p": stats.ttest_rel(wide[fs], wide[BASE]).pvalue if not tied.all() else 1.0,
            "wilcoxon_p": stats.wilcoxon(d).pvalue if not tied.all() else 1.0,
        })
    return out


def verdict(row):
    if row["ci95_low"] > 0:
        return "beats 2D"
    if row["ci95_high"] < 0:
        return "loses to 2D"
    return "no detectable difference from 2D"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="run2")
    run_dir = HERE / ap.parse_args().run
    res_dir = run_dir / "results"
    out = run_dir / "fold_compare"
    out.mkdir(exist_ok=True)

    # Step 1: fold scores
    scores = pd.DataFrame([r for fs in xgbfs.FEATURESETS for r in fold_scores(fs, res_dir)])
    scores.round(4).to_csv(out / "fold_scores.csv", index=False)

    comp = pd.read_csv(res_dir / "comparison.csv").set_index("featureset")
    checks = []
    for fs, g in scores.groupby("featureset"):
        mean_saved = comp.loc[fs, "cv_val_acc_mean"]
        refit_gap = (g["accuracy_refit"] - g["accuracy"]).abs().max()
        checks.append(f"- {short(fs)}: fold mean {g['accuracy'].mean():.4f} vs comparison.csv "
                      f"{mean_saved:.4f} -> {'OK' if abs(g['accuracy'].mean() - mean_saved) < 5e-5 else 'MISMATCH'}; "
                      f"refit vs saved fold accuracy max gap {refit_gap:.4f}")
    print("\n".join(checks))

    # Step 2: paired comparison against fs1
    pv = pd.DataFrame(paired(scores, "accuracy") + paired(scores, "f1_macro"))
    pv.round(4).to_csv(out / "paired_vs_2D.csv", index=False)

    # Plot: one line per fold across the six sets
    wide = scores.pivot(index="featureset", columns="fold", values="accuracy")
    fig, ax = plt.subplots(figsize=(7, 4))
    for k in wide.columns:
        ax.plot([short(f) for f in wide.index], wide[k], marker="o", label=f"fold {k}")
    ax.plot([short(f) for f in wide.index], wide.mean(axis=1), "k--", lw=1, label="mean")
    ax.set_ylabel("validation accuracy")
    ax.set_title(f"{run_dir.name}: per-fold accuracy (best params per set)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / "fold_accuracy.png", dpi=150)

    # README
    n = scores["fold"].nunique()
    min_w = 2 / 2 ** n   # smallest two-sided exact Wilcoxon p with n non-zero differences
    t = pv.assign(featureset=pv["featureset"].map(short)).round(4).astype(str)
    table = "\n".join(["| " + " | ".join(t.columns) + " |", "|" + "---|" * t.shape[1]]
                      + ["| " + " | ".join(r) + " |" for r in t.to_numpy()])
    acc = pv[pv["metric"] == "accuracy"].set_index("featureset")
    f1 = pv[pv["metric"] == "f1_macro"].set_index("featureset")
    lines = [f"- **{short(fs)}**: {verdict(acc.loc[fs])} on accuracy "
             f"(mean diff {acc.loc[fs, 'mean_diff']:+.4f}, 95% CI {acc.loc[fs, 'ci95_low']:+.4f} to "
             f"{acc.loc[fs, 'ci95_high']:+.4f}, better/worse/tied {acc.loc[fs, 'better']}/"
             f"{acc.loc[fs, 'worse']}/{acc.loc[fs, 'tied']}); macro F1: {verdict(f1.loc[fs])}."
             for fs in acc.index]
    (out / "README.md").write_text(f"""# {run_dir.name}: feature sets vs fs1 (2D), fold by fold

Each set uses its own best hyperparameters from {run_dir.name} and the same {n} ligand-level folds
(StratifiedShuffleSplit by SMILES_class, random_state=14). Accuracy per fold is read from the saved
cv_results; macro F1 comes from refitting the best params on the same folds (no search).
Difference = set − fs1. "Beats/loses" means the 95% t-interval of the mean difference excludes 0.

## Checks

{chr(10).join(checks)}

## Paired comparison

{table}

## Verdict

{chr(10).join(lines)}

## Caveats

- Only {n} folds. The t-interval has {n - 1} degrees of freedom, so it is wide.
- With {n} folds the exact two-sided Wilcoxon p-value cannot go below {min_w:.4f}, so it can never
  reach 0.05 here.
- The folds come from StratifiedShuffleSplit, so validation ligands can overlap between folds; the
  folds are not independent, and the t-test treats them as if they were.

![per-fold accuracy](fold_accuracy.png)
""")
    print(f"\nsaved to {out}")


if __name__ == "__main__":
    main()
