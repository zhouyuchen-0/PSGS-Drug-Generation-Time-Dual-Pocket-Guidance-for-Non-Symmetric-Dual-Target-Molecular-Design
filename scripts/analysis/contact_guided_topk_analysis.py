#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analyze Top-k enrichment in the contact-guided final-generation set.

Outputs full-set metrics, dock_sum-ranked Top-k metrics, and
integrated-priority-ranked Top-k metrics.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Contact-guided final set Top-k analysis.")
    parser.add_argument("--input_csv", required=True, help="Contact-guided final set CSV.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--topk", nargs="+", type=int, default=[10, 20, 50, 100])
    return parser.parse_args()


def read_contact_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path).copy()
    required = ["smiles", "dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for col in ["dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa", "lip_viol", "is_novel"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "lip_viol" not in df.columns:
        df["lip_viol"] = 0
    if "is_novel" not in df.columns:
        df["is_novel"] = 1

    return df.dropna(subset=["dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa"]).copy()


def dual_hit_rate(df: pd.DataFrame, threshold: float) -> float:
    return float(((df["dock_3fap"] <= threshold) & (df["dock_7pqv"] <= threshold)).mean() * 100.0)


def minmax_score(values: pd.Series, higher_is_better: bool = True) -> pd.Series:
    v = pd.to_numeric(values, errors="coerce")
    v_min = v.min()
    v_max = v.max()
    if np.isclose(v_min, v_max):
        return pd.Series(np.full(len(v), 0.5), index=v.index)
    score = (v - v_min) / (v_max - v_min)
    return score if higher_is_better else 1.0 - score


def add_priority_scores(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["dock_balance_gap"] = (df["dock_3fap"] - df["dock_7pqv"]).abs()

    df["S_dock_sum"] = minmax_score(df["dock_sum"], higher_is_better=False)
    df["S_dock_3fap"] = minmax_score(df["dock_3fap"], higher_is_better=False)
    df["S_dock_7pqv"] = minmax_score(df["dock_7pqv"], higher_is_better=False)
    df["S_balance"] = minmax_score(df["dock_balance_gap"], higher_is_better=False)

    cond_85 = (df["dock_3fap"] <= -8.5) & (df["dock_7pqv"] <= -8.5)
    cond_90 = (df["dock_3fap"] <= -9.0) & (df["dock_7pqv"] <= -9.0)
    cond_95 = (df["dock_3fap"] <= -9.5) & (df["dock_7pqv"] <= -9.5)
    df["S_dual_hit"] = 0.0
    df.loc[cond_85, "S_dual_hit"] = 0.50
    df.loc[cond_90, "S_dual_hit"] = 0.75
    df.loc[cond_95, "S_dual_hit"] = 1.00

    df["S_qed"] = df["qed"].clip(lower=0, upper=1)
    df["S_sa"] = minmax_score(df["sa"], higher_is_better=False)
    df["S_lipinski"] = (pd.to_numeric(df["lip_viol"], errors="coerce").fillna(99) <= 0).astype(float)
    df["S_novelty"] = pd.to_numeric(df["is_novel"], errors="coerce").fillna(1).astype(float)

    df["S_target"] = (
        0.45 * df["S_dock_sum"]
        + 0.20 * df["S_dual_hit"]
        + 0.15 * df["S_dock_3fap"]
        + 0.15 * df["S_dock_7pqv"]
        + 0.05 * df["S_balance"]
    )
    df["S_chem"] = (
        0.45 * df["S_qed"]
        + 0.35 * df["S_sa"]
        + 0.15 * df["S_lipinski"]
        + 0.05 * df["S_novelty"]
    )
    df["S_priority"] = 0.65 * df["S_target"] + 0.35 * df["S_chem"]
    return df


def summarize_full_set(df: pd.DataFrame) -> dict:
    return {
        "n": len(df),
        "unique_smiles": df["smiles"].nunique(),
        "mean_dock_3fap": df["dock_3fap"].mean(),
        "mean_dock_7pqv": df["dock_7pqv"].mean(),
        "mean_dock_sum": df["dock_sum"].mean(),
        "dual_hit_le_-8.5": dual_hit_rate(df, -8.5),
        "dual_hit_le_-9.0": dual_hit_rate(df, -9.0),
        "dual_hit_le_-9.5": dual_hit_rate(df, -9.5),
        "qed_mean": df["qed"].mean(),
        "sa_mean": df["sa"].mean(),
    }


def summarize_topk(df: pd.DataFrame, k: int, rank_by: str) -> dict:
    ascending = rank_by == "dock_sum"
    sub = df.sort_values(rank_by, ascending=ascending).head(k)
    row = {
        "rank_by": rank_by,
        "topk": k,
        "n": len(sub),
        "mean_dock_3fap": sub["dock_3fap"].mean(),
        "mean_dock_7pqv": sub["dock_7pqv"].mean(),
        "mean_dock_sum": sub["dock_sum"].mean(),
        "dual_hit_le_-8.5": dual_hit_rate(sub, -8.5),
        "dual_hit_le_-9.0": dual_hit_rate(sub, -9.0),
        "dual_hit_le_-9.5": dual_hit_rate(sub, -9.5),
        "qed_mean": sub["qed"].mean(),
        "sa_mean": sub["sa"].mean(),
    }
    for col in ["S_target", "S_chem", "S_priority", "dock_balance_gap"]:
        if col in sub.columns:
            row[f"mean_{col}"] = sub[col].mean()
    return row


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_contact_csv(Path(args.input_csv))
    df = add_priority_scores(df)

    pd.DataFrame([summarize_full_set(df)]).to_csv(out_dir / "contact_full_set_summary.csv", index=False)

    rows = []
    for k in args.topk:
        rows.append(summarize_topk(df, k, "dock_sum"))
        rows.append(summarize_topk(df, k, "S_priority"))
    pd.DataFrame(rows).to_csv(out_dir / "contact_topk_summary.csv", index=False)

    # Save full table with priority scores for reproducibility.
    df.to_csv(out_dir / "contact_with_priority_scores.csv", index=False)
    print(f"Saved contact-guided analysis outputs to: {out_dir}")


if __name__ == "__main__":
    main()
