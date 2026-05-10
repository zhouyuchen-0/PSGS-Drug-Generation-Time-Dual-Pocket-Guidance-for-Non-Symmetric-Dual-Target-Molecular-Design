#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analyze A1 versus A2 Top-k behavior.

A1 is the dual-pocket-prior-only ablation setting.
A2 is the dual-pocket-prior plus contact-guided prefix ablation setting.

This script summarizes Top-k molecules under a specified ranking criterion.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze A1 versus A2 Top-k metrics.")
    parser.add_argument("--a1_csv", required=True, help="A1 prior-only CSV file.")
    parser.add_argument("--a2_csv", required=True, help="A2 prior-plus-seed CSV file.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument(
        "--sort_by",
        default="dock_sum",
        choices=["dock_sum", "reward", "qed", "sa", "r_dock"],
        help="Ranking column. dock_sum and sa are sorted ascending; all others descending.",
    )
    parser.add_argument("--topk", nargs="+", type=int, default=[10, 20, 50, 100])
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    df = pd.read_csv(path)
    for col in [
        "dock_3fap",
        "dock_7pqv",
        "dock_sum",
        "reward",
        "r_dock",
        "qed",
        "sa",
        "lip_viol",
        "is_novel",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def dual_hit_rate(df: pd.DataFrame, threshold: float) -> float:
    return float(((df["dock_3fap"] <= threshold) & (df["dock_7pqv"] <= threshold)).mean() * 100.0)


def summarize_subset(df: pd.DataFrame, group: str, topk: int, sort_by: str) -> dict:
    ascending = sort_by in {"dock_sum", "sa", "lip_viol"}
    sub = df.sort_values(sort_by, ascending=ascending).head(topk).copy()

    row = {
        "group": group,
        "topk": topk,
        "sort_by": sort_by,
        "n": len(sub),
        "mean_dock_3fap": sub["dock_3fap"].mean(),
        "mean_dock_7pqv": sub["dock_7pqv"].mean(),
        "mean_dock_sum": sub["dock_sum"].mean(),
        "mean_qed": sub["qed"].mean(),
        "mean_sa": sub["sa"].mean(),
        "mean_balance_gap": (sub["dock_3fap"] - sub["dock_7pqv"]).abs().mean(),
        "dual_hit_le_-8.5": dual_hit_rate(sub, -8.5),
        "dual_hit_le_-9.0": dual_hit_rate(sub, -9.0),
        "dual_hit_le_-9.5": dual_hit_rate(sub, -9.5),
    }
    if "lip_viol" in sub.columns:
        row["lipinski_pass_rate"] = float((sub["lip_viol"] <= 0).mean() * 100.0)
    if "is_novel" in sub.columns:
        row["novelty_rate"] = float(pd.to_numeric(sub["is_novel"], errors="coerce").fillna(1).mean() * 100.0)
    return row


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_a1 = read_table(Path(args.a1_csv))
    df_a2 = read_table(Path(args.a2_csv))

    rows = []
    for k in args.topk:
        rows.append(summarize_subset(df_a1, "A1_prior_only", k, args.sort_by))
        rows.append(summarize_subset(df_a2, "A2_prior_plus_seed", k, args.sort_by))

    out = pd.DataFrame(rows)
    out_path = out_dir / "A1_A2_topk_summary.csv"
    out.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
