#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Summarize higher-exhaustiveness repeated docking results.

This script does not perform docking. It summarizes a raw repeated-docking CSV
with original and higher-exhaustiveness scores.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize higher-exhaustiveness redocking results.")
    parser.add_argument("--input_csv", required=True, help="Raw repeated-docking CSV.")
    parser.add_argument("--out_summary", required=True, help="Output summary CSV.")
    parser.add_argument("--out_per_molecule", default=None, help="Optional per-molecule output CSV.")
    parser.add_argument("--strict_threshold", type=float, default=-9.5)
    return parser.parse_args()


def dual_hit_rate(df: pd.DataFrame, c1: str, c2: str, threshold: float) -> float:
    return float(((df[c1] <= threshold) & (df[c2] <= threshold)).mean() * 100.0)


def spearman_rank(x: pd.Series, y: pd.Series) -> float:
    return float(pd.Series(x).rank().corr(pd.Series(y).rank()))


def summarize(df: pd.DataFrame, name: str, threshold: float) -> dict:
    return {
        "Subset": name,
        "N": len(df),
        "Original exhaustiveness": 4,
        "Repeated exhaustiveness": 16,
        "Original mean dock_sum": df["orig_dock_sum"].mean(),
        "Repeated mean dock_sum": df["high_exh_dock_sum"].mean(),
        "Delta mean dock_sum high_minus_original": (df["high_exh_dock_sum"] - df["orig_dock_sum"]).mean(),
        "Spearman rho original_vs_repeated_dock_sum": spearman_rank(df["orig_dock_sum"], df["high_exh_dock_sum"]) if len(df) > 2 else np.nan,
        "Original strict dual-hit <= -9.5 (%)": dual_hit_rate(df, "orig_dock_3FAP", "orig_dock_7PQV", threshold),
        "Repeated strict dual-hit <= -9.5 (%)": dual_hit_rate(df, "high_exh_dock_3FAP", "high_exh_dock_7PQV", threshold),
    }


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input_csv)

    required = [
        "orig_rank_by_dock_sum",
        "selected_group",
        "orig_dock_3FAP",
        "orig_dock_7PQV",
        "orig_dock_sum",
        "high_exh_dock_3FAP",
        "high_exh_dock_7PQV",
        "high_exh_dock_sum",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for col in required:
        if col != "selected_group":
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["high_exh_dock_3FAP", "high_exh_dock_7PQV", "high_exh_dock_sum"]).copy()

    rows = [
        summarize(df[df["orig_rank_by_dock_sum"] <= 20], "Top-20 by original dock_sum", args.strict_threshold),
        summarize(df[df["orig_rank_by_dock_sum"] <= 50], "Top-50 by original dock_sum", args.strict_threshold),
    ]

    dbr = df[df["selected_group"].astype(str).eq("DBR-11")]
    if len(dbr) > 0:
        rows.append(summarize(dbr, "DBR-11", args.strict_threshold))

    out_summary = Path(args.out_summary)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_summary, index=False)

    if args.out_per_molecule:
        out_per = Path(args.out_per_molecule)
        out_per.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_per, index=False)

    print(f"Saved: {out_summary}")


if __name__ == "__main__":
    main()
