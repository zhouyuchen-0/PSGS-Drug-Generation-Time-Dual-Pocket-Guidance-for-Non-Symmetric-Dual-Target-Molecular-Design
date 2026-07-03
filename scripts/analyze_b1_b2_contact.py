#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analyze B1/B2 single-pocket conditioning-source asymmetry controls and the
contact-guided final-generation set under the shared downstream dual-target
docking protocol.

B1 denotes 3FAP-only conditioning during generation, whereas B2 denotes
7PQV-only conditioning during generation. Both libraries are subsequently
evaluated by the same retrospective two-context docking workflow. This script
does not test a sequential first-target/second-target generation order.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze B1/B2 single-pocket conditioning-source asymmetry controls "
            "and the contact-guided final-generation set."
        )
    )
    parser.add_argument(
        "--b1_csv",
        required=True,
        help=(
            "B1 CSV: 3FAP-only conditioned library evaluated by the shared "
            "retrospective two-context docking workflow."
        ),
    )
    parser.add_argument(
        "--b2_csv",
        required=True,
        help=(
            "B2 CSV: 7PQV-only conditioned library evaluated by the shared "
            "retrospective two-context docking workflow."
        ),
    )
    parser.add_argument(
        "--contact_csv",
        required=True,
        help="Contact-guided final-generation set CSV.",
    )
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--topk", nargs="+", type=int, default=[50])
    return parser.parse_args()


def load_and_standardize(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_csv(path).copy()
    if "dock_sum" not in df.columns:
        if "dock_combined" in df.columns:
            df["dock_sum"] = df["dock_combined"]
        else:
            df["dock_sum"] = pd.to_numeric(df["dock_3fap"], errors="coerce") + pd.to_numeric(df["dock_7pqv"], errors="coerce")

    if "qed" not in df.columns:
        df["qed"] = df["QED"] if "QED" in df.columns else np.nan
    if "sa" not in df.columns:
        if "SA_score" in df.columns:
            df["sa"] = df["SA_score"]
        elif "Synth" in df.columns:
            df["sa"] = df["Synth"]
        else:
            df["sa"] = np.nan
    if "lip_viol" not in df.columns:
        df["lip_viol"] = np.nan

    for col in ["dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa", "lip_viol"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["method"] = label
    return df


def dual_hit_rate(df: pd.DataFrame, threshold: float) -> float:
    return float(((df["dock_3fap"] <= threshold) & (df["dock_7pqv"] <= threshold)).mean() * 100.0)


def summarize_method(df: pd.DataFrame, method_name: str) -> dict:
    out = {
        "method": method_name,
        "n": len(df),
        "mean_dock_3fap": df["dock_3fap"].mean(),
        "mean_dock_7pqv": df["dock_7pqv"].mean(),
        "mean_dock_sum": df["dock_sum"].mean(),
        "dock_balance_gap_abs": (df["dock_3fap"] - df["dock_7pqv"]).abs().mean(),
        "dual_hit_le_-8.0": dual_hit_rate(df, -8.0),
        "dual_hit_le_-8.5": dual_hit_rate(df, -8.5),
        "dual_hit_le_-9.0": dual_hit_rate(df, -9.0),
        "dual_hit_le_-9.5": dual_hit_rate(df, -9.5),
        "qed_mean": df["qed"].mean(),
        "sa_mean": df["sa"].mean(),
    }
    if df["lip_viol"].notna().any():
        out["lipinski_pass_rate"] = float((df["lip_viol"] <= 0).mean() * 100.0)
    return out


def topk_summary(df: pd.DataFrame, method: str, k: int) -> dict:
    sub = df.sort_values("dock_sum", ascending=True).head(k)
    return {
        "method": method,
        "topk": k,
        "topk_mean_dock_sum": sub["dock_sum"].mean(),
        "topk_strict_dual_hit_le_-9.5": dual_hit_rate(sub, -9.5),
        "topk_qed_mean": sub["qed"].mean(),
        "topk_sa_mean": sub["sa"].mean(),
    }


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    datasets = {
        "B1 (3FAP-only)": load_and_standardize(Path(args.b1_csv), "B1 (3FAP-only)"),
        "B2 (7PQV-only)": load_and_standardize(Path(args.b2_csv), "B2 (7PQV-only)"),
        "Contact-guided final-generation set": load_and_standardize(
            Path(args.contact_csv),
            "Contact-guided final-generation set",
        ),
    }

    full_rows = [summarize_method(df, name) for name, df in datasets.items()]
    full_df = pd.DataFrame(full_rows)
    full_df.to_csv(out_dir / "B1_B2_conditioning_source_asymmetry_full_summary.csv", index=False)

    topk_rows = []
    for name, df in datasets.items():
        for k in args.topk:
            topk_rows.append(topk_summary(df, name, k))
    pd.DataFrame(topk_rows).to_csv(out_dir / "B1_B2_conditioning_source_asymmetry_topk_summary.csv", index=False)

    print(f"Saved B1/B2 conditioning-source-asymmetry summaries to: {out_dir}")


if __name__ == "__main__":
    main()
