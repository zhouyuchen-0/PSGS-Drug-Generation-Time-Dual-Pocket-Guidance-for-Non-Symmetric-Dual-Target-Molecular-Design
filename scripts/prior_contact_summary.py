#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Summarize formal prior-guided and contact-guided PSGS-Drug result sets.

Outputs:
- full_quality_summary_prior_vs_contact.csv
- topk_by_docksum_prior_vs_contact.csv
- topk_by_integrated_priority_prior_vs_contact.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize prior-guided versus contact-guided result sets.")
    parser.add_argument("--prior_csv", required=True, help="Prior-guided population set CSV.")
    parser.add_argument("--contact_csv", required=True, help="Contact-guided final set CSV.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--topk", nargs="+", type=int, default=[10, 20, 50, 100, 200])
    return parser.parse_args()


def load_result(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    df = pd.read_csv(path).copy()
    required = ["smiles", "dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name} missing columns: {missing}")
    for col in ["dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa", "lip_viol", "is_novel"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["set_name"] = label
    return df


def dual_hit_rate(df: pd.DataFrame, threshold: float) -> float:
    return float(((df["dock_3fap"] <= threshold) & (df["dock_7pqv"] <= threshold)).mean() * 100.0)


def summarize(df: pd.DataFrame, label: str) -> dict:
    return {
        "set_name": label,
        "n": len(df),
        "unique_smiles": df["smiles"].nunique(),
        "mean_dock_3fap": df["dock_3fap"].mean(),
        "mean_dock_7pqv": df["dock_7pqv"].mean(),
        "mean_dock_sum": df["dock_sum"].mean(),
        "dual_hit_le_-8.0": dual_hit_rate(df, -8.0),
        "dual_hit_le_-8.5": dual_hit_rate(df, -8.5),
        "dual_hit_le_-9.0": dual_hit_rate(df, -9.0),
        "dual_hit_le_-9.5": dual_hit_rate(df, -9.5),
        "qed_mean": df["qed"].mean(),
        "sa_mean": df["sa"].mean(),
    }


def summarize_topk(df: pd.DataFrame, label: str, k: int) -> dict:
    sub = df.sort_values("dock_sum", ascending=True).head(k)
    return {
        "set_name": label,
        "topk": k,
        "mean_dock_3fap": sub["dock_3fap"].mean(),
        "mean_dock_7pqv": sub["dock_7pqv"].mean(),
        "mean_dock_sum": sub["dock_sum"].mean(),
        "strict_dual_hit_le_-9.5": dual_hit_rate(sub, -9.5),
        "qed_mean": sub["qed"].mean(),
        "sa_mean": sub["sa"].mean(),
    }


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prior = load_result(Path(args.prior_csv), "prior_guided_population_set")
    contact = load_result(Path(args.contact_csv), "contact_guided_final_set")

    pd.DataFrame([summarize(prior, "prior_guided_population_set"), summarize(contact, "contact_guided_final_set")]).to_csv(
        out_dir / "full_quality_summary_prior_vs_contact.csv", index=False
    )

    rows = []
    for label, df in [("prior_guided_population_set", prior), ("contact_guided_final_set", contact)]:
        for k in args.topk:
            rows.append(summarize_topk(df, label, k))
    pd.DataFrame(rows).to_csv(out_dir / "topk_by_docksum_prior_vs_contact.csv", index=False)

    print(f"Saved prior/contact summaries to: {out_dir}")


if __name__ == "__main__":
    main()
