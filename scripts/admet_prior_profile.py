#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute medicinal-chemistry and ADMET pass-rate summaries.

This script is intended for the prior-guided population set or any
ADMETLab-exported molecule table used in the PSGS-Drug manuscript.

Inputs
------
A CSV file containing molecule-level descriptors and/or ADMET endpoints.

Outputs
-------
1. descriptor_summary.csv
2. pass_rate_summary.csv

The script only summarizes processed values. It does not perform ADMET
prediction by itself.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute descriptor and ADMET pass-rate summaries."
    )
    parser.add_argument("--input_csv", required=True, help="Input ADMET/descriptor CSV file.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument(
        "--manual_min_freq_gt_100",
        type=float,
        default=None,
        help="Optional manually verified ring-system pass rate if not present in the CSV.",
    )
    parser.add_argument(
        "--manual_bm_scaffold_in_zinc20",
        type=float,
        default=None,
        help="Optional manually verified Bemis-Murcko scaffold-in-ZINC20 rate if not present in the CSV.",
    )
    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")


def to_numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        raise KeyError(f"Required column is missing: {col}")
    return pd.to_numeric(df[col], errors="coerce")


def find_column(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower_map = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def summarize_numeric(df: pd.DataFrame, col: str, label: str | None = None) -> dict:
    s = to_numeric(df, col).dropna()
    if s.empty:
        return {
            "metric": label or col,
            "n_nonmissing": 0,
            "mean": np.nan,
            "sd": np.nan,
            "median": np.nan,
            "min": np.nan,
            "q1": np.nan,
            "q3": np.nan,
            "max": np.nan,
        }
    return {
        "metric": label or col,
        "n_nonmissing": int(s.shape[0]),
        "mean": float(s.mean()),
        "sd": float(s.std(ddof=1)),
        "median": float(s.median()),
        "min": float(s.min()),
        "q1": float(s.quantile(0.25)),
        "q3": float(s.quantile(0.75)),
        "max": float(s.max()),
    }


def percentage(mask: pd.Series) -> float:
    return float(mask.fillna(False).mean() * 100.0)


def threshold_pass_rate(
    df: pd.DataFrame,
    col: str,
    op: str,
    value: float,
    label: str,
) -> dict:
    s = to_numeric(df, col)
    if op == ">":
        mask = s > value
    elif op == ">=":
        mask = s >= value
    elif op == "<":
        mask = s < value
    elif op == "<=":
        mask = s <= value
    else:
        raise ValueError(f"Unsupported operator: {op}")

    return {
        "metric": label,
        "pass_rate_percent": percentage(mask),
        "n_pass": int(mask.fillna(False).sum()),
        "n_total": int(len(df)),
        "source": "calculated_from_file",
    }


def range_pass_rate(df: pd.DataFrame, col: str, low: float, high: float, label: str) -> dict:
    s = to_numeric(df, col)
    mask = (s >= low) & (s <= high)
    return {
        "metric": label,
        "pass_rate_percent": percentage(mask),
        "n_pass": int(mask.fillna(False).sum()),
        "n_total": int(len(df)),
        "source": "calculated_from_file",
    }


def structural_alert_clean_rate(df: pd.DataFrame, col: str, label: str) -> dict:
    if col not in df.columns:
        return {
            "metric": label,
            "pass_rate_percent": np.nan,
            "n_pass": np.nan,
            "n_total": len(df),
            "source": "missing_column",
        }

    s = df[col].astype(str).str.strip()
    clean = (s == "['-']") | (s == "-") | (s.str.lower() == "nan") | (s == "")
    return {
        "metric": label,
        "pass_rate_percent": percentage(clean),
        "n_pass": int(clean.sum()),
        "n_total": int(len(df)),
        "source": "calculated_from_file",
    }


def manual_rate(label: str, value: float | None, n_total: int) -> dict:
    return {
        "metric": label,
        "pass_rate_percent": np.nan if value is None else float(value),
        "n_pass": np.nan,
        "n_total": int(n_total),
        "source": "manual_external" if value is not None else "not_provided",
    }


def main() -> None:
    args = parse_args()
    input_csv = Path(args.input_csv)
    out_dir = Path(args.out_dir)
    require_file(input_csv)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv)

    descriptor_aliases = {
        "QED": ["qed", "QED"],
        "Synth": ["Synth", "synth", "sa", "SA", "SA_score"],
        "MW": ["MW", "mw", "MolWt", "molecular_weight"],
        "logP": ["logP", "LogP", "mol_logp"],
        "TPSA": ["TPSA", "tpsa"],
    }

    descriptor_rows = []
    for label, aliases in descriptor_aliases.items():
        col = find_column(df, aliases)
        if col is not None:
            descriptor_rows.append(summarize_numeric(df, col, label))

    descriptor_df = pd.DataFrame(descriptor_rows)
    descriptor_df.to_csv(out_dir / "descriptor_summary.csv", index=False)

    pass_rows = []

    qed_col = find_column(df, ["qed", "QED"])
    synth_col = find_column(df, ["Synth", "synth", "sa", "SA", "SA_score"])
    logp_col = find_column(df, ["logP", "LogP", "mol_logp"])
    lip_col = find_column(df, ["lip_viol", "Lipinski_violations", "lipinski_violations"])

    if qed_col:
        pass_rows.append(threshold_pass_rate(df, qed_col, ">", 0.50, "QED > 0.50"))
        pass_rows.append(threshold_pass_rate(df, qed_col, ">", 0.67, "QED > 0.67"))
    if logp_col:
        pass_rows.append(range_pass_rate(df, logp_col, 0.0, 5.0, "0 <= logP <= 5"))
    if synth_col:
        pass_rows.append(threshold_pass_rate(df, synth_col, "<", 4.0, "Synth < 4"))
    if lip_col:
        pass_rows.append(threshold_pass_rate(df, lip_col, "<=", 0.0, "Lipinski compliant"))

    pass_rows.append(structural_alert_clean_rate(df, "PAINS", "PAINS-clean"))
    pass_rows.append(structural_alert_clean_rate(df, "BMS", "BMS-clean"))
    pass_rows.append(structural_alert_clean_rate(df, "Chelating", "Chelating-clean"))

    # Ring/scaffold metrics may be precomputed externally.
    pass_rows.append(manual_rate("min_freq > 100", args.manual_min_freq_gt_100, len(df)))
    pass_rows.append(manual_rate("BM scaffold in ZINC20", args.manual_bm_scaffold_in_zinc20, len(df)))

    admet_candidates = [
        ("Promiscuous", "Promiscuous < 0.5"),
        ("hERG", "hERG < 0.5"),
        ("Aggregators", "Aggregators < 0.5"),
        ("Carcinogenicity", "Carcinogenicity < 0.5"),
        ("Ames", "Ames < 0.5"),
        ("SkinSen", "SkinSen < 0.5"),
        ("DILI", "DILI < 0.5"),
        ("Genotoxicity", "Genotoxicity < 0.5"),
    ]
    for col, label in admet_candidates:
        if col in df.columns:
            pass_rows.append(threshold_pass_rate(df, col, "<", 0.5, label))

    pass_df = pd.DataFrame(pass_rows)
    pass_df.to_csv(out_dir / "pass_rate_summary.csv", index=False)

    print(f"Saved descriptor summary to: {out_dir / 'descriptor_summary.csv'}")
    print(f"Saved pass-rate summary to: {out_dir / 'pass_rate_summary.csv'}")


if __name__ == "__main__":
    main()
