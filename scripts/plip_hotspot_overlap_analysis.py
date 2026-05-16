#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Summarize candidate-level PLIP hotspot-overlap analysis.

This script expects a precomputed PLIP summary table and writes a compact
Table S8-ready file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare Table S8 from PLIP summary data.")
    parser.add_argument("--plip_summary_csv", required=True, help="Input plip_summary_top20.csv.")
    parser.add_argument("--out_csv", required=True, help="Output compact Table S8 CSV.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.plip_summary_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = pd.read_csv(input_path)
    required = [
        "candidate_label",
        "target",
        "hbond",
        "hydrophobic",
        "pistacking",
        "pication",
        "halogen",
        "all_contact_residues_n",
        "key_recovery",
        "region_recovery",
        "key_hit_residues",
        "region_hit_residues",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out = df[required].copy()
    out = out.rename(columns={
        "candidate_label": "Candidate",
        "target": "Target",
        "hbond": "H-bond",
        "hydrophobic": "Hydrophobic",
        "pistacking": "Pi-pi",
        "pication": "Pi-cation",
        "halogen": "Halogen",
        "all_contact_residues_n": "Residues_n",
        "key_recovery": "Key_recovery",
        "region_recovery": "Region_recovery",
        "key_hit_residues": "Key_residues_recovered",
        "region_hit_residues": "Reference_region_residues_recovered",
    })

    out["Candidate"] = out["Candidate"].replace({
        "PLIP_Candidate_1": "Candidate 1 (PLIP Candidate 1; Fig. 5a)",
        "PLIP_Candidate_2": "Candidate 2 (PLIP Candidate 2; Fig. 5a)",
        "PLIP_Candidate_3": "Candidate 3 (PLIP Candidate 3; Fig. 5a)",
        "PLIP_Candidate_11": "DBR-11 (PLIP Candidate 11; Fig. 5a)",
    })
    out = out.fillna("None")

    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
