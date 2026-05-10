#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Select balanced multi-criteria representative molecules from the contact-guided final set.

The selected molecules are intended for representative candidate cards and
follow-up PLIP interaction analysis. Ranking is not dock_sum-only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Select balanced multi-criteria representative molecules.")
    parser.add_argument("--input_csv", required=True, help="Contact-guided final set CSV.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--top_n", type=int, default=6)
    parser.add_argument("--dock_threshold", type=float, default=-9.0)
    parser.add_argument("--qed_min", type=float, default=0.65)
    parser.add_argument("--sa_max", type=float, default=4.5)
    parser.add_argument("--max_per_scaffold", type=int, default=1)
    parser.add_argument("--max_tanimoto", type=float, default=0.70)
    return parser.parse_args()


def minmax_scale(values, higher_is_better=True):
    x = np.asarray(values, dtype=float)
    if len(x) == 0:
        return x
    xmin, xmax = np.nanmin(x), np.nanmax(x)
    if np.isclose(xmin, xmax):
        return np.ones_like(x) * 0.5
    s = (x - xmin) / (xmax - xmin)
    return s if higher_is_better else 1.0 - s


def mol_from_smiles(smiles):
    return Chem.MolFromSmiles(str(smiles)) if pd.notna(smiles) else None


def get_scaffold(smiles):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return ""


def get_fp(smiles):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    required = ["smiles", "dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for col in ["dock_3fap", "dock_7pqv", "dock_sum", "qed", "sa", "lip_viol"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "lip_viol" not in df.columns:
        df["lip_viol"] = 0
    if "pains_clean_bool" not in df.columns:
        df["pains_clean_bool"] = True

    df = df.dropna(subset=required).copy()
    df = df[
        (df["dock_3fap"] <= args.dock_threshold)
        & (df["dock_7pqv"] <= args.dock_threshold)
        & (df["qed"] >= args.qed_min)
        & (df["sa"] <= args.sa_max)
        & (df["lip_viol"] <= 0)
        & (df["pains_clean_bool"].astype(bool))
    ].copy()

    df["dock_balance_gap"] = (df["dock_3fap"] - df["dock_7pqv"]).abs()
    df["S_dual_docking"] = minmax_scale(df["dock_sum"], higher_is_better=False)
    df["S_qed"] = minmax_scale(df["qed"], higher_is_better=True)
    df["S_sa"] = minmax_scale(df["sa"], higher_is_better=False)
    df["S_lipinski"] = 1.0
    df["S_priority"] = (
        0.50 * df["S_dual_docking"]
        + 0.30 * df["S_qed"]
        + 0.15 * df["S_sa"]
        + 0.05 * df["S_lipinski"]
    )

    df["scaffold"] = df["smiles"].apply(get_scaffold)
    df["fp"] = df["smiles"].apply(get_fp)
    df = df.sort_values("S_priority", ascending=False)

    selected = []
    scaffold_counts = {}
    selected_fps = []

    for _, row in df.iterrows():
        scaf = row["scaffold"]
        if scaffold_counts.get(scaf, 0) >= args.max_per_scaffold:
            continue

        fp = row["fp"]
        if fp is not None and selected_fps:
            max_sim = max(DataStructs.TanimotoSimilarity(fp, sfp) for sfp in selected_fps if sfp is not None)
            if max_sim > args.max_tanimoto:
                continue

        selected.append(row.drop(labels=["fp"]).to_dict())
        scaffold_counts[scaf] = scaffold_counts.get(scaf, 0) + 1
        selected_fps.append(fp)

        if len(selected) >= args.top_n:
            break

    out = pd.DataFrame(selected)
    out.insert(0, "candidate_label", [f"Candidate {i+1}" for i in range(len(out))])

    out_csv = out_dir / "top6_balanced_multi_criteria_candidates.csv"
    out_smi = out_dir / "top6_balanced_multi_criteria_candidates.smi"
    out.to_csv(out_csv, index=False)
    out[["smiles", "candidate_label"]].to_csv(out_smi, sep="\t", header=False, index=False)

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_smi}")


if __name__ == "__main__":
    main()
