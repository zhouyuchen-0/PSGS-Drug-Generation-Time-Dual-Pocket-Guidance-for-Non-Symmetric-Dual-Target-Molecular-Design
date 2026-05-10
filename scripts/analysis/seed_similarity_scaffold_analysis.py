#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Seed-fragment similarity and scaffold-overlap analysis for PSGS-Drug.

This script tests whether contact-guided Top-k molecules are trivial
seed-fragment analogues by calculating seed-fragment Tanimoto similarity,
substructure containment, and Bemis-Murcko scaffold overlap.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


DEFAULT_SEEDS = {
    "Seed_1": "C=C(F)C(F)=C(C)Nc1ccccc1",
    "Seed_2": "Fc1ccccc1",
    "Seed_3": "C=C(F)C(F)=C(C)NC(=O)C",
    "Seed_4": "CCC(=O)N(C)CNc1ccccc1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze seed-fragment similarity and scaffold overlap.")
    parser.add_argument("--input_csv", required=True, help="Contact-guided final set CSV.")
    parser.add_argument("--out_dir", required=True, help="Output directory.")
    parser.add_argument("--topk", nargs="+", type=int, default=[20, 50, 100])
    return parser.parse_args()


def mol_from_smiles(smiles: str):
    if pd.isna(smiles):
        return None
    return Chem.MolFromSmiles(str(smiles))


def canonicalize(smiles: str):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def fingerprint(mol, radius: int = 2, nbits: int = 2048):
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def scaffold_smiles(mol):
    if mol is None:
        return ""
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaf, canonical=True) if scaf and scaf.GetNumAtoms() else ""
    except Exception:
        return ""


def generic_scaffold_smiles(mol):
    if mol is None:
        return ""
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        if scaf is None or scaf.GetNumAtoms() == 0:
            return ""
        generic = MurckoScaffold.MakeScaffoldGeneric(scaf)
        return Chem.MolToSmiles(generic, canonical=True)
    except Exception:
        return ""


def prepare_seed_table() -> pd.DataFrame:
    rows = []
    for seed_id, smi in DEFAULT_SEEDS.items():
        mol = mol_from_smiles(smi)
        rows.append({
            "seed_id": seed_id,
            "seed_smiles": smi,
            "seed_canonical_smiles": canonicalize(smi),
            "seed_scaffold": scaffold_smiles(mol),
            "seed_generic_scaffold": generic_scaffold_smiles(mol),
        })
    return pd.DataFrame(rows)


def enrich_molecules(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")
    if "dock_sum" not in df.columns:
        raise ValueError("Input CSV must contain a 'dock_sum' column for Top-k analysis.")

    seed_mols = {k: mol_from_smiles(v) for k, v in DEFAULT_SEEDS.items()}
    seed_fps = {k: fingerprint(m) for k, m in seed_mols.items()}
    seed_scaffolds = {k: scaffold_smiles(m) for k, m in seed_mols.items()}
    seed_generic_scaffolds = {k: generic_scaffold_smiles(m) for k, m in seed_mols.items()}

    max_sims = []
    contains_seed = []
    same_scaffold = []
    scaffolds = []
    generic_scaffolds = []

    for smi in df["smiles"]:
        mol = mol_from_smiles(smi)
        fp = fingerprint(mol)
        sims = [DataStructs.TanimotoSimilarity(fp, sfp) for sfp in seed_fps.values() if fp is not None and sfp is not None]
        max_sims.append(max(sims) if sims else np.nan)

        contains = False
        for smol in seed_mols.values():
            if mol is not None and smol is not None and mol.HasSubstructMatch(smol):
                contains = True
                break
        contains_seed.append(contains)

        scaf = scaffold_smiles(mol)
        gscaf = generic_scaffold_smiles(mol)
        scaffolds.append(scaf)
        generic_scaffolds.append(gscaf)

        same = scaf in set(seed_scaffolds.values()) or gscaf in set(seed_generic_scaffolds.values())
        same_scaffold.append(same)

    df["max_seed_tanimoto"] = max_sims
    df["contains_any_seed"] = contains_seed
    df["BM_scaffold"] = scaffolds
    df["generic_scaffold"] = generic_scaffolds
    df["same_seed_scaffold"] = same_scaffold
    return df


def summarize_subset(df: pd.DataFrame, label: str) -> dict:
    return {
        "subset": label,
        "n": len(df),
        "mean_max_seed_tanimoto": df["max_seed_tanimoto"].mean(),
        "median_max_seed_tanimoto": df["max_seed_tanimoto"].median(),
        "maximum_seed_tanimoto": df["max_seed_tanimoto"].max(),
        "unique_bm_scaffolds": df["BM_scaffold"].nunique(),
        "unique_generic_scaffolds": df["generic_scaffold"].nunique(),
        "scaffold_molecule_ratio": df["BM_scaffold"].nunique() / len(df) if len(df) else np.nan,
        "contains_any_seed_percent": float(df["contains_any_seed"].mean() * 100.0),
        "same_seed_scaffold_percent": float(df["same_seed_scaffold"].mean() * 100.0),
        "max_seed_tanimoto_gt_0.4_percent": float((df["max_seed_tanimoto"] > 0.4).mean() * 100.0),
        "max_seed_tanimoto_gt_0.6_percent": float((df["max_seed_tanimoto"] > 0.6).mean() * 100.0),
        "max_seed_tanimoto_gt_0.8_percent": float((df["max_seed_tanimoto"] > 0.8).mean() * 100.0),
    }


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)
    df["dock_sum"] = pd.to_numeric(df["dock_sum"], errors="coerce")
    df = df.dropna(subset=["smiles", "dock_sum"]).copy()
    df = enrich_molecules(df)

    df.to_csv(out_dir / "seed_similarity_molecule_level.csv", index=False)
    prepare_seed_table().to_csv(out_dir / "seed_fragment_scaffolds.csv", index=False)

    rows = [summarize_subset(df, "Full set")]
    sorted_df = df.sort_values("dock_sum", ascending=True)
    for k in args.topk:
        rows.append(summarize_subset(sorted_df.head(k), f"Top-{k}"))

    pd.DataFrame(rows).to_csv(out_dir / "seed_similarity_scaffold_summary.csv", index=False)
    print(f"Saved seed similarity outputs to: {out_dir}")


if __name__ == "__main__":
    main()
