#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run PLIP analysis for Top-20 candidate molecules.

This script assumes candidate docking poses or complex PDB files are available.
If docking and complex construction are performed separately, use this script
only for PLIP execution and parsing.
"""

from __future__ import annotations

import argparse
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PLIP for Top-20 candidate complex PDB files.")
    parser.add_argument("--top20_csv", required=True, help="Top-20 candidate CSV with candidate_label and smiles.")
    parser.add_argument("--complex_dir", required=True, help="Directory containing candidate/target complex PDB files.")
    parser.add_argument("--out_dir", required=True, help="PLIP output directory.")
    parser.add_argument("--plip_cmd", default="plip", help="PLIP command, e.g., plip or python -m plip.plipcmd.")
    return parser.parse_args()


def run_command(command: list[str]) -> None:
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with code {result.returncode}\n"
            f"Command: {' '.join(command)}\n"
            f"STDOUT: {result.stdout}\n"
            f"STDERR: {result.stderr}"
        )


def parse_plip_xml(xml_file: Path) -> dict:
    if not xml_file.exists():
        return {}
    root = ET.parse(xml_file).getroot()
    counts = {
        "hbond": 0,
        "hydrophobic": 0,
        "pistacking": 0,
        "pication": 0,
        "saltbridge": 0,
        "waterbridge": 0,
        "halogen": 0,
        "metal": 0,
    }
    residues = set()

    for elem in root.iter():
        tag = elem.tag.lower()
        if "hydrophobic" in tag:
            counts["hydrophobic"] += 1
        elif "hbond" in tag or "hydrogen" in tag:
            counts["hbond"] += 1
        elif "pistacking" in tag:
            counts["pistacking"] += 1
        elif "pication" in tag:
            counts["pication"] += 1
        elif "saltbridge" in tag:
            counts["saltbridge"] += 1
        elif "waterbridge" in tag:
            counts["waterbridge"] += 1
        elif "halogen" in tag:
            counts["halogen"] += 1
        elif "metal" in tag:
            counts["metal"] += 1

        if elem.tag.lower().endswith("restype") or elem.tag.lower().endswith("resnr"):
            pass

    # Residue parsing is PLIP-version-dependent. If detailed residue extraction is
    # required, use the supplied plip_summary_top20.csv generated during the study.
    counts["all_contact_residues"] = ""
    counts["all_contact_residues_n"] = 0
    return counts


def main() -> None:
    args = parse_args()
    top20 = pd.read_csv(args.top20_csv)
    complex_dir = Path(args.complex_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for _, row in top20.iterrows():
        candidate = str(row["candidate_label"])
        for target in ["3FAP", "7PQV"]:
            complex_pdb = complex_dir / candidate / target / f"{candidate}_{target}_complex.pdb"
            if not complex_pdb.exists():
                continue

            target_out = out_dir / candidate / target
            target_out.mkdir(parents=True, exist_ok=True)
            command = args.plip_cmd.split() + ["-f", str(complex_pdb), "-o", str(target_out), "-x"]
            run_command(command)

            xml_files = list(target_out.glob("*.xml"))
            xml_file = xml_files[0] if xml_files else None
            parsed = parse_plip_xml(xml_file) if xml_file else {}

            parsed.update({
                "candidate_label": candidate,
                "smiles": row.get("smiles", ""),
                "target": target,
                "complex_pdb": str(complex_pdb),
                "plip_xml": str(xml_file) if xml_file else "",
            })
            rows.append(parsed)

    out_csv = out_dir / "plip_summary_top20_raw.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
