#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unified PSGS-Drug generation runner.

This script provides a project-root entry point for running different PSGS-Drug
generation configurations with project-relative paths.

Examples
--------
Default contact-guided final-generation demo:
    python run_generation_demo.py

Run a specific mode:
    python run_generation_demo.py --mode prior
    python run_generation_demo.py --mode contact
    python run_generation_demo.py --mode a0
    python run_generation_demo.py --mode a1
    python run_generation_demo.py --mode a2
    python run_generation_demo.py --mode b1
    python run_generation_demo.py --mode b2

Run a custom configuration file:
    python run_generation_demo.py --cfg configs/setting_prior_contactfrag.yaml
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent


MODE_TO_CONFIG = {
    "prior": "configs/setting_prior.yaml",
    "contact": "configs/setting_prior_contactfrag.yaml",
    "a0": "configs/setting_ablation_A0_noprior_noseed.yaml",
    "a1": "configs/setting_ablation_A1_prior_only.yaml",
    "a2": "configs/setting_ablation_A2_prior_plus_seed.yaml",
    "b1": "configs/setting_baseline_B1_3fap_only.yaml",
    "b2": "configs/setting_baseline_B2_7pqv_only.yaml",
}


MODE_TO_SCRIPT = {
    "prior": "psgs_model/generation/sbmolgen_prior.py",
    "contact": "psgs_model/generation/sbmolgen_contactfrag.py",
    "a0": "psgs_model/generation/sbmolgen_A0A1A2.py",
    "a1": "psgs_model/generation/sbmolgen_A0A1A2.py",
    "a2": "psgs_model/generation/sbmolgen_A0A1A2.py",
    "b1": "psgs_model/generation/sbmolgen_B1B2.py",
    "b2": "psgs_model/generation/sbmolgen_B1B2.py",
}


def resolve_path(path_like: str | Path) -> Path:
    """Resolve a project-relative path."""
    path = Path(path_like)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def infer_mode_from_cfg(cfg_path: Path) -> str:
    """Infer mode from config filename when --cfg is used."""
    name = cfg_path.name.lower()

    if "contactfrag" in name or "contact" in name:
        return "contact"
    if "prior" in name and "ablation" not in name:
        return "prior"
    if "a0" in name:
        return "a0"
    if "a1" in name:
        return "a1"
    if "a2" in name:
        return "a2"
    if "b1" in name:
        return "b1"
    if "b2" in name:
        return "b2"

    raise ValueError(
        f"Cannot infer generation mode from config filename: {cfg_path.name}. "
        "Please specify --mode explicitly."
    )


def check_required_files(mode: str, cfg_path: Path, script_path: Path) -> None:
    """Check whether key files exist before running."""
    required = [
        cfg_path,
        script_path,
        PROJECT_ROOT / "psgs_model" / "generation" / "RNN-model" / "model.json",
        PROJECT_ROOT / "psgs_model" / "generation" / "RNN-model" / "model.h5",
        PROJECT_ROOT / "data" / "250k_rndm_zinc_drugs_clean.smi",
    ]

    # Modes using pocket priors need the prior service inputs.
    if mode in {"prior", "contact", "a1", "a2", "b1", "b2"}:
        required.extend(
            [
                PROJECT_ROOT / "receptors" / "3fap.pkl",
                PROJECT_ROOT / "receptors" / "7pqv.pkl",
            ]
        )

    missing = [str(p.relative_to(PROJECT_ROOT)) for p in required if not p.exists()]

    if missing:
        raise FileNotFoundError(
            "Required files are missing:\n"
            + "\n".join(f"  - {m}" for m in missing)
            + "\n\nPlease keep the repository structure unchanged."
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run PSGS-Drug generation with project-relative paths."
    )
    parser.add_argument(
        "--mode",
        choices=sorted(MODE_TO_CONFIG.keys()),
        default="contact",
        help="Generation mode. Default: contact.",
    )
    parser.add_argument(
        "--cfg",
        type=str,
        default=None,
        help="Optional custom YAML configuration path.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print resolved paths without running generation.",
    )

    args = parser.parse_args()

    if args.cfg:
        cfg_path = resolve_path(args.cfg)
        mode = args.mode if args.mode else infer_mode_from_cfg(cfg_path)
        if args.mode == "contact":
            # If user supplies --cfg but does not explicitly change --mode,
            # infer mode from file name to avoid using the wrong script.
            mode = infer_mode_from_cfg(cfg_path)
    else:
        mode = args.mode
        cfg_path = resolve_path(MODE_TO_CONFIG[mode])

    script_path = resolve_path(MODE_TO_SCRIPT[mode])

    check_required_files(mode=mode, cfg_path=cfg_path, script_path=script_path)

    env = os.environ.copy()
    env["SFGDRUG_CFG_PATH"] = str(cfg_path)
    env["PSGS_PROJECT_ROOT"] = str(PROJECT_ROOT)

    print("========== PSGS-Drug generation runner ==========")
    print(f"Project root : {PROJECT_ROOT}")
    print(f"Mode         : {mode}")
    print(f"Config       : {cfg_path}")
    print(f"Script       : {script_path}")
    print("=================================================")

    if args.dry_run:
        print("[DRY RUN] Path check completed. No generation was launched.")
        return

    cmd = [sys.executable, str(script_path), "--cfg", str(cfg_path)]

    subprocess.run(
        cmd,
        cwd=str(PROJECT_ROOT),
        env=env,
        check=True,
    )


if __name__ == "__main__":
    main()