#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Start the PSGS-Drug Token-Mol-derived pocket-prior service."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
server = ROOT / "psgs_model" / "token_prior" / "tokenmol_prior_server.py"
if not server.exists():
    # fallback for flat reviewer package
    server = ROOT / "tokenmol_prior_server.py"
if not server.exists():
    raise FileNotFoundError(f"Cannot find tokenmol_prior_server.py at {server}")
subprocess.run([sys.executable, str(server)], cwd=str(server.parent), check=True)
