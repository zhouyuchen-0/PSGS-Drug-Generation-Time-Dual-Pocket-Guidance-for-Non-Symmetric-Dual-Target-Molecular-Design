#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Formal contact-guided PSGS-Drug molecular generation script.

This script starts from contact-derived seed prefixes and extends them using MCTS-GRU generation under dual-pocket prior guidance.
Configuration is supplied through --cfg or the SFGDRUG_CFG_PATH environment variable.
"""

import json
import urllib.request
import urllib.error

"""
sbmolgen_seed_dual_prior.py
------------------------------------------------------------
       (seeded generation)   :
-        /     seed   (         )
-         seed   (MCTS       )
-   Token-Mol   pocket       token   
-      docking + QED/SA +      reward
-     SFG-Drug RNN rollout(    CLM)

  :
- prior_client.py should be placed in the PSGS-Drug generation module directory.
- utils/rdock_test_MP.py(        )
-   SFG-Drug utils/*
"""

import os
import re
import sys
import time
import yaml
import csv
import math
import traceback
from typing import Any, Dict, List, Tuple, Optional

import numpy as np

from rdkit import Chem

# ===================== [fixcfg8] fallback + diag =====================
def _fallback_balance_smiles(smiles: str,
                             max_append_close: int = 6,
                             drop_unmatched_close: bool = True,
                             drop_unpaired_ring_digit: bool = True) -> str:
    """
    Very light SMILES repair (NOT full sanitization):
      - drop unmatched ')'
      - append missing ')' (capped)
      - drop unpaired ring digits (0-9)
    Goal: reduce trivial RDKit parse failures so docking becomes reachable.
    """
    if not smiles:
        return smiles

    out = []
    open_cnt = 0
    digit_cnt = {str(i): 0 for i in range(10)}

    i = 0
    while i < len(smiles):
        ch = smiles[i]

        # keep bracket atoms intact: [NH+], [C@H], etc.
        if ch == '[':
            j = smiles.find(']', i+1)
            if j == -1:
                break
            out.append(smiles[i:j+1])
            i = j + 1
            continue

        if ch == '(':
            open_cnt += 1
            out.append(ch); i += 1; continue

        if ch == ')':
            if open_cnt > 0:
                open_cnt -= 1
                out.append(ch)
            else:
                if not drop_unmatched_close:
                    out.append(ch)
            i += 1
            continue

        if ch.isdigit():
            digit_cnt[ch] += 1
            out.append(ch); i += 1; continue

        out.append(ch); i += 1

    if open_cnt > 0:
        out.append(')' * min(open_cnt, max_append_close))

    repaired = "".join(out)

    if drop_unpaired_ring_digit:
        for d, c in digit_cnt.items():
            if c % 2 == 1:
                repaired = repaired.replace(d, "")

    return repaired
# ================================================================

from rdkit.Chem import QED, rdMolDescriptors

from utils.add_node_type_zinc import (
    predict_smile, make_input_smile, expanded_node,
    node_to_add, chem_kn_simulation, check_node_type
)
from utils.load_model import loaded_model
from utils.make_smile import zinc_processed_with_bracket, zinc_data_with_bracket_original

#            docking   (  )
from utils.rdock_test_MP import vinadock_score

# Token-Mol prior client
# Token-Mol prior client
try:
    pass  #   ,      
except Exception as e:
    print(f"[ERROR]    prior_client   : {e}")
    sys.exit(1)

if DEBUG_DIAG:
    print("Debug mode enabled")
# =========================================================
#     debug   (     + YAML debug )
# =========================================================
def env_flag(name: str, default: bool = False) -> bool:
    v = os.environ.get(name, "")
    if v == "":
        return default
    return str(v).strip().lower() in {"1", "true", "yes", "y", "on"}


DEBUG_DIAG = env_flag("SFG_DEBUG_DIAG", False)

# ============================

# ===== Prior client (Token-Mol) =====
def _post_json(url: str, payload: dict, timeout_s: float = 5.0) -> dict:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            raw = resp.read().decode("utf-8", errors="ignore")
        return json.loads(raw) if raw else {}
    except urllib.error.HTTPError as e:
        try:
            raw = e.read().decode("utf-8", errors="ignore")
        except Exception:
            raw = ""
        return {"status": "error", "message": f"HTTPError {e.code}", "detail": raw}
    except Exception as e:
        return {"status": "error", "message": str(e)}

def _combine_dual(pr1, pr2, mode: str = "min"):
    if not pr1 or not pr2 or len(pr1) != len(pr2):
        return None
    if mode == "mean":
        comb = [(a + b) * 0.5 for a, b in zip(pr1, pr2)]
    elif mode == "geom":
        comb = [(max(a, 1e-12) * max(b, 1e-12)) ** 0.5 for a, b in zip(pr1, pr2)]
    else:
        comb = [a if a < b else b for a, b in zip(pr1, pr2)]
    s = float(sum(comb)) if comb else 0.0
    if s <= 0:
        return None
    return [float(x) / s for x in comb]

def call_prior(prefix_tokens, candidates, prior_url, pocket1_path, pocket2_path):
    """Token-Mol prior server schema:
    POST /prior {"prefix_tokens":[...], "candidates":[...], "protein_path":"..."}
    Dual pocket = two calls + elementwise min + renormalize.
    """
    bad_prefix = {"<|beginoftext|>", "<|mask:0|>", "<bos>", "<pad>", "<eos>", "<unk>"}
    prefix = [t for t in prefix_tokens if t and (t not in bad_prefix)]
    payload_base = {"prefix_tokens": prefix[-64:], "candidates": candidates}

    r1 = _post_json(prior_url, {**payload_base, "protein_path": pocket1_path})
    r2 = _post_json(prior_url, {**payload_base, "protein_path": pocket2_path})

    if r1.get("status") != "ok" or r2.get("status") != "ok":
        return [1.0 / max(1, len(candidates))] * len(candidates)

    pr1 = r1.get("prior", [])
    pr2 = r2.get("prior", [])
    comb = _combine_dual(pr1, pr2, mode="min")
    if comb is None:
        return [1.0 / max(1, len(candidates))] * len(candidates)
    return comb

# Seed normalization (input-only fix)
# - avoid opening a branch on a terminal halogen like "F(" / "Cl(" / "Br(" / "I(" which causes RDKit valence errors
# - keep everything else unchanged
# ============================
_HALOGEN_OPEN_PAT = re.compile(r"(Cl|Br|F|I)\($")

def normalize_seed_prefix(seed_smiles: str) -> str:
    """  seed    "      "     Token-Mol/VAL      [*].

        (  "   "   ):
    -       seedmap         "("(    "   "     SMILES     )
    -          "("(  "...F(" / "...Cl(" / "...Br(" / "...I("),         [*]
      (    =1,           valence error)
    -       "   '(' -> '[*]'"     
    """
    s = (seed_smiles or "").strip()

    # A)      "(":    ,           
    #     :...cc1F(  -> ...cc1[*]
    s = re.sub(r"(Cl|Br|I|F)\($", "(", s)

    # B)    "(":    [*](    )
    if s.endswith("("):
        s = s[:-1] + "[*]"

    # C)   :   seed      [*]      ;       '*'     [*]
    if s.endswith("*") and not s.endswith("[*]"):
        s = s[:-1] + "[*]"

    return s

def sanitize_smiles_for_eval(smi: str) -> str:
    """    SMILES      ,   RDKit / docking        [*]        .
      :    "    "  ,        .
    """
    if smi is None:
        return ""
    s = smi.strip()

    # 1)    Token-Mol/VAL   dummy    (          )
    s = s.replace("[*]", "")

    # 2)      vocab        '*',   
    s = s.replace("*", "")

    # 3)           
    s = s.replace("()", "")

    return s
DEBUG_TRACEBACK = env_flag("SFG_DEBUG_TRACEBACK", False)


def dprint(*args, **kwargs):
    if DEBUG_DIAG:
        print(*args, **kwargs)


# =========================================================
#  RDKit   
# =========================================================

def _quick_smiles_guard(s: str) -> bool:
    """Fast syntax guard for SMILES: only enforces 3 minimal constraints.
    1) Parentheses must be balanced and never go negative while scanning.
    2) Ring digits 1-9 must be paired (even count); disallow consecutive digits.
    3) No other chemistry/grammar rules are imposed (avoid behavior changes).
    """
    if not s:
        return False
    open_par = 0
    ring_parity = {d: 0 for d in "123456789"}  # 1 means currently unpaired
    prev_is_digit = False
    for ch in s:
        if ch == "(":
            open_par += 1
            prev_is_digit = False
            continue
        if ch == ")":
            open_par -= 1
            if open_par < 0:
                return False
            prev_is_digit = False
            continue
        if ch.isdigit():
            if prev_is_digit:
                return False
            if ch in ring_parity:
                ring_parity[ch] ^= 1
            prev_is_digit = True
        else:
            prev_is_digit = False
    if open_par != 0:
        return False
    if any(v == 1 for v in ring_parity.values()):
        return False
    return True

def safe_mol_from_smiles(smi):
    """
    RDKit parse wrapper with lightweight fallback.
    Enable diag via: export DIAG_SMILES=1
    """
    import os
    diag = (os.environ.get("DIAG_SMILES", "0") == "1")

    try:
        m = Chem.MolFromSmiles(smi)
    except Exception:
        m = None

    if m is not None:
        if diag:
            try: print(f"[diag] mol_ok=1 len={len(smi)} head={smi[:80]}", flush=True)
            except Exception: pass
        return m

    # fallback repair (parentheses / ring digits)
    smi_fb = _fallback_balance_smiles(smi)
    if smi_fb != smi:
        try:
            m2 = Chem.MolFromSmiles(smi_fb)
        except Exception:
            m2 = None
        if m2 is not None:
            if diag:
                try: print(f"[diag] mol_ok=1 fallback=1 len={len(smi_fb)} head={smi_fb[:80]}", flush=True)
                except Exception: pass
            return m2

    if diag:
        try: print(f"[diag] mol_ok=0 len={len(smi)} head={smi[:80]}", flush=True)
        except Exception: pass
    return None
def safe_qed(smi: str) -> float:
    try:
        m = safe_mol_from_smiles(smi)
        return float(QED.qed(m)) if m else 0.0
    except Exception:
        return 0.0


def safe_sa(smi: str) -> float:
    """
       SA   (        ,       )
         1~10,    
    """
    try:
        m = safe_mol_from_smiles(smi)
        if not m:
            return 6.0
        heavy = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)
        rings = int(Chem.GetSSSR(m))
        sa = 2.0 + 0.20 * max(0, heavy - 20) + 0.25 * max(0, rings - 2)
        return float(min(max(sa, 1.0), 10.0))
    except Exception:
        return 6.0


def lipinski_violations(smi: str) -> int:
    try:
        m = safe_mol_from_smiles(smi)
        if m is None:
            return 99
        mw = rdMolDescriptors.CalcExactMolWt(m)
        logp = rdMolDescriptors.CalcCrippenDescriptors(m)[0]
        hbd = rdMolDescriptors.CalcNumHBD(m)
        hba = rdMolDescriptors.CalcNumHBA(m)
        viol = 0
        if mw > 500:
            viol += 1
        if logp > 5:
            viol += 1
        if hbd > 5:
            viol += 1
        if hba > 10:
            viol += 1
        return int(viol)
    except Exception:
        return 99


def morgan_fp_bytes(smi: str, nbits: int = 2048) -> bytes:
    try:
        m = safe_mol_from_smiles(smi)
        if m is None:
            return smi.encode("utf-8", errors="ignore")
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius=2, nBits=nbits)
        return fp.ToBitString().encode("ascii")
    except Exception:
        return smi.encode("utf-8", errors="ignore")


# =========================================================
#      (     mcts.hours / hours_per_seed)
# =========================================================
def load_conf() -> Dict[str, Any]:
    """
       :
      1)       --cfg
      2)      SFGDRUG_CFG_PATH
      3)       
    """
    import sys
    dprint("load_conf started")

    # 1)       --cfg
    for i, arg in enumerate(sys.argv):
        if arg == "--cfg" and i + 1 < len(sys.argv):
            cfg_path = sys.argv[i + 1]
            dprint(f"config path from --cfg: {cfg_path}")
            if os.path.exists(cfg_path):
                with open(cfg_path, "r", encoding="utf-8") as f:
                    data = yaml.load(f, Loader=yaml.SafeLoader) or {}
                data["_cfg_path"] = cfg_path
                dprint("config loaded from command line")
                return data
            else:
                print(f"DEBUG: --cfg file not found: {cfg_path}")

    # 2)     
    env_cfg = os.environ.get("SFGDRUG_CFG_PATH", "").strip()
    if env_cfg:
        dprint(f"config path from SFGDRUG_CFG_PATH: {env_cfg}")
        if os.path.exists(env_cfg):
            with open(env_cfg, "r", encoding="utf-8") as f:
                data = yaml.load(f, Loader=yaml.SafeLoader) or {}
            data["_cfg_path"] = env_cfg
            dprint("config loaded from environment variable")
            return data

    # 3)     
    candidates = [
        "configs/setting_prior_contactfrag.yaml",
        os.path.join(os.getcwd(), "setting_prior_contactfrag.yaml"),
        "setting_prior_contactfrag.yaml",
    ]
    for p in candidates:
        if p and os.path.exists(p):
            print(f"DEBUG: found config at {p}")
            with open(p, "r", encoding="utf-8") as f:
                data = yaml.load(f, Loader=yaml.SafeLoader) or {}
            data["_cfg_path"] = p
            print("DEBUG: config loaded from default")
            return data

    print("DEBUG: no config file found, about to raise error")
    raise FileNotFoundError("    seed     .        SFGDRUG_CFG_PATH     --cfg   .")


CONF = load_conf()
MCONF = CONF.get("mcts", {}) or {}
RCONF = CONF.get("reward", {}) or {}
PCONF = CONF.get("prefilter", {}) or {}
DCONF = CONF.get("debug", {}) or {}

#    debug   :      
DEBUG_DIAG = env_flag("SFG_DEBUG_DIAG", bool(DCONF.get("diag", False)))
DEBUG_TRACEBACK = env_flag("SFG_DEBUG_TRACEBACK", bool(DCONF.get("traceback", False)))

TARGET = str(CONF.get("target", "3fap-7pqv-prior-seed"))
TRIAL = int(CONF.get("trial", 1))
TARGET_TOTAL_MOLS = int(CONF.get("target_total_mols", 200))
LOOP_NODE_EXP = int(CONF.get("loop_num_nodeExpansion", 82))
SIM_NUM = int(CONF.get("simulation_num", 3))

# hours   :   hours_per_seed,   mcts.hours,     hours
HOURS_PER_SEED = float(
    CONF.get(
        "hours_per_seed",
        MCONF.get("hours", CONF.get("hours", 1.0))
    )
)

SA_THRESHOLD = float(CONF.get("sa_threshold", 5.0))
RULE5 = int(CONF.get("rule5", 1))
RADICAL_CHECK = bool(CONF.get("radical_check", True))
HASHIMOTO_FILTER = bool(CONF.get("hashimoto_filter", True))
MODEL_NAME = str(CONF.get("model_name", "model"))

BASE_VINA_3FAP = float(CONF.get("base_vina_3fap", CONF.get("base_vinadock_score", -7.0)))
BASE_VINA_7PQV = float(CONF.get("base_vina_7pqv", CONF.get("base_vinadock_score", -7.0)))

# MCTS
C_VAL = float(MCONF.get("c_val", 1.3))
EXPANDED_K = int(MCONF.get("expanded_k", 24))
DIVERSITY_BONUS = float(MCONF.get("diversity_bonus", 0.05))

# reward weights
W_DOCK = float(RCONF.get("w_docking", 0.75))
W_QED = float(RCONF.get("w_qed", 0.20))
W_SA = float(RCONF.get("w_sa", 0.05))
W_LIP = float(RCONF.get("w_lip_penalty", 0.03))
W_PAINS = float(RCONF.get("w_pains_penalty", 0.00))  #         PAINS,    

# prefilter
PF_QED_MIN = float(PCONF.get("qed_min", 0.0))
PF_SA_MAX = float(PCONF.get("sa_max", 10.0))
PF_LIP_MAX_VIOL = int(PCONF.get("lipinski_max_viol", 99))
DOCK_MAX_PER_ROLLOUT = int(PCONF.get("dock_max_per_rollout", 999999))

# prior
PRIOR_URL = CONF.get("prior_url", None)
P1 = CONF.get("prior_protein_path_1", None)
P2 = CONF.get("prior_protein_path_2", None)

# seed
SEED_ENABLE = bool(CONF.get("seed_enable", True))
SEED_MAX_NUM = int(CONF.get("seed_max_num", 8))
SEED_SCAFFOLDS_OVERRIDE = CONF.get("seed_scaffolds_override", []) or []
SEED_SMILES_OVERRIDE = CONF.get("seed_smiles_override", []) or []  #         (    SMILES)
SEED_DRUGS = CONF.get("seed_drugs", []) or []

#   
#   
OUTPUT_CSV = str(CONF.get("output", f"result_{TARGET}_trial{TRIAL}.csv"))
RAW_OUTPUT_CSV = OUTPUT_CSV.replace(".csv", ".raw.csv")
SEED_STATS_CSV = OUTPUT_CSV.replace(".csv", ".seed_stats.csv")
DEBUG_LOG = str(CONF.get("debug_log", f"debug_{TARGET}_trial{TRIAL}.log"))
SEEDMAP_CSV = str(CONF.get("seedmap", f"seedmap_{TARGET}_trial{TRIAL}.csv"))

# smoke bypass(      )
SMOKE_BYPASS_RADICAL = env_flag("SFG_SMOKE_BYPASS_RADICAL", bool(DCONF.get("smoke_bypass_radical", False)))
SMOKE_BYPASS_HASHIMOTO = env_flag("SFG_SMOKE_BYPASS_HASHIMOTO", bool(DCONF.get("smoke_bypass_hashimoto", False)))

print("====== SEEDED PRIOR CONFIG ======")
print("cfg_path:", CONF.get("_cfg_path", "unknown"))
print("target:", TARGET, "trial:", TRIAL)
print("target_total_mols:", TARGET_TOTAL_MOLS)
print("hours_per_seed:", HOURS_PER_SEED)
print("mcts: c_val=", C_VAL, "expanded_k=", EXPANDED_K, "sim_num=", SIM_NUM)
print("prefilter:", {
    "qed_min": PF_QED_MIN, "sa_max": PF_SA_MAX,
    "lip_max_viol": PF_LIP_MAX_VIOL, "dock_max_per_rollout": DOCK_MAX_PER_ROLLOUT
})
print("reward:", {"w_dock": W_DOCK, "w_qed": W_QED, "w_sa": W_SA, "w_lip": W_LIP})
print("prior_url:", PRIOR_URL)
print("pocket1:", P1)
print("pocket2:", P2)
print("seed_enable:", SEED_ENABLE, "seed_max_num:", SEED_MAX_NUM)
print("seed_smiles_override:", len(SEED_SMILES_OVERRIDE))
print("seed_scaffolds_override:", len(SEED_SCAFFOLDS_OVERRIDE))
print("seed_drugs:", len(SEED_DRUGS))
print("radical_check:", RADICAL_CHECK, "hashimoto_filter:", HASHIMOTO_FILTER)
print("SMOKE_BYPASS_RADICAL:", SMOKE_BYPASS_RADICAL)
print("SMOKE_BYPASS_HASHIMOTO:", SMOKE_BYPASS_HASHIMOTO)
print("DEBUG_DIAG:", DEBUG_DIAG, "DEBUG_TRACEBACK:", DEBUG_TRACEBACK)
print("===============================")

#   :       fragment prefix     
SEED_FRAGMENT_PREFIXES = CONF.get("seed_prefixes_override", []) or []
#   :    extract_contact_seed_fragments_v2.py      (            )
# - seed_scaffolds_no_star_override:       prefix(   ,   '*'         )
# - seed_scaffolds_override:   '*'    (     vocab    '*',   ;          skip)
if not SEED_FRAGMENT_PREFIXES:
    SEED_FRAGMENT_PREFIXES = CONF.get("seed_scaffolds_no_star_override", []) or []
if not SEED_FRAGMENT_PREFIXES:
    SEED_FRAGMENT_PREFIXES = CONF.get("seed_scaffolds_override", []) or []

SEED_FRAGMENT_MODE = CONF.get("seed_fragment_mode", "prefer")        #    prefer
SEED_FRAGMENT_MAX_LEN = int(CONF.get("seed_fragment_max_len", 120))
SEED_FRAGMENT_MAX_LEN_TOKENS = int(CONF.get("seed_fragment_max_len_tokens", 80))
#       :     sp3   (     "     "   /      )
SEED_FRAGMENT_PREFER_SP3 = bool(CONF.get("seed_fragment_prefer_sp3", True))
#   :  smarts     (    /  ,      ;  "    "   )
SEED_SCAFFOLDS_SMARTS_OVERRIDE = CONF.get("seed_scaffolds_smarts_override", []) or []
DEFAULT_FRAGMENT_PREFIXES = []  #       (   )

# =========================================================
#       / MCTS   
# =========================================================
class ChemicalState:
    __slots__ = ('position',)
    def __init__(self, init_tokens: Optional[List[str]] = None):
        self.position = ["&"]
        if init_tokens:
            self.position.extend(list(init_tokens))
    def clone(self):
        # Fast clone: bypass __init__ (no semantic change)
        st = ChemicalState.__new__(ChemicalState)
        st.position = self.position[:]
        return st

    def select_position(self, tok: str):
        self.position.append(tok)


class Node:
    __slots__ = ("position", "parentNode", "childNodes", "wins", "visits", "P")
    def __init__(self, position=None, parent=None):
        self.position = position
        self.parentNode: Optional["Node"] = parent
        self.childNodes: List["Node"] = []
        self.wins = 0.0
        self.visits = 0
        self.P = 0.0  # prior prob

    def select_node(self):
        # UCB selection (same formula), optimized to avoid numpy scalar overhead
        children = self.childNodes
        if not children:
            return None

        parent_visits = self.visits if self.visits > 0 else 1
        nchild = len(children) or 1

        sqrt_parent = math.sqrt(parent_visits)
        inv_nchild = 1.0 / nchild
        c = C_VAL

        best = -1e18
        sel = None

        for ch in children:
            v = ch.visits
            q = (ch.wins / v) if v > 0 else 0.0
            p = ch.P if ch.P > 0.0 else inv_nchild
            u = q + c * p * sqrt_parent / (1.0 + v)
            if u > best:
                best = u
                sel = ch
        return sel

    def add_child(self, token: str):
        n = Node(position=token, parent=self)
        self.childNodes.append(n)
        return n

    def update(self, r: float):
        self.visits += 1
        self.wins += float(r)


# =========================================================
#  seed   (  )
# =========================================================
def build_vocab_and_model():
    print("DEBUG: build_vocab_and_model started")
    zinc_candidates = [
        os.environ.get("SFGDRUG_ZINC_PATH", "").strip(),
        "data/250k_rndm_zinc_drugs_clean.smi",
        "data/250k_rndm_zinc_drugs_clean.smi",
    ]
    zinc_path = None
    for zp in zinc_candidates:
        if zp and os.path.exists(zp):
            zinc_path = zp
            break
    if zinc_path is None:
        raise FileNotFoundError("    ZINC     ,    SFGDRUG_ZINC_PATH")
    print(f"DEBUG: zinc_path = {zinc_path}")

    smile_old = zinc_data_with_bracket_original(zinc_path)
    print("DEBUG: after reading zinc file")

    val, _ = zinc_processed_with_bracket(smile_old)
    print("DEBUG: after processing zinc vocab, VAL size =", len(val))

    print("DEBUG: about to load model from RNN-model/" + MODEL_NAME)
    model = loaded_model("RNN-model/" + MODEL_NAME)
    print("DEBUG: model loaded successfully")
    return val, model


def greedy_tokenize_by_val(smiles: str, val: List[str]) -> Optional[List[str]]:
    """
      VAL          .
       SFG-Drug      token(  [C@H], Cl, Br  ).
    """
    if not smiles:
        return None

    #        ,    [C@H]    '[' 'C' '@' ...
    vocab = sorted(set(val), key=lambda x: (-len(x), x))
    i = 0
    toks: List[str] = []
    n = len(smiles)
    while i < n:
        matched = None
        for t in vocab:
            if smiles.startswith(t, i):
                matched = t
                break
        if matched is None:
            return None
        toks.append(matched)
        i += len(matched)
    return toks


def normalize_seed_string(seed: str) -> List[str]:
    """
       dummy '*'   seed scaffold          .
      :        ,  "    "  .
    """
    cands = []
    s = seed.strip()

    # 1)   
    cands.append(s)

    # 2) [*] -> C
    s2 = s.replace("[*]", "C")
    cands.append(s2)

    # 3)   * -> C
    s3 = s.replace("*", "C")
    cands.append(s3)

    # 4)    *
    s4 = s.replace("[*]", "").replace("*", "")
    cands.append(s4)

    #     
    out = []
    seen = set()
    for x in cands:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def pick_seeds() -> List[Tuple[str, str]]:
    """
       [(seed_name, seed_smiles_or_prefix), ...]
       (  ):
      seed_smiles_override > seed_scaffolds_override > seed_drugs
      seed_fragment_mode=prefer/only  ,        seed_fragment_prefixes(           )
    """
    seeds: List[Tuple[str, str]] = []

    fragment_pool = [str(s) for s in (SEED_FRAGMENT_PREFIXES or []) if str(s).strip()]

    #     :    "   sp3   "        (    /       RDKit   )
    #   :         ,     /     ;         .
    if SEED_FRAGMENT_PREFER_SP3 and fragment_pool:
        def _sp3_score(s: str) -> int:
            ss = str(s).strip()
            #      
            bad_tail = ("(", "=", "#", "/", "\\", "%", "[", "+", "-", "@")
            if ss.endswith(bad_tail):
                return -10
            #        
            for t in ("Cl", "Br", "Si"):
                if ss.endswith(t):
                    return 6
            if ss.endswith(tuple("CNOSPFI")):
                return 6
            if ss.endswith(tuple("cnos")):  #     :      kekulize/     
                return 1
            if ss and ss[-1].isdigit():
                return 0
            return 3

        fragment_pool = sorted(fragment_pool, key=_sp3_score, reverse=True)
    if not fragment_pool and SEED_FRAGMENT_MODE in ("prefer", "only"):
        fragment_pool = list(DEFAULT_FRAGMENT_PREFIXES)

    if SEED_FRAGMENT_MODE == "only" and fragment_pool:
        for i, s in enumerate(fragment_pool, 1):
            seeds.append((f"seed_fragment_{i}", s))
    elif SEED_FRAGMENT_MODE == "prefer" and fragment_pool:
        for i, s in enumerate(fragment_pool, 1):
            seeds.append((f"seed_fragment_{i}", s))
        #           (     )
        if SEED_SMILES_OVERRIDE:
            for i, s in enumerate(SEED_SMILES_OVERRIDE, 1):
                seeds.append((f"seed_smiles_{i}", str(s)))
        elif SEED_SCAFFOLDS_OVERRIDE:
            for i, s in enumerate(SEED_SCAFFOLDS_OVERRIDE, 1):
                seeds.append((f"seed_scaffold_{i}", str(s)))
        elif SEED_DRUGS:
            for i, s in enumerate(SEED_DRUGS, 1):
                seeds.append((f"seed_drug_{i}", str(s)))
    else:
        if SEED_SMILES_OVERRIDE:
            for i, s in enumerate(SEED_SMILES_OVERRIDE, 1):
                seeds.append((f"seed_smiles_{i}", str(s)))
        elif SEED_SCAFFOLDS_OVERRIDE:
            for i, s in enumerate(SEED_SCAFFOLDS_OVERRIDE, 1):
                seeds.append((f"seed_scaffold_{i}", str(s)))
        elif SEED_DRUGS:
            for i, s in enumerate(SEED_DRUGS, 1):
                seeds.append((f"seed_drug_{i}", str(s)))
        else:
            seeds.append(("seed_empty", ""))

    #     (        )
    dedup = []
    seen = set()
    for n, s in seeds:
        key = str(s).strip()
        if not key or key in seen:
            continue
        seen.add(key)
        #     ,         
        if len(key) > SEED_FRAGMENT_MAX_LEN:
            key = key[:SEED_FRAGMENT_MAX_LEN]
        dedup.append((n, key))

    if SEED_MAX_NUM > 0:
        dedup = dedup[:SEED_MAX_NUM]
    return dedup



def _is_balanced_prefix_smiles(s: str) -> bool:
    """
             (       ,            ):
    -   /       "    "
    -       /  /  /        (        )
    -     '%'         (      '%'   '%1')
    """
    if not s:
        return True

    # 1)         
    par = 0
    brk = 0
    for ch in s:
        if ch == '(':
            par += 1
        elif ch == ')':
            par -= 1
            if par < 0:
                return False
        elif ch == '[':
            brk += 1
        elif ch == ']':
            brk -= 1
            if brk < 0:
                return False
    #      par/brk > 0(   )
    # 2)        "    "    
    if s[-1] in ['=', '#', '/', '\\']:
        return False
    # 3)           %   
    if s.endswith('%') or re.search(r'%\d$', s):
        return False
    return True


def _is_good_growth_endpoint_token(tok: str) -> bool:
    """
        :   "       "   ,          /   /  .
    """
    if tok in ('=', '#', '/', '\\', '%'):
        return False
    #    token(     VAL     )
    if tok in ('C','N','O','S','P','F','Cl','Br','I','B','Si','Se','c','n','o','s'):
        return True
    if tok.startswith('[') and tok.endswith(']'):
        return True
    if tok.isdigit():
        return True
    if tok == ')':
        return True
    return False


def _sanitize_seed_prefix(seed_raw: str, val: List[str], max_len_tokens: int = 120) -> Tuple[Optional[List[str]], str]:
    """
       "seed     ":            ,  seed_raw    
    1)    val   
    2)           
    3)    token        
       (tokens or None, used_prefix_str)
    """
    s = str(seed_raw).strip()
    if not s:
        return [], ""

    #      
    toks = greedy_tokenize_by_val(s, val)
    if toks is not None and len(toks) > 0 and len(toks) <= max_len_tokens and _is_balanced_prefix_smiles(s) and _is_good_growth_endpoint_token(toks[-1]):
        return toks, s

    #       (   ),       +      +    token   
    #   "    "   :   prefix             RDKit      .
    for cut in range(len(s)-1, 0, -1):
        ss = s[:cut].strip()
        if not ss:
            break
        if not _is_balanced_prefix_smiles(ss):
            continue
        tt = greedy_tokenize_by_val(ss, val)
        if tt is None or len(tt) == 0:
            continue
        if len(tt) > max_len_tokens:
            continue
        #    token        
        if not _is_good_growth_endpoint_token(tt[-1]):
            continue
        return tt, ss

    return None, s

def resolve_seed_to_tokens(seed_name: str, seed_raw: str, val: List[str]) -> Tuple[Optional[List[str]], str]:
    """
      seed(     / scaffold /      )      SFG-Drug   token   
      : (tokens or None, chosen_seed_string)
    """
    if not seed_raw:
        return [], ""

    seed_raw = str(seed_raw).strip()

    #   /    :  "   "        ,     RDKit     
    if seed_name.startswith("seed_fragment_") or SEED_FRAGMENT_MODE in ("prefer", "only"):
        toks, used = _sanitize_seed_prefix(seed_raw, val, max_len_tokens=SEED_FRAGMENT_MAX_LEN_TOKENS)
        if toks is not None:
            return toks, used
        return None, seed_raw

    #     :          (   * dummy)
    tries = normalize_seed_string(seed_raw)
    for s in tries:
        toks = greedy_tokenize_by_val(s, val)
        if toks is not None and len(toks) > 0:
            return toks, s

    #     :          make_input_smile
    try:
        x = make_input_smile([seed_raw])
        if isinstance(x, list) and len(x) > 0:
            first = x[0]
            if isinstance(first, list) and all(isinstance(t, str) for t in first) and len(first) > 0:
                return first, seed_raw
    except Exception:
        pass

    return None, seed_raw


# =========================================================
#          
# =========================================================
def check_prior_service() -> bool:
    """Lightweight health-check for Token-Mol prior service."""
    try:
        if not PRIOR_URL:
            print("[ERROR] prior_url    (PRIOR_URL=None)")
            return False
        if not P1 or not P2:
            print("[ERROR] prior_protein_path_1 / prior_protein_path_2    ")
            return False

        test_prefix = ["<|beginoftext|>", "<|mask:0|>", "<|mask:0|>", "C"]
        test_candidates = ["C", "N", "O", "(", ")", "[", "]"]

        # reuse the same lightweight candidate filter as in main loop
        filtered = _filter_candidates(test_prefix, test_candidates)
        priors = call_prior(test_prefix, filtered, PRIOR_URL, P1, P2)

        if not isinstance(priors, (list, tuple)) or len(priors) != len(filtered):
            print("[ERROR]           ")
            return False

        # priors should be numeric
        try:
            pmin = float(min(priors))
            pmax = float(max(priors))
        except Exception:
            print("[ERROR]                ")
            return False

        print("[OK]         ")
        print(f"          : {pmin:.6f} - {pmax:.6f}")
        return True

    except Exception as e:
        print(f"[ERROR]         : {e}")
        if DEBUG_TRACEBACK:
            traceback.print_exc()
        return False


# =========================================================
#        (  docking)
# =========================================================
def cheap_prefilter(smi: str) -> Tuple[bool, Dict[str, float]]:
    qed_v = safe_qed(smi)
    sa_v = safe_sa(smi)
    lip_v = lipinski_violations(smi)

    ok = True
    if qed_v < PF_QED_MIN:
        ok = False
    if sa_v > PF_SA_MAX:
        ok = False
    if lip_v > PF_LIP_MAX_VIOL:
        ok = False

    meta = {"qed": qed_v, "sa": sa_v, "lip_viol": float(lip_v)}
    return ok, meta


# =========================================================
#  docking reward
# =========================================================
def docking_to_reward(scores) -> Tuple[float, float, float]:
    """
         :
      [s1, s2, 0] / [s1, s2] /    
      :
      r_dock, s1, s2
    """
    s1, s2 = 1e10, 1e10

    if isinstance(scores, (list, tuple)):
        if len(scores) >= 2:
            s1 = float(scores[0])
            s2 = float(scores[1])
        elif len(scores) == 1:
            s1 = float(scores[0])
            s2 = float(scores[0])
    else:
        s1 = float(scores)
        s2 = float(scores)

    #      (     )
    if abs(s1) > 1e8 or abs(s2) > 1e8:
        return -1.0, s1, s2

    #     ,      
    delta1 = s1 - BASE_VINA_3FAP
    delta2 = s2 - BASE_VINA_7PQV

    # ---     docking reward      ---
    #   :
    #   delta < 0     base   (  );delta > 0     
    #   combined = delta1 + delta2   "    "
    #   worst = max(delta1, delta2)   "  "(    )
    combined = delta1 + delta2
    worst = max(delta1, delta2)
    imbalance = abs(delta1 - delta2)

    def _sat(x: float, k: float) -> float:
        #     :x<0(  )->    ;x>0(  )->    
        return (-(x) * k) / (1.0 + abs(x) * k)

    mode = str(CONF.get("dock_combine_mode", "balanced_dual")).strip().lower()

    #     (   YAML     )
    k_sum = float(CONF.get("dock_k_sum", 0.10))
    k_worst = float(CONF.get("dock_k_worst", 0.10))
    k_imb = float(CONF.get("dock_k_imbalance", 0.10))

    if mode in ("baseline_like", "sum_heavy", "sum-heavy"):
        #    baseline:      (       )
        w_sum = float(CONF.get("dock_w_sum", 0.75))
        w_worst = float(CONF.get("dock_w_worst", 0.25))
        r_sum = _sat(combined, k_sum)
        r_worst = _sat(worst, k_worst)
        r_dock = w_sum * r_sum + w_worst * r_worst

    elif mode in ("balanced_dual", "balanced", "dual"):
        #     :    +    +      
        w_sum = float(CONF.get("dock_w_sum", 0.50))
        w_worst = float(CONF.get("dock_w_worst", 0.50))
        w_imb = float(CONF.get("dock_w_imbalance", 0.30))
        r_sum = _sat(combined, k_sum)
        r_worst = _sat(worst, k_worst)
        r_imb = _sat(imbalance, k_imb)
        r_dock = w_sum * r_sum + w_worst * r_worst - w_imb * r_imb

    elif mode in ("strict_gate", "gate", "threshold"):
        #    :         ,      (   docking       )
        t1 = float(CONF.get("dock_gate_s1", -8.0))
        t2 = float(CONF.get("dock_gate_s2", -8.0))
        fail_value = float(CONF.get("dock_gate_fail_value", -0.20))
        if (s1 > t1) or (s2 > t2):
            r_dock = fail_value
        else:
            #       balanced_dual
            w_sum = float(CONF.get("dock_w_sum", 0.50))
            w_worst = float(CONF.get("dock_w_worst", 0.50))
            w_imb = float(CONF.get("dock_w_imbalance", 0.30))
            r_sum = _sat(combined, k_sum)
            r_worst = _sat(worst, k_worst)
            r_imb = _sat(imbalance, k_imb)
            r_dock = w_sum * r_sum + w_worst * r_worst - w_imb * r_imb
    else:
        #      
        r_sum = _sat(combined, 0.10)
        r_min = _sat(worst, 0.20)
        r_dock = 0.5 * r_sum + 0.5 * r_min

    return float(r_dock), float(s1), float(s2)


# =========================================================
#  MCTS(   seed)
# =========================================================

# ===== Lightweight SMILES constraints (decode-time pruning) =====
def _smiles_balance_state(tokens):
    open_paren = 0
    ring = {}
    last = ""
    for t in tokens:
        if not t or t.startswith("<"):
            continue
        last = t
        if t == "(":
            open_paren += 1
        elif t == ")":
            open_paren = max(0, open_paren - 1)
        elif t.isdigit():
            ring[t] = 1 - ring.get(t, 0)
    open_rings = sum(v for v in ring.values())
    return open_paren, open_rings, last

def _filter_candidates(prefix_tokens, candidates):
    open_paren, open_rings, last = _smiles_balance_state(prefix_tokens)
    out = []
    halogens = {"F", "Cl", "Br", "I"}
    bondish = {"=", "#", "/", "\\", "%"}
    for c in candidates:
        if c is None:
            continue
        if c == "\n" and (open_paren > 0 or open_rings > 0):
            continue
        if c == ")" and open_paren <= 0:
            continue
        if last in halogens and c not in {")", "\n"}:
            continue
        if last in bondish or last == "(":
            if c in {")", "\n"}:
                continue
        out.append(c)
    return out if out else candidates

def mcts_for_seed(
    seed_name: str,
    seed_tokens: List[str],
    seed_smiles_used: str,
    VAL: List[str],
    MODEL,
    global_seen_fps: set,
    global_valid_records: List[Dict[str, Any]],
    max_need: int,
    cfg: Dict[str, Any],
) -> Dict[str, Any]:
    """
        seed        (HOURS_PER_SEED)
         seed       seed_stat
    """
    start_time = time.time()
    seed_start_time = start_time
    seed_records_local: List[Dict[str, Any]] = []
    SAVE_EVERY_MIN = float(cfg.get("save_every_minutes", 10.0))
    SAVE_EVERY_MIN = max(0.5, SAVE_EVERY_MIN)
    last_save = start_time
    last_saved_n = 0

    cnt_depth_guard = 0
    cnt_expand_empty = 0
    cnt_expand_ok = 0
    cnt_prior_call = 0
    cnt_dock_try = 0
    cnt_pushed = 0
    beat_interval = 30.0
    last_beat = start_time
    deadline = start_time + max(1e-6, 3600.0 * HOURS_PER_SEED)

    root_state = ChemicalState(init_tokens=seed_tokens)
    rootnode = Node()

    iter_idx = 0
    dict_id = 1

    print(f"\n===== RUN SEED: {seed_name} =====")
    print(f"[seed] raw/used = {seed_smiles_used}")
    print(f"[seed] token_len = {len(seed_tokens)} | tokens = {seed_tokens[:50]}{'...' if len(seed_tokens) > 50 else ''}")

    while time.time() <= deadline:
        if len(global_valid_records) >= max_need:
            print("[stop] reached target_total_mols")
            break

        iter_idx += 1
        node = rootnode
        state = root_state.clone()


        # ===== loop-level heartbeat (not dependent on selection) =====
        now = time.time()
        if now - last_beat >= beat_interval:
            last_beat = now
            elapsed_m = (now - start_time) / 60.0
            left_m = (HOURS_PER_SEED * 60.0) - elapsed_m
            print(
                f"[heartbeat][{seed_name}] iter={iter_idx} phase=loop "
                f"root_visits={rootnode.visits} elapsed={elapsed_m:.2f}m left={left_m:.2f}m "
                f"total_valid={len(global_valid_records)} unique={len(global_seen_fps)} "
                f"depth_guard={cnt_depth_guard} exp_ok={cnt_expand_ok} exp_empty={cnt_expand_empty} "
                f"prior={cnt_prior_call} dock_try={cnt_dock_try} pushed={cnt_pushed}",
                flush=True
            )
        # ===== periodic checkpoint save (overwrite CSV atomically) =====
        if (now - last_save) >= SAVE_EVERY_MIN * 60.0 and len(global_valid_records) > last_saved_n:
            try:
                save_csv(global_valid_records, OUTPUT_CSV)
                last_saved_n = len(global_valid_records)
                last_save = now
            except Exception as _e:
                print(f"[WARN] periodic save_csv failed: {_e}", flush=True)

        # -------- Selection --------
        while node.childNodes:
            node = node.select_node()
            if node is None:
                break
            state.select_position(node.position)
        if node is None:
            continue
        #     :      /       .
        #   :    token   >=82 ,    token(     )     .
        while (len(state.position) >= 82) or (node.position == "\n"):
            cnt_depth_guard += 1
            if node.parentNode is None or len(state.position) <= 1:
                cur = node
                while cur is not None:
                    cur.update(-1.0)
                    cur = cur.parentNode
                node = None
                break
            node = node.parentNode
            try:
                state.position.pop()
            except Exception:
                state.position = state.position[:-1]

        if node is None:
            continue
        # -------- Expansion --------
        try:
            t0 = time.time()
            expanded_idx = expanded_node(MODEL, state.position, VAL, LOOP_NODE_EXP)
            cnt_expand_ok += 1
            dt = time.time() - t0
            if dt > 2.0:
                print(f"[timing][{seed_name}] expanded_node took {dt:.2f}s (iter={iter_idx})", flush=True)
        except Exception as e:
            print(f"[WARN] expanded_node  : {e}")
            if DEBUG_TRACEBACK:
                traceback.print_exc()
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
            continue
        cnt_expand_ok += 1


        if not expanded_idx:
            cnt_expand_empty += 1
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
            continue

        expanded_idx = expanded_idx[:EXPANDED_K]

        candidates = []
        seen_tok = set()
        for idx in expanded_idx:
            if 0 <= idx < len(VAL):
                t = VAL[idx]
                if t not in seen_tok:
                    seen_tok.add(t)
                    candidates.append(t)

        if not candidates:
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
            continue


        # -------- Candidate sanitization (lightweight grammar/chem guards) --------
        # Goal: improve SMILES validity without retraining; keep main loop logic unchanged.
        pos_tokens = [t for t in state.position if t != "&"]
        last_tok = pos_tokens[-1] if pos_tokens else ""
        paren_bal = pos_tokens.count("(") - pos_tokens.count(")")

        # --- extra SMILES guards (only affects candidate token filtering) ---
        MAX_OPEN_PARENS = int((cfg or {}).get("max_open_parens", 3))
        RING_GUARD = bool((cfg or {}).get("ring_digit_pair_guard", True))
        MAX_OPEN_RINGS = int((cfg or {}).get("max_open_rings", 2))
        ring_cnt = {}
        if RING_GUARD:
            for d in "123456789":
                ring_cnt[d] = pos_tokens.count(d)
            open_rings = sum(1 for d,c in ring_cnt.items() if (c % 2)==1)
        else:
            open_rings = 0

        HALO = {"F", "Cl", "Br", "I"}
        BAD_AFTER_BOND = {"=", "#", "(", ")", "\n", "/", "\\", "%"}
        BAD_TAIL = {"(", "=", "#", "/", "\\", "%", "[", "+", "-", "@"}

        filtered = []
        for tok in candidates:
            # prevent premature end unless balanced and not in a bond-like tail
            if tok == "\n":
                if paren_bal != 0 or (last_tok in BAD_TAIL):
                    continue
                filtered.append(tok)
                continue
            # do not close more than opened
            if tok == ")" and paren_bal <= 0:
                continue
            # limit '(' to avoid runaway branch openings
            if tok == "(" and paren_bal >= MAX_OPEN_PARENS:
                continue

            # simple ring digit pairing guard (optional)
            if RING_GUARD and tok in ring_cnt:
                # avoid placing ring digits right after bond/branch open/end markers
                if last_tok in BAD_AFTER_BOND or last_tok in {"=", "#", "/", "\\"}:
                    continue
                # opening a new ring uses an unused digit; cap number of simultaneously open rings
                if (ring_cnt.get(tok, 0) % 2) == 0:
                    if open_rings >= MAX_OPEN_RINGS:
                        continue
            # after halogens: only allow closing branch or end
            if last_tok in HALO:
                if tok not in {")", "\n"}:
                    continue
            # after bond symbols: next must be an atom (not another bond/paren/end)
            if last_tok in {"=", "#", "/", "\\"}:
                if tok in BAD_AFTER_BOND or tok.isdigit():
                    continue
            # after open paren: avoid immediate close or end or bond direction markers
            if last_tok == "(":
                if tok in {")", "\n", "/", "\\", "%"}:
                    continue
            filtered.append(tok)

        candidates = filtered if filtered else candidates

        # -------- Token-Mol   pocket    --------
        prefix = ["<|beginoftext|>", "<|mask:0|>", "<|mask:0|>"] + [t for t in state.position if t != "&"]
        try:
            cnt_prior_call += 1
            if iter_idx <= 3:
                print(f"[prior_in] iter={iter_idx} prefix_tail={prefix[-12:]} candidates_head={candidates[:12]}", flush=True)
            priors = call_prior(prefix, candidates, PRIOR_URL, P1, P2)
            if not isinstance(priors, (list, tuple)) or len(priors) != len(candidates):
                raise RuntimeError(f"prior len mismatch: {len(priors) if hasattr(priors, '__len__') else 'NA'} vs {len(candidates)}")
            prior_map = {tok: float(p) for tok, p in zip(candidates, priors)}
            if iter_idx <= 10:
                pv = [prior_map.get(t, 0.0) for t in candidates]
                try:
                    import numpy as _np
                    print(f"[prior_stat] iter={iter_idx} min={min(pv):.6g} max={max(pv):.6g} uniq={len(set([round(float(x),6) for x in pv]))}", flush=True)
                except Exception:
                    print(f"[prior_stat] iter={iter_idx} min={min(pv)} max={max(pv)}", flush=True)
                top5 = sorted(candidates, key=lambda t: -prior_map.get(t, 0.0))[:5]
                print("[prior_top5] " + ", ".join([f"{t}:{prior_map.get(t,0.0):.4g}" for t in top5]), flush=True)

        except Exception as e:
            print(f"[ERROR]         (seed={seed_name} iter={iter_idx}): {e}")
            if DEBUG_TRACEBACK:
                traceback.print_exc()
            #      prior,         
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
            continue

        if DEBUG_DIAG:
            top_preview = sorted(candidates, key=lambda t: -prior_map.get(t, 0.0))[:10]
            print(f"[diag] iter={iter_idx} top_candidates_preview={top_preview}")

        # Progressive widening
        max_children = int(1 + (node.visits + 1) ** 0.5)
        quota = max_children - len(node.childNodes)
        allow = max(1, quota)

        # -----------------      prior-guided child sampling -----------------
        #         prior top-N,          .    :
        #   1) temperature     
        #   2) top_k / top_p   
        #   3) Dirichlet noise (  )     
        #   4)       allow   child
        def _sample_children(_cands, _prior_map, _n, _temp, _top_k, _top_p, _dir_alpha, _dir_eps):
            import numpy as _np

            # --- SMILES guards (only 3 rules; avoid broader behavior changes) ---
            # 1) Dynamic mask for ')' when no open '('
            # 2) Upper-bound open '(' count (paren_max_open)
            # 3) Simple ring digit pairing guard (ring_pair_guard; ring_max_open)
            try:
                _pos = [t for t in state.position if t != "&"]
                _banned = set()
                _open_par = _pos.count("(") - _pos.count(")")
                if _open_par <= 0:
                    _banned.add(")")
                _paren_max = int(cfg.get("paren_max_open", 4))
                if _open_par >= _paren_max:
                    _banned.add("(")

                if bool(cfg.get("ring_pair_guard", True)):
                    _ring_max = int(cfg.get("ring_max_open", 3))
                    _last = _pos[-1] if _pos else ""
                    _open_rings = sum((_pos.count(d) % 2) for d in "123456789")
                    if _last.isdigit():
                        _banned.update(list("123456789"))
                    else:
                        if _open_rings >= _ring_max:
                            for d in "123456789":
                                if (_pos.count(d) % 2) == 0:  # would OPEN a new ring
                                    _banned.add(d)

                if _banned:
                    _filtered = [t for t in _cands if t not in _banned]
                    if _filtered:
                        _cands = _filtered
            except Exception:
                pass
            if _n <= 0:
                return []
            ps = _np.array([max(0.0, float(_prior_map.get(t, 0.0))) for t in _cands], dtype=_np.float64)
            if (not _np.isfinite(ps).all()) or ps.sum() <= 0:
                ps = _np.ones(len(_cands), dtype=_np.float64)

            T = max(1e-6, float(_temp))
            ps = _np.power(ps, 1.0 / T)

            if _top_k is not None and int(_top_k) > 0 and int(_top_k) < len(_cands):
                k = int(_top_k)
                idx = _np.argpartition(-ps, k-1)[:k]
                mask = _np.zeros(len(_cands), dtype=bool)
                mask[idx] = True
                ps = ps * mask

            s = ps.sum()
            if s <= 0:
                ps = _np.ones(len(_cands), dtype=_np.float64)
                s = ps.sum()
            ps = ps / s

            if _top_p is not None and 0 < float(_top_p) < 1.0:
                order = _np.argsort(-ps)
                cum = _np.cumsum(ps[order])
                keep_n = int(_np.searchsorted(cum, float(_top_p), side="left")) + 1
                keep_idx = order[:keep_n]
                mask = _np.zeros(len(_cands), dtype=bool)
                mask[keep_idx] = True
                ps = ps * mask
                ps = ps / ps.sum()

            eps = float(_dir_eps) if _dir_eps is not None else 0.0
            if eps > 1e-9:
                alpha = float(_dir_alpha) if _dir_alpha is not None else 0.3
                alpha = max(1e-6, alpha)
                noise = _np.random.dirichlet([alpha] * len(_cands))
                ps = (1.0 - eps) * ps + eps * noise
                ps = ps / ps.sum()

            n = min(int(_n), len(_cands))
            try:
                idx = _np.random.choice(len(_cands), size=n, replace=False, p=ps)
            except Exception:
                idx = _np.random.choice(len(_cands), size=n, replace=False)
            return [_cands[i] for i in idx]

        PRIOR_TEMP = float(cfg.get("prior_temperature", 1.15))
        PRIOR_TOP_K = cfg.get("prior_top_k", 64)
        PRIOR_TOP_P = float(cfg.get("prior_top_p", 0.95))
        DIR_ALPHA = float(cfg.get("prior_dirichlet_alpha", 0.30))
        DIR_EPS = float(cfg.get("prior_dirichlet_eps", 0.15))

        sampled = _sample_children(
            candidates, prior_map, allow,
            PRIOR_TEMP, PRIOR_TOP_K, PRIOR_TOP_P,
            DIR_ALPHA, DIR_EPS
        )
        for tok in sampled:
            ch = node.add_child(tok)
            ch.P = prior_map.get(tok, 1.0 / max(1, len(candidates)))

        # -------- Simulation / Rollout(RNN)--------
        new_compounds: List[str] = []
        nodeadded_all: List[str] = []
        try:
            for _ in range(SIM_NUM):
                nodeadded = node_to_add(expanded_idx, VAL)
                all_possible = chem_kn_simulation(MODEL, state.position, VAL, nodeadded)
                gen = predict_smile(all_possible, VAL)
                new_compounds.extend(make_input_smile(gen))
                nodeadded_all.extend(nodeadded)
        except Exception as e:
            print(f"[WARN] rollout  : {e}")
            if DEBUG_TRACEBACK:
                traceback.print_exc()
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
            continue

        # -------- check_node_type(     )--------
        generated_dict = {}
        use_radical = (False if SMOKE_BYPASS_RADICAL else RADICAL_CHECK)
        use_hashimoto = (False if SMOKE_BYPASS_HASHIMOTO else HASHIMOTO_FILTER)

        try:
            node_index, vinadock_score_list, valid_smile, generated_dict = check_node_type(
                new_compounds, generated_dict,
                sa_threshold=SA_THRESHOLD, rule=RULE5,
                radical=use_radical, hashimoto_filter=use_hashimoto,
                dict_id=dict_id, trial=TRIAL
            )
        except Exception as e:
            print(f"[WARN] check_node_type  : {e}")
            if DEBUG_TRACEBACK:
                traceback.print_exc()
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
            continue

        dict_id += 1

        # --------        + docking(      rdock_test_MP)--------
        #   :       check_node_type     dock  ;       ,      rdock_test_MP         
        #       check_node_type         ,      "   /    "
        candidates_for_dock = []
        for smi in valid_smile:
            ok, meta = cheap_prefilter(smi)
            if ok:
                candidates_for_dock.append((smi, meta))

        #      docking   (   )
        if len(candidates_for_dock) > DOCK_MAX_PER_ROLLOUT:
            #     QED    (     prior/   )
            candidates_for_dock = sorted(candidates_for_dock, key=lambda x: x[1]["qed"], reverse=True)[:DOCK_MAX_PER_ROLLOUT]

        #    atom -> child     (      )
        child_by_atom = {}
        for c in node.childNodes:
            child_by_atom[c.position] = c

        rewards_for_child: List[Tuple[Node, float]] = []
        n_docked = 0

        cnt_dock_try += 1

        for smi, meta in candidates_for_dock:
            if len(global_valid_records) >= max_need:
                break

            smi_eval = sanitize_smiles_for_eval(smi)
            if not smi_eval:
                continue
            #            ,      docking(     +    unique)
            if cfg.get("skip_dock_if_duplicate", True):
                try:
                    fpb0 = morgan_fp_bytes(smi_eval)
                    if fpb0 in global_seen_fps:
                        continue
                except Exception:
                    pass


            try:
                scores = vinadock_score(smi_eval)  #   :       
            except Exception as e:
                print(f"[WARN] vinadock_score  : {smi} -> {e}")
                if DEBUG_TRACEBACK:
                    traceback.print_exc()
                continue

            r_dock, s1, s2 = docking_to_reward(scores)
            if n_docked < 3:
                print(f"[dock] iter={iter_idx} smi={smi_eval} s1={s1:.3f} s2={s2:.3f} r_dock={r_dock:.3f} qed={meta.get('qed', float('nan')):.3f} sa={meta.get('sa', float('nan')):.3f}", flush=True)

            qed_v = meta["qed"]
            sa_v = meta["sa"]
            lip_v = int(meta["lip_viol"])

            sa_term = max(0.0, 1.0 - (sa_v - 2.5) / 3.5)
            lip_pen = float(max(0, lip_v - PF_LIP_MAX_VIOL))

            fpb = morgan_fp_bytes(smi_eval)
            is_novel = fpb not in global_seen_fps
            if is_novel:
                global_seen_fps.add(fpb)

            r = (
                W_DOCK * r_dock
                + W_QED * qed_v
                + W_SA * sa_term
                - W_LIP * lip_pen
                + (DIVERSITY_BONUS if is_novel else 0.0)
            )

            #     
            rec = {
                "seed_name": seed_name,
                "seed_smiles_used": seed_smiles_used,
                "smiles": smi_eval,
                "dock_3fap": s1,
                "dock_7pqv": s2,
                "dock_sum": (s1 + s2) if (abs(s1) < 1e8 and abs(s2) < 1e8) else float("nan"),
                "reward": float(r),
                "r_dock": float(r_dock),
                "qed": float(qed_v),
                "sa": float(sa_v),
                "lip_viol": int(lip_v),
                "is_novel": int(is_novel),
                "depth": len(state.position),
                "iter_idx": iter_idx,
                "elapsed_sec": round(time.time() - start_time, 3),
            }
            cnt_pushed += 1
            global_valid_records.append(rec)
            append_raw_record(RAW_OUTPUT_CSV, rec)
            seed_records_local.append(rec)
            n_docked += 1

            #             child(     :    nodeadded_all   token)
            #         ,        node
            assigned = False
            for atom in nodeadded_all:
                ch = child_by_atom.get(atom, None)
                if ch is not None:
                    rewards_for_child.append((ch, r))
                    assigned = True
                    break
            if not assigned:
                rewards_for_child.append((node, r))  # fallback

        # -------- Backprop --------
        if rewards_for_child:
            for target_node, rr in rewards_for_child:
                cur = target_node
                while cur is not None:
                    cur.update(rr)
                    cur = cur.parentNode
        else:
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode

        if iter_idx % 1 == 0:
            elapsed_m = (time.time() - start_time) / 60.0
            left_m = max((deadline - time.time()) / 60.0, 0.0)
            print(
                f"[progress][{seed_name}] iter={iter_idx} "
                f"elapsed={elapsed_m:.2f}m left={left_m:.2f}m "
                f"docked={n_docked} total_valid={len(global_valid_records)} unique={len(global_seen_fps)}"
            )

    seed_stat = make_seed_stat(
        seed_name=seed_name,
        seed_smiles_used=seed_smiles_used,
        seed_records=seed_records_local,
        elapsed_sec=(time.time() - seed_start_time),
    )
    return seed_stat
# =========================================================
#      
# =========================================================
RAW_COLS = [
    "seed_name", "seed_smiles_used", "smiles",
    "dock_3fap", "dock_7pqv", "dock_sum",
    "reward", "r_dock", "qed", "sa", "lip_viol", "is_novel",
    "depth", "iter_idx", "elapsed_sec",
]

SEED_STATS_COLS = [
    "seed_name", "seed_smiles_used",
    "n_valid", "n_unique", "mean_dock_3fap", "mean_dock_7pqv", "mean_dock_sum",
    "mean_reward", "mean_qed", "mean_sa", "elapsed_min",
]


def ensure_csv_header(csv_path: str, cols: List[str]):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True) if os.path.dirname(csv_path) else None
    if not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=cols)
            writer.writeheader()


def append_csv_row(csv_path: str, row: Dict[str, Any], cols: List[str]):
    ensure_csv_header(csv_path, cols)
    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writerow({k: row.get(k, "") for k in cols})
        f.flush()
        os.fsync(f.fileno())


def append_raw_record(raw_csv: str, rec: Dict[str, Any]):
    append_csv_row(raw_csv, rec, RAW_COLS)


def make_seed_stat(seed_name: str, seed_smiles_used: str, seed_records: List[Dict[str, Any]], elapsed_sec: float) -> Dict[str, Any]:
    if not seed_records:
        return {
            "seed_name": seed_name,
            "seed_smiles_used": seed_smiles_used,
            "n_valid": 0,
            "n_unique": 0,
            "mean_dock_3fap": "",
            "mean_dock_7pqv": "",
            "mean_dock_sum": "",
            "mean_reward": "",
            "mean_qed": "",
            "mean_sa": "",
            "elapsed_min": round(elapsed_sec / 60.0, 3),
        }

    n = len(seed_records)
    uniq = len({r.get("smiles", "") for r in seed_records if r.get("smiles", "")})
    mean = lambda key: round(sum(float(r.get(key, 0.0)) for r in seed_records) / max(n, 1), 6)
    return {
        "seed_name": seed_name,
        "seed_smiles_used": seed_smiles_used,
        "n_valid": n,
        "n_unique": uniq,
        "mean_dock_3fap": mean("dock_3fap"),
        "mean_dock_7pqv": mean("dock_7pqv"),
        "mean_dock_sum": mean("dock_sum"),
        "mean_reward": mean("reward"),
        "mean_qed": mean("qed"),
        "mean_sa": mean("sa"),
        "elapsed_min": round(elapsed_sec / 60.0, 3),
    }
def save_csv(records: List[Dict[str, Any]], out_csv: str):
    import csv
    os.makedirs(os.path.dirname(out_csv), exist_ok=True) if os.path.dirname(out_csv) else None

    cols = [
        "seed_name", "seed_smiles_used", "smiles",
        "dock_3fap", "dock_7pqv", "dock_sum",
        "reward", "r_dock", "qed", "sa", "lip_viol", "is_novel",
        "depth", "iter_idx", "elapsed_sec",
    ]
    with open(out_csv + ".tmp", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in records:
            row = {k: r.get(k, "") for k in cols}
            w.writerow(row)
    os.replace(out_csv + ".tmp", out_csv)
    print(f"[save] results -> {out_csv}  (n={len(records)})")


def save_seedmap(seed_rows: List[Dict[str, Any]], out_csv: str):
    import csv
    if not out_csv:
        return
    os.makedirs(os.path.dirname(out_csv), exist_ok=True) if os.path.dirname(out_csv) else None

    cols = ["seed_name", "seed_raw", "seed_used", "tokenizable", "token_len", "note"]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in seed_rows:
            w.writerow({k: r.get(k, "") for k in cols})
    print(f"[save] seedmap -> {out_csv}  (n={len(seed_rows)})")


# =========================================================
#     
# =========================================================
def main():
    print("DEBUG: main started")
    VAL, MODEL = build_vocab_and_model()
    cfg = CONF  # use global loaded config

    print("     Token-Mol     ...")
    if not check_prior_service():
        print("[ERROR]        ,  ")
        sys.exit(1)

    seeds = pick_seeds()
    if not seeds:
        print("[ERROR]      seed")
        sys.exit(1)

    seed_rows = []
    resolved_seeds: List[Tuple[str, List[str], str]] = []

    for seed_name, seed_raw in seeds:
        seed_norm = normalize_seed_prefix(seed_raw)
        toks, used_smi = resolve_seed_to_tokens(seed_name, seed_norm, VAL)
        if toks is None:
            print(f"[seed][SKIP] {seed_name}     : {seed_raw}")
            seed_rows.append({
                "seed_name": seed_name,
                "seed_raw": seed_raw,
                "seed_used": used_smi,
                "tokenizable": 0,
                "token_len": 0,
                "note": "tokenize_failed",
            })
            continue

        print(f"[seed][OK] {seed_name} token_len={len(toks)} used={used_smi}")
        seed_rows.append({
            "seed_name": seed_name,
            "seed_raw": seed_raw,
            "seed_used": used_smi,
            "tokenizable": 1,
            "token_len": len(toks),
            "note": ("norm_halogen_open" if seed_norm != seed_raw else ""),
        })
        resolved_seeds.append((seed_name, toks, used_smi))

    save_seedmap(seed_rows, SEEDMAP_CSV)

    if not resolved_seeds:
        print("[ERROR]    seed      ,  ")
        sys.exit(1)

    #     
    global_seen_fps = set()
    global_valid_records: List[Dict[str, Any]] = []
    ensure_csv_header(RAW_OUTPUT_CSV, RAW_COLS)
    ensure_csv_header(SEED_STATS_CSV, SEED_STATS_COLS)
    print("     Seeded MCTS + Token-Mol   pocket     ")
    for seed_name, seed_toks, seed_used in resolved_seeds:
        if len(global_valid_records) >= TARGET_TOTAL_MOLS:
            break
        seed_stat = mcts_for_seed(
            seed_name=seed_name,
            seed_tokens=seed_toks,
            seed_smiles_used=seed_used,
            VAL=VAL,
            MODEL=MODEL,
            global_seen_fps=global_seen_fps,
            global_valid_records=global_valid_records,
            max_need=TARGET_TOTAL_MOLS,
            cfg=cfg,
        )
        append_csv_row(SEED_STATS_CSV, seed_stat, SEED_STATS_COLS)
        save_csv(global_valid_records, OUTPUT_CSV)
    #     (  dock_sum     )
    def _sort_key(r):
        s = r.get("dock_sum", float("inf"))
        try:
            if math.isnan(float(s)):
                return (1, float("inf"))
            return (0, float(s))
        except Exception:
            return (1, float("inf"))

    global_valid_records = sorted(global_valid_records, key=_sort_key)
    save_csv(global_valid_records, OUTPUT_CSV)

    print("[OK]   ")
    print("total_valid_records:", len(global_valid_records))
    if global_valid_records:
        print("top1:", global_valid_records[0])


if __name__ == "__main__":
    main()