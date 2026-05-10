#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FastAPI server for Token-Mol pocket-prior inference.

The server exposes /health and /prior endpoints used by the PSGS-Drug generation scripts.
"""

# tokenmol_prior_server.py
import os
import pickle
from typing import List, Dict
import numpy as np
import torch
from fastapi import FastAPI
from pydantic import BaseModel

from tokenmol_policy import load_tokenmol_model

import traceback

# -----------------   (            ) -----------------
MODEL_CKPT = os.environ.get("TOKENMOL_MODEL_CKPT", "./Trained_model/pocket_generation.pt")
VOCAB_PATH = os.environ.get("TOKENMOL_VOCAB_PATH", "./data/torsion_version/torsion_voc_pocket.csv")
PRETRAIN_DIR = os.environ.get("TOKENMOL_PRETRAIN_DIR", "./Pretrained_model")

app = FastAPI(title="Token-Mol Dual-Pocket Prior Server")

print("   Token-Mol protein-aware   ...")
POLICY, TOKENIZER, DEVICE = load_tokenmol_model(
    MODEL_CKPT,
    VOCAB_PATH,
    PRETRAIN_DIR,
)
print("[OK]       ,  :", DEVICE)

# pocket   ,        
_pocket_cache: Dict[str, torch.Tensor] = {}



def load_pocket(path: str) -> torch.Tensor:
    """
         pocket    :
            .pkl   :
      -     np.ndarray,   [L, D]   [1, L, D]
      - list / tuple:      ,       ndarray   
      - dict:    'pocket' / 'protein' / 'data'    key,         value
    """
    if path in _pocket_cache:
        return _pocket_cache[path]

    if not os.path.exists(path):
        raise FileNotFoundError(f"Pocket      : {path}")

    with open(path, "rb") as f:
        obj = pickle.load(f)

    pocket = None

    #    1:    ndarray
    if isinstance(obj, np.ndarray):
        pocket = obj

    #    2:list   tuple(   [ndarray, ndarray, ...])
    elif isinstance(obj, (list, tuple)):
        if len(obj) == 0:
            raise ValueError(f"Pocket    {path}    list.")
        first = obj[0]
        if isinstance(first, np.ndarray):
            pocket = first
        elif isinstance(first, (list, tuple)):
            pocket = np.array(first, dtype=np.float32)
        else:
            #       ,    
            print(f"[load_pocket] {path} list[0]     {type(first)},    np.array")
            pocket = np.array(first, dtype=np.float32)

    #    3:dict(        {'pocket': ndarray, ...})
    elif isinstance(obj, dict):
        for key in ["pocket", "protein", "data", "feature", "pocket_feature"]:
            if key in obj:
                pocket = obj[key]
                print(f"[load_pocket]   dict    key='{key}'    pocket")
                break
        if pocket is None:
            #      :     value
            first_key = next(iter(obj.keys()))
            pocket = obj[first_key]
            print(f"[load_pocket] dict     key,     key='{first_key}'    pocket")

    else:
        #       ,       np.array
        print(f"[load_pocket] {path}     {type(obj)},   np.array   ")
        pocket = np.array(obj, dtype=np.float32)

    #     
    if pocket is None:
        raise ValueError(f"    {path}    pocket   ,    {type(obj)}")

    pocket = np.asarray(pocket, dtype=np.float32)
    pocket_tensor = torch.as_tensor(pocket, dtype=torch.float32, device=DEVICE)

    #       [1, L, D]
    if pocket_tensor.dim() == 2:
        pocket_tensor = pocket_tensor.unsqueeze(0)
    elif pocket_tensor.dim() == 3:
        #       [1, L, D]   [N, L, D],        pocket
        if pocket_tensor.size(0) > 1:
            pocket_tensor = pocket_tensor[0:1]
    else:
        raise ValueError(
            f"Pocket tensor     : {pocket_tensor.dim()},shape={tuple(pocket_tensor.shape)}"
        )

    _pocket_cache[path] = pocket_tensor
    print(f"[OK]    pocket: {path}, shape={tuple(pocket_tensor.shape)}")
    return pocket_tensor



class PriorRequest(BaseModel):
    prefix_tokens: List[str]
    candidates: List[str]
    protein_path: str


class HealthResponse(BaseModel):
    status: str
    models_loaded: bool
    pockets_loaded: bool


@app.get("/health", response_model=HealthResponse)
def health():
    return HealthResponse(
        status="ok",
        models_loaded=True,
        pockets_loaded=True,
    )


@app.post("/prior")
def prior(req: PriorRequest):
    """
      :
    - prefix_tokens: token   
    - candidates:       token
    - protein_path: pocket.pkl    

      :
    - prior:   candidates        
    - status: "ok"   "error"
    - message:     ,      
    """
    try:
        pocket_tensor = load_pocket(req.protein_path)

        prior_probs = POLICY.prior_for_candidates(
            prefix_tokens=req.prefix_tokens,
            candidates=req.candidates,
            protein_tensor=pocket_tensor,
        )

        return {
            "status": "ok",
            "prior": prior_probs,
            "message": "",
        }

    except Exception as e:
        #            traceback,    
        print("[/prior] Error when handling protein_path =", req.protein_path)
        traceback.print_exc()

        #        +     (HTTP 200,  status=error)
        n = len(req.candidates)
        prior = [1.0 / max(n, 1)] * n
        return {
            "status": "error",
            "prior": prior,
            "message": str(e),
        }


