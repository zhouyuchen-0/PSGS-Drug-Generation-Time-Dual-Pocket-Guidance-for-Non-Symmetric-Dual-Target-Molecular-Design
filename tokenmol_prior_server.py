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

# ----------------- 配置（根据你服务器路径适当修改） -----------------
MODEL_CKPT = "./Trained_model/pocket_generation.pt"
VOCAB_PATH = "./data/torsion_version/torsion_voc_pocket.csv"
PRETRAIN_DIR = "./Pretrained_model"

app = FastAPI(title="Token-Mol Dual-Pocket Prior Server")

print("加载 Token-Mol protein-aware 模型...")
POLICY, TOKENIZER, DEVICE = load_tokenmol_model(
    MODEL_CKPT,
    VOCAB_PATH,
    PRETRAIN_DIR,
)
print("✅ 模型加载完成，设备:", DEVICE)

# pocket 缓存，避免反复从磁盘读
_pocket_cache: Dict[str, torch.Tensor] = {}



def load_pocket(path: str) -> torch.Tensor:
    """
    加载单个 pocket 的表示:
    兼容几种可能的 .pkl 结构：
      - 直接是 np.ndarray，形状 [L, D] 或 [1, L, D]
      - list / tuple：取第一个元素，如果里面还是 ndarray 再取
      - dict：尝试用 'pocket' / 'protein' / 'data' 这些 key，找不到就取第一个 value
    """
    if path in _pocket_cache:
        return _pocket_cache[path]

    if not os.path.exists(path):
        raise FileNotFoundError(f"Pocket 文件不存在: {path}")

    with open(path, "rb") as f:
        obj = pickle.load(f)

    pocket = None

    # 情况 1：直接是 ndarray
    if isinstance(obj, np.ndarray):
        pocket = obj

    # 情况 2：list 或 tuple（例如 [ndarray, ndarray, ...]）
    elif isinstance(obj, (list, tuple)):
        if len(obj) == 0:
            raise ValueError(f"Pocket 文件 {path} 为空 list。")
        first = obj[0]
        if isinstance(first, np.ndarray):
            pocket = first
        elif isinstance(first, (list, tuple)):
            pocket = np.array(first, dtype=np.float32)
        else:
            # 打印一下类型，便于调试
            print(f"[load_pocket] {path} list[0] 类型为 {type(first)}，尝试转 np.array")
            pocket = np.array(first, dtype=np.float32)

    # 情况 3：dict（有些脚本会存成 {'pocket': ndarray, ...}）
    elif isinstance(obj, dict):
        for key in ["pocket", "protein", "data", "feature", "pocket_feature"]:
            if key in obj:
                pocket = obj[key]
                print(f"[load_pocket] 从 dict 使用 key='{key}' 作为 pocket")
                break
        if pocket is None:
            # 退而求其次：取第一个 value
            first_key = next(iter(obj.keys()))
            pocket = obj[first_key]
            print(f"[load_pocket] dict 无标准 key，取第一个 key='{first_key}' 作为 pocket")

    else:
        # 其他奇怪类型，直接尝试转成 np.array
        print(f"[load_pocket] {path} 类型为 {type(obj)}，尝试 np.array 转换")
        pocket = np.array(obj, dtype=np.float32)

    # 最终检查
    if pocket is None:
        raise ValueError(f"无法从 {path} 解析 pocket 数据，类型为 {type(obj)}")

    pocket = np.asarray(pocket, dtype=np.float32)
    pocket_tensor = torch.as_tensor(pocket, dtype=torch.float32, device=DEVICE)

    # 期望最终是 [1, L, D]
    if pocket_tensor.dim() == 2:
        pocket_tensor = pocket_tensor.unsqueeze(0)
    elif pocket_tensor.dim() == 3:
        # 假设已经是 [1, L, D] 或 [N, L, D]，我们只用第一个 pocket
        if pocket_tensor.size(0) > 1:
            pocket_tensor = pocket_tensor[0:1]
    else:
        raise ValueError(
            f"Pocket tensor 维度异常: {pocket_tensor.dim()}，shape={tuple(pocket_tensor.shape)}"
        )

    _pocket_cache[path] = pocket_tensor
    print(f"✅ 加载 pocket: {path}, shape={tuple(pocket_tensor.shape)}")
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
    给定：
    - prefix_tokens: token 序列
    - candidates: 候选下一个 token
    - protein_path: pocket.pkl 的路径

    返回：
    - prior: 与 candidates 等长的概率列表
    - status: "ok" 或 "error"
    - message: 如果出错，附带错误信息
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
        # 在服务器终端打印完整 traceback，方便你看
        print("[/prior] Error when handling protein_path =", req.protein_path)
        traceback.print_exc()

        # 返回均匀分布 + 错误提示（HTTP 200，但 status=error）
        n = len(req.candidates)
        prior = [1.0 / max(n, 1)] * n
        return {
            "status": "error",
            "prior": prior,
            "message": str(e),
        }


