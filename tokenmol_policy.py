# tokenmol_policy.py
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import List, Tuple, Optional

from bert_tokenizer import ExpressionBertTokenizer
from ada_model import Token3D
from pocket_fine_tuning_rmse import Ada_config


class TokenMolProteinPolicy(nn.Module):
    """
    真正 protein-aware 的 Token-Mol 策略网络：
    - 使用 Token3D 结构
    - forward(x, protein_matrix) 里显式融合 pocket 表示
    """
    def __init__(
        self,
        vocab_path: str,
        pretrain_dir: str,
        model_ckpt: str,
        device: Optional[torch.device] = None,
    ):
        super().__init__()
        self.device = device or torch.device(
            "cuda" if torch.cuda.is_available() else "cpu"
        )

        # 和原项目 gen.py 一致的 tokenizer
        self.tokenizer = ExpressionBertTokenizer.from_pretrained(vocab_path)

        # 和 pocket_fine_tuning_rmse.py 一致的 Token3D
        self.model = Token3D(pretrain_path=pretrain_dir, config=Ada_config)

        # 兼容 DataParallel 的 state_dict
        state = torch.load(model_ckpt, map_location=self.device)
        if isinstance(state, dict) and "state_dict" in state:
            state = state["state_dict"]
        state = {k.replace("module.", ""): v for k, v in state.items()}

        missing, unexpected = self.model.load_state_dict(state, strict=False)
        if missing:
            print("[TokenMolProteinPolicy] Missing keys:", missing)
        if unexpected:
            print("[TokenMolProteinPolicy] Unexpected keys:", unexpected)

        self.model.to(self.device).eval()

    # --------- 编码 / 计算 logits 的小工具 ---------
    def encode_tokens(self, tokens: List[str]) -> torch.Tensor:
        """
        直接用 convert_tokens_to_ids 保证和原项目一致，不走自动分词。
        """
        ids = self.tokenizer.convert_tokens_to_ids(tokens)
        return torch.tensor([ids], dtype=torch.long, device=self.device)

    def step_logits(
        self,
        prefix_tokens: List[str],
        protein_tensor: torch.Tensor,
    ) -> torch.Tensor:
        """
        给定 token 前缀和一个蛋白 pocket 表示，返回最后一个位置的 logits
        protein_tensor: shape [L, D] 或 [1, L, D]
        返回: logits_last (vocab_size,)
        """
        input_ids = self.encode_tokens(prefix_tokens)  # [1, T]

        if protein_tensor.dim() == 2:
            protein_batch = protein_tensor.unsqueeze(0)  # [1, L, D]
        else:
            protein_batch = protein_tensor  # [1, L, D]

        protein_batch = protein_batch.to(self.device)

        with torch.no_grad():
            outputs = self.model(input_ids, protein_batch)
            logits = outputs.logits  # [1, T, V]
            last_logits = logits[:, -1, :].squeeze(0)  # [V]

        return last_logits

    # --------- 计算候选 token 的先验概率 ---------
    def prior_for_candidates(
        self,
        prefix_tokens: List[str],
        candidates: List[str],
        protein_tensor: torch.Tensor,
        temperature: float = 1.0,
    ) -> List[float]:
        logits = self.step_logits(prefix_tokens, protein_tensor) / max(
            temperature, 1e-6
        )

        # softmax 得到全 vocab 概率
        probs = F.softmax(logits, dim=-1)

        # 抽取我们关心的候选 token 概率
        cand_probs = []
        for tok in candidates:
            tok_id = self.tokenizer.convert_tokens_to_ids(tok)
            if tok_id is None:
                cand_probs.append(0.0)
            else:
                cand_probs.append(float(probs[tok_id]))

        # 再归一化一遍，以防候选子集中总和 < 1
        s = sum(cand_probs)
        if s > 0:
            cand_probs = [p / s for p in cand_probs]
        else:
            # 极端情况全 0，回退成均匀分布
            cand_probs = [1.0 / len(candidates)] * len(candidates)

        return cand_probs


def load_tokenmol_model(
    model_ckpt: str,
    vocab_path: str,
    pretrain_dir: str,
) -> Tuple[TokenMolProteinPolicy, ExpressionBertTokenizer, torch.device]:
    """
    供 tokenmol_prior_server.py 调用的加载函数，接口保持不变。
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    policy = TokenMolProteinPolicy(
        vocab_path=vocab_path,
        pretrain_dir=pretrain_dir,
        model_ckpt=model_ckpt,
        device=device,
    )
    return policy, policy.tokenizer, device

