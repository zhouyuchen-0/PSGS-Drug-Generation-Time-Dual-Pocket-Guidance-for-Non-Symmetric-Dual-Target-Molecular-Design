# sbmolgen.py  （改进版：双靶 + Token-Mol 双 pocket 先验 + 双靶点结果输出）
# 说明：
# 1）Token-Mol 先验通过 prior_client.call_prior() 完成
# 2）prior_client.py 放在 SFG-Drug 根目录（与本文件同一目录）
# 3）Token-Mol 服务器端只跑 tokenmol_prior_server.py 即可

import os
import time
import yaml
from math import sqrt
from typing import Any, Dict, List

import numpy as np
from rdkit import Chem
from rdkit.Chem import QED, rdMolDescriptors

from utils.add_node_type_zinc import (
    predict_smile, make_input_smile, expanded_node,
    node_to_add, chem_kn_simulation, check_node_type
)
from utils.load_model import loaded_model
from utils.make_smile import zinc_processed_with_bracket, zinc_data_with_bracket_original

# ⭐⭐ 导入 prior_client 中的函数和异常类 ⭐⭐
try:
    from prior_client import call_prior, PriorServiceError
except ImportError as e:
    print(f"❌ 导入 prior_client 失败: {e}")
    print("❌ 请确保 prior_client.py 文件在正确的位置")
    import sys
    sys.exit(1)


# =========================================================
#  RDKit 辅助
# =========================================================
def safe_qed(smi: str) -> float:
    try:
        m = Chem.MolFromSmiles(smi)
        return float(QED.qed(m)) if m else 0.0
    except Exception:
        return 0.0


def safe_sa(smi: str) -> float:
    """简化 SA 打分，与原项目风格一致。"""
    try:
        m = Chem.MolFromSmiles(smi)
        if not m:
            return 6.0
        heavy = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1)
        rings = Chem.GetSSSR(m)
        sa = 2.0 + 0.2 * max(0, heavy - 20) + 0.2 * max(0, rings - 2)
        return float(min(max(sa, 1.0), 10.0))
    except Exception:
        return 6.0


def morgan_fp(smi: str, nbits: int = 2048):
    """把分子转成 bit string，用于多样性判断。"""
    try:
        m = Chem.MolFromSmiles(smi)
        if not m:
            return smi.encode("utf-8")
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius=2, nBits=nbits)
        return fp.ToBitString().encode("ascii")
    except Exception:
        return smi.encode("utf-8")


# =========================================================
#  读取配置（双靶 + 先验）
# =========================================================
def load_conf() -> Dict[str, Any]:
    cfg_path = r"E:\graduate\meeting\mwx\SFG-Drug-code\SFG-Drug\example_ligand_design\setting_prior.yaml"
    with open(cfg_path, "r", encoding="utf-8") as f:
        return yaml.load(f, Loader=yaml.SafeLoader) or {}


conf = load_conf()
TARGET = conf.get("target", "3fap-7pqv-prior")
trial = conf.get("trial", 1)
LOOP_NODE_EXP = conf.get("loop_num_nodeExpansion", 82)
HOURS = conf.get("hours", 12)

SA_THRESHOLD = conf.get("sa_threshold", 5.0)
RULE5 = conf.get("rule5", 1)
RADICAL_CHECK = conf.get("radical_check", True)
HASHIMOTO_FILTER = conf.get("hashimoto_filter", True)
SIM_NUM = conf.get("simulation_num", 3)
MODEL_NAME = conf.get("model_name", "model")
BASE_VINA = conf.get("base_vinadock_score", -7.0)

MCONF = conf.get("mcts", {})
C_VAL = float(MCONF.get("c_val", 1.3))
EXPANDED_K = int(MCONF.get("expanded_k", 24))
DIVERSITY_BONUS = float(MCONF.get("diversity_bonus", 0.05))

RCONF = conf.get("reward", {})
W_DOCK = RCONF.get("w_docking", 0.70)
W_QED = RCONF.get("w_qed", 0.20)
W_SA = RCONF.get("w_sa", 0.10)

# 先验相关（只用于打印，实际逻辑在 prior_client 里）
PRIOR_URL = conf.get("prior_url", None)
P1 = conf.get("prior_protein_path_1", None)
P2 = conf.get("prior_protein_path_2", None)

print("====== PRIOR CONFIG ======")
print("target:", TARGET)
print("trial:", trial)
print("c_val:", C_VAL, "expanded_k:", EXPANDED_K, "hours:", HOURS)
print("prior_url:", PRIOR_URL)
print("pocket1:", P1)
print("pocket2:", P2)
print("==========================")

OUTPUT_PATH = f"result_{TARGET}_C{C_VAL}_trial{trial}.txt"
with open(OUTPUT_PATH, "w", encoding="utf-8") as g:
    # ⭐ 修改：结果文件头改为双靶点格式，一行一个分子
    g.write("smiles,dock_3fap,dock_7pqv,dock_combined,depth,used_time\n")


# =========================================================
#  状态与节点定义
# =========================================================
class ChemicalState:
    def __init__(self):
        self.position = ["&"]

    def clone(self):
        st = ChemicalState()
        st.position = self.position[:]
        return st

    def select_position(self, tok):
        self.position.append(tok)


class Node:
    def __init__(self, position=None, parent=None, state=None):
        self.position = position
        self.parentNode = parent
        self.childNodes: List["Node"] = []
        self.wins = 0.0
        self.visits = 0
        self.P = 0.0  # 先验概率

    def is_root(self):
        return self.parentNode is None

    def select_node(self):
        best, sel = -1e9, None
        for ch in self.childNodes:
            Q = (ch.wins / ch.visits) if ch.visits > 0 else 0.0
            prior = ch.P if ch.P > 0 else (1.0 / max(1, len(self.childNodes)))
            U = Q + C_VAL * prior * np.sqrt(max(1, self.visits)) / (1 + ch.visits)
            if U > best:
                best, sel = U, ch
        return sel

    def add_child(self, token, state):
        n = Node(position=token, parent=self, state=state)
        self.childNodes.append(n)
        return n

    def update(self, r):
        self.visits += 1
        self.wins += float(r)


# =========================================================
#  MCTS 主流程（这里调用 Token-Mol 双 pocket 先验）
# =========================================================
def mcts(root: ChemicalState):
    start_time = time.time()
    deadline = start_time + 3600 * HOURS

    rootnode = Node(state=root)

    # 多样性缓存
    seen_fps = set()

    dict_id = 1
    valid_compounds: List[str] = []
    out_f = open(OUTPUT_PATH, "a", encoding="utf-8")

    while time.time() <= deadline:
        node = rootnode
        state = root.clone()

        # -------- Selection --------
        while node.childNodes:
            node = node.select_node()
            state.select_position(node.position)

        # 深度 / 终止判断
        if len(state.position) >= 70 or node.position == "\n":
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
            continue

        # -------- Expansion：得到一批候选 token --------
        expanded_idx = expanded_node(MODEL, state.position, VAL, LOOP_NODE_EXP)
        if not expanded_idx:
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

        # -------- 调 Token-Mol 双 pocket 先验 --------
        prefix = ["<|beginoftext|>", "<|mask:0|>", "<|mask:0|>"]
        prefix.extend([t for t in state.position if t != "&"])

        try:
            priors = call_prior(prefix, candidates, PRIOR_URL, P1, P2)
            prior_map = {tok: float(p) for tok, p in zip(candidates, priors)}
            print(f"✅ 成功获取双靶点先验，候选token: {candidates}")
        except PriorServiceError as e:
            print(f"❌ 先验服务调用失败: {e}")
            print("❌ 终止运行：必须使用Token-Mol先验服务")
            import sys
            sys.exit(1)
        except Exception as e:
            print(f"❌ 未知错误: {e}")
            import sys
            sys.exit(1)

        # Progressive Widening
        max_children = int(1 + (node.visits + 1) ** 0.5)
        quota = max_children - len(node.childNodes)
        allow = max(1, quota)

        ranked = sorted(candidates, key=lambda t: -prior_map.get(t, 0.0))[:allow]
        for tok in ranked:
            ch = node.add_child(tok, state)
            ch.P = prior_map.get(tok, 1.0 / len(candidates))

        # -------- Simulation（沿用原 SFG-Drug）--------
        new_compounds: List[str] = []
        nodeadded_all: List[str] = []
        for _ in range(SIM_NUM):
            nodeadded = node_to_add(expanded_idx, VAL)
            all_possible = chem_kn_simulation(MODEL, state.position, VAL, nodeadded)
            gen = predict_smile(all_possible, VAL)
            new_compounds.extend(make_input_smile(gen))
            nodeadded_all.extend(nodeadded)

        generated_dict = {}
        node_index, vinadock_score, valid_smile, generated_dict = check_node_type(
            new_compounds, generated_dict,
            sa_threshold=SA_THRESHOLD, rule=RULE5,
            radical=RADICAL_CHECK, hashimoto_filter=HASHIMOTO_FILTER,
            dict_id=dict_id, trial=trial
        )
        dict_id += 1
        valid_compounds.extend(valid_smile)

        # ⭐ 新增：逐个分子写入双靶点分数
        now = time.time()
        for smi, s in zip(valid_smile, vinadock_score):
            # s 可能是单个分数，也可能是 [s1, s2]
            if isinstance(s, (list, tuple)) and len(s) >= 2:
                s1 = float(s[0])
                s2 = float(s[1])
                combined = s1 + s2
            else:
                # 兼容旧格式：只有一个 docking 分数
                s1 = float(s)
                s2 = float("nan")
                combined = s1
            out_f.write(
                f"{smi},{s1:.3f},{s2:.3f},{combined:.3f},{len(state.position)},{now - start_time:.3f}\n"
            )
        out_f.flush()

        # -------- Backprop / Reward 计算 --------
        if not node_index:
            cur = node
            while cur is not None:
                cur.update(-1.0)
                cur = cur.parentNode
        else:
            node_pool: List[Node] = []
            rewards: List[float] = []
            seen_atom = set()

            for i, idx in enumerate(node_index):
                if idx < 0 or idx >= len(nodeadded_all):
                    continue
                atom = nodeadded_all[idx]

                if atom not in seen_atom:
                    seen_atom.add(atom)
                    ch = None
                    for c in node.childNodes:
                        if c.position == atom:
                            ch = c
                            break
                    if ch is None:
                        ch = node.add_child(atom, state)
                        ch.P = 1.0 / max(1, len(node.childNodes))
                    node_pool.append(ch)

                # (1) docking reward：combined = s1 + s2
                s = vinadock_score[i]
                if isinstance(s, (list, tuple)) and len(s) >= 2:
                    s1 = float(s[0])
                    s2 = float(s[1])
                    combined = s1 + s2
                else:
                    combined = float(s[0]) if isinstance(s, (list, tuple)) else float(s)

                # 基准是 2 * BASE_VINA（两个靶点都在 base_vinadock_score 附近）
                # B1(min) 融合版：同时考虑「总和」与「最差靶点」以促使双靶点折中
                base = float(BASE_VINA)

                # sum 项：整体表现
                delta_sum = combined - 2.0 * base
                r_sum = (-(delta_sum) * 0.1) / (1.0 + abs(delta_sum) * 0.1)

                # min 项：短板靶点（越差惩罚越大）
                if isinstance(s, (list, tuple)) and len(s) >= 2:
                    s_min = min(s1, s2)
                    delta_min = s_min - base
                    r_min = (-(delta_min) * 0.2) / (1.0 + abs(delta_min) * 0.2)
                    r_dock = 0.5 * r_sum + 0.5 * r_min
                else:
                    # 兼容单分数（例如只跑单靶点时）
                    r_dock = r_sum
# (2) QED / SA / 多样性
                smi = valid_smile[i] if i < len(valid_smile) else None
                r = r_dock
                if smi:
                    qed_v = safe_qed(smi)
                    sa_v = safe_sa(smi)
                    sa_term = max(0.0, 1.0 - (sa_v - 2.5) / 3.5)

                    fp = morgan_fp(smi)
                    is_novel = (fp not in seen_fps)
                    if is_novel:
                        seen_fps.add(fp)

                    r = (
                        W_DOCK * r_dock
                        + W_QED * qed_v
                        + W_SA * sa_term
                        + (DIVERSITY_BONUS if is_novel else 0.0)
                    )

                rewards.append(float(r))

            # 反向传播到祖先节点
            for ch, r in zip(node_pool, rewards):
                cur = ch
                while cur is not None:
                    cur.update(r)
                    cur = cur.parentNode

        elapsed = (time.time() - start_time) / 60.0
        left = max((deadline - time.time()) / 60.0, 0.0)
        print(f"[progress] {elapsed:.2f}m elapsed, {left:.2f}m left, valid={len(valid_compounds)}")

    out_f.close()
    print("=== PRIOR MCTS finished ===")
    print("total valid compounds:", len(valid_compounds))
    print("results saved to:", OUTPUT_PATH)
    return valid_compounds


# =========================================================
#  主入口 & 先验健康检查
# =========================================================
def check_prior_service():
    """启动前检查先验服务是否可用：能成功调用就算通过，不再禁止均匀分布。"""
    try:
        from prior_client import call_prior, PriorServiceError

        test_prefix = ["<|beginoftext|>", "<|mask:0|>", "<|mask:0|>", "C"]
        test_candidates = ["C", "N", "O", "(", ")", "[", "]"]

        priors = call_prior(test_prefix, test_candidates, PRIOR_URL, P1, P2)

        print("✅ 先验服务检查通过")
        print(f"   检查用先验范围: {min(priors):.4f} - {max(priors):.4f}")
        return True

    except PriorServiceError as e:
        print(f"❌ 先验服务检查失败: {e}")
        return False
    except Exception as e:
        print(f"❌ 先验服务检查异常: {e}")
        return False


if __name__ == "__main__":
    # ZINC 词表
    zinc_path = r"E:\graduate\meeting\mwx\SFG-Drug-code\SFG-Drug\data\250k_rndm_zinc_drugs_clean.smi"
    smile_old = zinc_data_with_bracket_original(zinc_path)
    VAL, _ = zinc_processed_with_bracket(smile_old)
    print("val size:", len(VAL))

    # RNN 片段模型
    MODEL = loaded_model("RNN-model/" + MODEL_NAME)

    # 检查先验服务
    print("🔍 检查Token-Mol先验服务...")
    if not check_prior_service():
        print("❌ 无法连接到Token-Mol先验服务，程序终止")
        import sys
        sys.exit(1)

    print("🚀 开始分子生成（使用Token-Mol双靶点先验）")
    init = ChemicalState()
    _ = mcts(init)
