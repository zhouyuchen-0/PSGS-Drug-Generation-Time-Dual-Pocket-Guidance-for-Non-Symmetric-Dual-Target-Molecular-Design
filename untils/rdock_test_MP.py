# -*- coding: utf-8 -*-
"""
utils/rdock_test_MP.py
缓存增强版（保持 vinadock_score(compound) 接口不变）

功能：
1) 双靶点 Vina 打分（3FAP + 7PQV）
2) 支持缓存（按 canonical SMILES 缓存，避免重复 docking）
3) 支持 debug 日志（SFG_DEBUG_DIAG=1）
4) 路径可通过环境变量覆盖
"""

import os
import re
import time
import json
import shutil
import tempfile
import traceback
import subprocess
from typing import Tuple, List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem

# -----------------------------
# 环境配置 / 默认路径
# -----------------------------
FAIL_SCORE = 10 ** 10
BEST_DOCKING_ID = 0

def _env_bool(name: str, default: bool = False) -> bool:
    v = os.environ.get(name, None)
    if v is None:
        return default
    return str(v).strip().lower() in ("1", "true", "yes", "y", "on")

DEBUG_DIAG = _env_bool("SFG_DEBUG_DIAG", False)
DEBUG_TRACEBACK = _env_bool("SFG_DEBUG_TRACEBACK", False)

VINA_BIN = os.environ.get("SFG_VINA_BIN", shutil.which("vina") or "/usr/bin/vina")
OBABEL_BIN = os.environ.get("SFG_OBABEL_BIN", shutil.which("obabel") or "/usr/bin/obabel")

CFG_3FAP = os.environ.get("SFG_CFG_3FAP", "/mnt/SFG-Drug/pro/3fapconfig.txt")
CFG_7PQV = os.environ.get("SFG_CFG_7PQV", "/mnt/SFG-Drug/pro/7pqvconfig.txt")

# 持久缓存文件（可改）
CACHE_PATH = os.environ.get("SFG_RDOCK_CACHE_PATH", "/tmp/sfgdock_vina_cache.jsonl")

# 内存缓存（当前进程）
_MEM_CACHE = {}


def _dbg(*args):
    if DEBUG_DIAG:
        print("[rdock_debug]", *args)


def _safe_print_tb():
    if DEBUG_TRACEBACK:
        traceback.print_exc()


# -----------------------------
# 缓存工具
# -----------------------------
def _canonical_smiles(smi: str) -> Optional[str]:
    try:
        m = Chem.MolFromSmiles(str(smi))
        if m is None:
            return None
        return Chem.MolToSmiles(m, canonical=True)
    except Exception:
        return None


def _cache_key_from_smiles(smi: str) -> Optional[str]:
    cs = _canonical_smiles(smi)
    return cs


def _load_cache_from_disk():
    global _MEM_CACHE
    if not os.path.exists(CACHE_PATH):
        return
    try:
        with open(CACHE_PATH, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                    k = obj.get("k")
                    v = obj.get("v")
                    if k and isinstance(v, list) and len(v) >= 3:
                        _MEM_CACHE[k] = v
                except Exception:
                    continue
        _dbg("cache loaded", "size=", len(_MEM_CACHE), "path=", CACHE_PATH)
    except Exception:
        _safe_print_tb()


def _append_cache_to_disk(k: str, v: List[float]):
    try:
        d = os.path.dirname(CACHE_PATH)
        if d:
            os.makedirs(d, exist_ok=True)
        with open(CACHE_PATH, "a", encoding="utf-8") as f:
            f.write(json.dumps({"k": k, "v": v}, ensure_ascii=False) + "\n")
    except Exception:
        _safe_print_tb()


_load_cache_from_disk()


# -----------------------------
# 兼容旧文件中的辅助函数（可保留）
# -----------------------------
def savefile(source_file, destination_dir):
    """兼容旧接口：保留但默认不使用"""
    try:
        os.makedirs(destination_dir, exist_ok=True)
        timestamp = int(time.time())
        filename, _ = os.path.splitext(os.path.basename(source_file))
        new_filename = f"{filename}_{timestamp}.pdbqt"
        destination_file = os.path.join(destination_dir, new_filename)
        shutil.copy2(source_file, destination_file)
    except Exception:
        _safe_print_tb()


def savesmile(source_file, destination_dir):
    """兼容旧接口：保留但默认不使用"""
    try:
        os.makedirs(destination_dir, exist_ok=True)
        timestamp = int(time.time())
        filename, _ = os.path.splitext(os.path.basename(source_file))
        new_filename = f"{filename}_{timestamp}.pdb"
        destination_file = os.path.join(destination_dir, new_filename)
        shutil.copy2(source_file, destination_file)
    except Exception:
        _safe_print_tb()


def docking_calculation(cmd):
    """兼容旧接口：保留"""
    try:
        subprocess.call(cmd, shell=True)
    except Exception:
        _safe_print_tb()


# -----------------------------
# Vina / OpenBabel 执行与解析
# -----------------------------
def _run_cmd(cmd: List[str], timeout: int = 600) -> Tuple[int, str, str]:
    _dbg("RUN:", " ".join(cmd))
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
    return p.returncode, p.stdout or "", p.stderr or ""


def _assert_file_exists(path: str, tag: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{tag} not found: {path}")


def _parse_vina_score_from_pdbqt(pdbqt_path: str) -> Optional[float]:
    # REMARK VINA RESULT:   -7.2      0.000      0.000
    try:
        if not os.path.exists(pdbqt_path):
            return None
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.strip().split()
                    # 形如：["REMARK","VINA","RESULT:","-7.2","0.000","0.000"]
                    if len(parts) >= 4:
                        try:
                            return float(parts[3])
                        except Exception:
                            pass
        return None
    except Exception:
        _safe_print_tb()
        return None


def _parse_vina_score_from_log(log_path: str) -> Optional[float]:
    # 解析表格第一行结果（最优 mode）
    try:
        if not os.path.exists(log_path):
            return None
        with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        # 找到表头后第一条模式行
        # 例如：   1         -2.4      0.000      0.000
        for line in lines:
            s = line.strip()
            if not s:
                continue
            if s.startswith("REMARK VINA RESULT:"):
                parts = s.split()
                if len(parts) >= 4:
                    try:
                        return float(parts[3])
                    except Exception:
                        pass
            if re.match(r"^\d+\s+-?\d+(\.\d+)?\s+", s):
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        return float(parts[1])
                    except Exception:
                        pass
        return None
    except Exception:
        _safe_print_tb()
        return None


def extract_scores(filename):
    """兼容旧接口：返回 PDBQT 中所有 Vina score 列表"""
    scores = []
    try:
        if not os.path.exists(filename):
            return scores
        with open(filename, "r", encoding="utf-8", errors="ignore") as file:
            for line in file:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            scores.append(float(parts[3]))
                        except Exception:
                            pass
    except Exception:
        _safe_print_tb()
    return scores


def _prepare_ligand_pdbqt(input_smiles: str, workdir: str) -> str:
    mol = Chem.MolFromSmiles(input_smiles)
    if mol is None:
        raise ValueError("MolFromSmiles failed")

    mol = Chem.AddHs(mol)
    # 嵌入 3D（固定随机种子，提升复现性）
    emb = AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
    if emb != 0:
        # 再尝试一次 ETKDG
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        emb = AllChem.EmbedMolecule(mol, params)
        if emb != 0:
            raise RuntimeError("EmbedMolecule failed")

    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=100)
    except Exception:
        # UFF 优化失败不致命
        pass

    pdb_path = os.path.join(workdir, "ligand.pdb")
    pdbqt_path = os.path.join(workdir, "ligand.pdbqt")
    Chem.MolToPDBFile(mol, pdb_path)

    # 使用更稳妥的显式输入/输出格式
    code, out, err = _run_cmd(
        [OBABEL_BIN, "-ipdb", pdb_path, "-opdbqt", "-O", pdbqt_path, "-h"],
        timeout=120
    )
    if code != 0:
        raise RuntimeError("obabel failed: " + (err[:500] if err else out[:500]))
    _assert_file_exists(pdbqt_path, "ligand.pdbqt")
    return pdbqt_path


def _run_vina(cfg_path: str, ligand_pdbqt: str, out_pdbqt: str, log_path: str) -> float:
    _assert_file_exists(cfg_path, "vina config")
    _assert_file_exists(ligand_pdbqt, "ligand pdbqt")

    cmd = [
        VINA_BIN,
        "--config", cfg_path,
        "--ligand", ligand_pdbqt,
        "--out", out_pdbqt,
        "--log", log_path,
    ]
    code, out, err = _run_cmd(cmd, timeout=1200)

    if DEBUG_DIAG and out:
        _dbg("vina stdout head:", out[:500])

    if code != 0:
        raise RuntimeError("vina failed: " + (err[:1000] if err else out[:1000]))

    score = _parse_vina_score_from_pdbqt(out_pdbqt)
    if score is None:
        score = _parse_vina_score_from_log(log_path)
    if score is None:
        raise RuntimeError(f"cannot parse vina score: out={out_pdbqt}, log={log_path}")
    return float(score)


def run_vina1(input_file, output_file):
    """兼容旧接口：仅运行 3FAP（返回 stdout, stderr 风格不再严格保证）"""
    log_path = output_file + ".log"
    try:
        score = _run_vina(CFG_3FAP, input_file, output_file, log_path)
        msg = f"vina1 score={score}"
        return msg.encode("utf-8"), b""
    except Exception as e:
        _safe_print_tb()
        return b"", str(e).encode("utf-8")


def run_vina2(input_file, output_file):
    """兼容旧接口：仅运行 7PQV"""
    log_path = output_file + ".log"
    try:
        score = _run_vina(CFG_7PQV, input_file, output_file, log_path)
        msg = f"vina2 score={score}"
        return msg.encode("utf-8"), b""
    except Exception as e:
        _safe_print_tb()
        return b"", str(e).encode("utf-8")


# -----------------------------
# 主接口（保持不变）
# -----------------------------
def vinadock_score(compound):
    """
    输入：SMILES（或可转为字符串）
    输出：[min_score1, min_score2, best_docking_id]
    与旧版接口保持一致
    """
    input_smiles = str(compound).strip()
    k = _cache_key_from_smiles(input_smiles)

    # 命中缓存
    if k is not None and k in _MEM_CACHE:
        v = _MEM_CACHE[k]
        _dbg("CACHE HIT", k, "->", v)
        return v

    # 基本环境检查
    try:
        _dbg("VINA_BIN =", VINA_BIN)
        _dbg("OBABEL_BIN =", OBABEL_BIN)
        _dbg("CFG_3FAP =", CFG_3FAP)
        _dbg("CFG_7PQV =", CFG_7PQV)

        if not shutil.which(os.path.basename(VINA_BIN)) and not os.path.exists(VINA_BIN):
            raise FileNotFoundError(f"vina not found: {VINA_BIN}")
        if not shutil.which(os.path.basename(OBABEL_BIN)) and not os.path.exists(OBABEL_BIN):
            raise FileNotFoundError(f"obabel not found: {OBABEL_BIN}")
        _assert_file_exists(CFG_3FAP, "CFG_3FAP")
        _assert_file_exists(CFG_7PQV, "CFG_7PQV")

        t0 = time.time()
        with tempfile.TemporaryDirectory(prefix="sfgdock_") as workdir:
            _dbg("workdir =", workdir)

            ligand_pdbqt = _prepare_ligand_pdbqt(input_smiles, workdir)

            out1 = os.path.join(workdir, "smile_out1.pdbqt")
            log1 = os.path.join(workdir, "vina_3fap.log")
            s1 = _run_vina(CFG_3FAP, ligand_pdbqt, out1, log1)

            out2 = os.path.join(workdir, "smile_out2.pdbqt")
            log2 = os.path.join(workdir, "vina_7pqv.log")
            s2 = _run_vina(CFG_7PQV, ligand_pdbqt, out2, log2)

            t1 = time.time()
            _dbg("docking time_used =", round(t1 - t0, 3), "s")

        ret = [float(s1), float(s2), BEST_DOCKING_ID]

        # 写缓存
        if k is not None:
            _MEM_CACHE[k] = ret
            _append_cache_to_disk(k, ret)

        return ret

    except Exception:
        print("error")
        _safe_print_tb()
        ret = [FAIL_SCORE, FAIL_SCORE, BEST_DOCKING_ID]
        # 失败结果不缓存（避免污染）
        return ret