# prior_client.py
import requests
import sys
# prior_client.py（本地 SFG-Drug 项目里）
DEFAULT_PRIOR_URL = "http://127.0.0.1:26974/prior"

DEFAULT_PROTEIN_PATH_1 = "/root/Token-Mol-main/example/3fap.pkl"
DEFAULT_PROTEIN_PATH_2 = "/root/Token-Mol-main/exampl/7pqv_mek1.pkl"



class PriorServiceError(Exception):
    """先验服务异常"""
    pass


def call_prior(prefix_tokens, candidates, prior_url, protein_path_1, protein_path_2):
    """
    请求两个靶点蛋白的先验。只要服务健康、返回格式正确，就认为先验有效，
    不再把“均匀分布”当作错误 —— 均匀只表示这一步模型对候选没有明显偏好。
    """
    try:
        # 1) 先做健康检查
        health_url = prior_url.replace('/prior', '/health')
        try:
            health_response = requests.get(health_url, timeout=10)
            if health_response.status_code != 200:
                raise PriorServiceError(f"先验服务健康检查失败: {health_response.status_code}")

            health_data = health_response.json()
            if not health_data.get('models_loaded', False) or not health_data.get('pockets_loaded', False):
                raise PriorServiceError("先验服务模型或口袋未正确加载")
        except requests.exceptions.RequestException:
            raise PriorServiceError("无法连接到先验服务健康检查端点")

        # 2) 请求第一个蛋白的先验
        response_1 = requests.post(
            prior_url,
            json={
                "prefix_tokens": prefix_tokens,
                "candidates": candidates,
                "protein_path": protein_path_1
            },
            timeout=30
        )
        if response_1.status_code != 200:
            raise PriorServiceError(f"3fap先验请求失败: {response_1.status_code}")

        data_1 = response_1.json()
        status_1 = data_1.get("status", None)
        msg_1 = data_1.get("message", "")
        if status_1 is not None and status_1 != "ok":
            raise PriorServiceError(f"3fap先验计算失败: {msg_1}")

        prior_1 = data_1.get("prior", None)
        if prior_1 is None or len(prior_1) != len(candidates):
            raise PriorServiceError("3fap先验返回格式不正确")

        # 3) 请求第二个蛋白的先验
        response_2 = requests.post(
            prior_url,
            json={
                "prefix_tokens": prefix_tokens,
                "candidates": candidates,
                "protein_path": protein_path_2
            },
            timeout=30
        )
        if response_2.status_code != 200:
            raise PriorServiceError(f"7pqv先验请求失败: {response_2.status_code}")

        data_2 = response_2.json()
        status_2 = data_2.get("status", None)
        msg_2 = data_2.get("message", "")
        if status_2 is not None and status_2 != "ok":
            raise PriorServiceError(f"7pqv先验计算失败: {msg_2}")

        prior_2 = data_2.get("prior", None)
        if prior_2 is None or len(prior_2) != len(candidates):
            raise PriorServiceError("7pqv先验返回格式不正确")

        # 这里不再把“均匀”当错误，只做个提示
        def is_uniform(arr):
            if not arr:
                return True
            vals = [round(float(x), 6) for x in arr]
            return len(set(vals)) == 1

        if is_uniform(prior_1) or is_uniform(prior_2):
            print("⚠️ [prior_client] 当前这一步某个靶点先验接近均匀分布，"
                  "可能是该局部对候选无明显偏好或 top-k 与候选重合较少，但仍继续使用。")

        # 4) 加权合并并归一化
        combined_prior = []
        for i in range(len(candidates)):
            p = 0.5 * float(prior_1[i]) + 0.5 * float(prior_2[i])
            combined_prior.append(p)

        total = sum(combined_prior)
        if total <= 0:
            raise PriorServiceError("合并先验概率总和为0")

        combined_prior = [p / total for p in combined_prior]

        # 额外打印一次 status/message/范围，便于排查是否真的使用了模型先验（而不是回退/异常）
        print(
            f"✅ 先验服务调用成功 | "
            f"3fap status={status_1}, msg={str(msg_1)[:120]} | range={min(prior_1):.4f}-{max(prior_1):.4f} ; "
            f"7pqv status={status_2}, msg={str(msg_2)[:120]} | range={min(prior_2):.4f}-{max(prior_2):.4f}"
        )
        return combined_prior

    except requests.exceptions.RequestException as e:
        raise PriorServiceError(f"网络连接失败: {e}")
    except PriorServiceError:
        # 已经是我们明确抛出的错误，直接往外传
        raise
    except Exception as e:
        raise PriorServiceError(f"先验服务异常: {e}")


# 确保这些类/函数可以被导入
__all__ = ['PriorServiceError', 'call_prior']