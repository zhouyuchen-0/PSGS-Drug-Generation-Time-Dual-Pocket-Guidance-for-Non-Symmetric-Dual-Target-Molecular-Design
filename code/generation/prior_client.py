#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
HTTP client for the Token-Mol pocket-prior service.

The client queries the local prior server for target-specific token priors and returns a fused prior distribution for candidate next tokens.
"""

# prior_client.py
import requests
import sys
# prior_client.py(   SFG-Drug    )
DEFAULT_PRIOR_URL = "http://127.0.0.1:26974/prior"

DEFAULT_PROTEIN_PATH_1 = "receptors/3fap.pkl"
DEFAULT_PROTEIN_PATH_2 = "receptors/7pqv.pkl"



class PriorServiceError(Exception):
    """      """
    pass


def call_prior(prefix_tokens, candidates, prior_url, protein_path_1, protein_path_2):
    """
               .      ,      ,       ,
       "    "     --                    .
    """
    try:
        # 1)       
        health_url = prior_url.replace('/prior', '/health')
        try:
            health_response = requests.get(health_url, timeout=10)
            if health_response.status_code != 200:
                raise PriorServiceError(f"          : {health_response.status_code}")

            health_data = health_response.json()
            if not health_data.get('models_loaded', False) or not health_data.get('pockets_loaded', False):
                raise PriorServiceError("              ")
        except requests.exceptions.RequestException:
            raise PriorServiceError("               ")

        # 2)           
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
            raise PriorServiceError(f"3fap      : {response_1.status_code}")

        data_1 = response_1.json()
        status_1 = data_1.get("status", None)
        msg_1 = data_1.get("message", "")
        if status_1 is not None and status_1 != "ok":
            raise PriorServiceError(f"3fap      : {msg_1}")

        prior_1 = data_1.get("prior", None)
        if prior_1 is None or len(prior_1) != len(candidates):
            raise PriorServiceError("3fap         ")

        # 3)           
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
            raise PriorServiceError(f"7pqv      : {response_2.status_code}")

        data_2 = response_2.json()
        status_2 = data_2.get("status", None)
        msg_2 = data_2.get("message", "")
        if status_2 is not None and status_2 != "ok":
            raise PriorServiceError(f"7pqv      : {msg_2}")

        prior_2 = data_2.get("prior", None)
        if prior_2 is None or len(prior_2) != len(candidates):
            raise PriorServiceError("7pqv         ")

        #      "  "   ,     
        def is_uniform(arr):
            if not arr:
                return True
            vals = [round(float(x), 6) for x in arr]
            return len(set(vals)) == 1

        if is_uniform(prior_1) or is_uniform(prior_2):
            print("[WARN] [prior_client]                  ,"
                  "                top-k        ,      .")

        # 4)         
        combined_prior = []
        for i in range(len(candidates)):
            p = 0.5 * float(prior_1[i]) + 0.5 * float(prior_2[i])
            combined_prior.append(p)

        total = sum(combined_prior)
        if total <= 0:
            raise PriorServiceError("         0")

        combined_prior = [p / total for p in combined_prior]

        #        status/message/  ,               (     /  )
        print(
            f"[OK]          | "
            f"3fap status={status_1}, msg={str(msg_1)[:120]} | range={min(prior_1):.4f}-{max(prior_1):.4f} ; "
            f"7pqv status={status_2}, msg={str(msg_2)[:120]} | range={min(prior_2):.4f}-{max(prior_2):.4f}"
        )
        return combined_prior

    except requests.exceptions.RequestException as e:
        raise PriorServiceError(f"      : {e}")
    except PriorServiceError:
        #             ,     
        raise
    except Exception as e:
        raise PriorServiceError(f"      : {e}")


#      /       
__all__ = ['PriorServiceError', 'call_prior']