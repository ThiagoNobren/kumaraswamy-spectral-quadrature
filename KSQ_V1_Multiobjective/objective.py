from __future__ import annotations

import _bootstrap  # noqa: F401

from dataclasses import dataclass

import numpy as np

from analytic_qoi import analytic_qoi
from cases import QUANTILE_LEVELS, KumaraswamyCase
from qoi_metrics import compute_qoi_metrics

EPS = 1e-30


def compute_objective_terms(
    case: KumaraswamyCase,
    nodes: np.ndarray,
    weights: np.ndarray,
    *,
    E_r: float,
) -> dict[str, float]:
    ref = analytic_qoi(case)
    metrics = compute_qoi_metrics(case, nodes, weights, ref, quantile_levels=QUANTILE_LEVELS)
    return {
        "E_r": float(E_r),
        "E_mean": float(metrics["rel_mean"]),
        "E_var": float(metrics["rel_var"]),
        "E_quant": float(metrics["rel_quant_mean"]),
    }


def compute_multiobjective_F(
    terms: dict[str, float],
    *,
    lambda_r: float,
    lambda_m: float,
    lambda_v: float,
    lambda_q: float,
) -> float:
    return (
        float(lambda_r) * terms["E_r"]
        + float(lambda_m) * terms["E_mean"]
        + float(lambda_v) * terms["E_var"]
        + float(lambda_q) * terms["E_quant"]
    )


def weighted_contributions(
    terms: dict[str, float],
    *,
    lambda_r: float,
    lambda_m: float,
    lambda_v: float,
    lambda_q: float,
) -> dict[str, float]:
    return {
        "contrib_r": float(lambda_r * terms["E_r"]),
        "contrib_m": float(lambda_m * terms["E_mean"]),
        "contrib_v": float(lambda_v * terms["E_var"]),
        "contrib_q": float(lambda_q * terms["E_quant"]),
    }


def eta_from_weights(weights: np.ndarray, z_ref: float) -> np.ndarray:
    w = np.maximum(np.asarray(weights, dtype=float).ravel(), EPS)
    z = max(float(z_ref), EPS)
    eta = np.log(w / z)
    eta -= float(np.mean(eta))
    return eta


def weights_from_eta(eta: np.ndarray, z_ref: float) -> np.ndarray:
    e = np.asarray(eta, dtype=float).ravel()
    ex = np.exp(e - np.max(e))
    w = ex / max(float(np.sum(ex)), EPS)
    return w * float(z_ref)
