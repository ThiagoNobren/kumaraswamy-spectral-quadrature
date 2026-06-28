from __future__ import annotations

import numpy as np

from analytic_qoi import AnalyticQoI
from cases import GRID_SIZE, KumaraswamyCase
from kumaraswamy import cdf

EPS = 1e-30


def _rel(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(b)), EPS)


def moments_from_rule(nodes: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    th = np.asarray(nodes, dtype=float).ravel()
    w = np.asarray(weights, dtype=float).ravel()
    z_hat = float(np.sum(w))
    if z_hat <= 0.0:
        raise ValueError("Soma dos pesos deve ser positiva.")
    mu = float(np.dot(w, th) / z_hat)
    m2 = float(np.dot(w, th * th) / z_hat)
    var = max(float(m2 - mu * mu), 0.0)
    return mu, var


def cdf_from_rule(a: float, nodes: np.ndarray, weights: np.ndarray) -> float:
    th = np.asarray(nodes, dtype=float).ravel()
    w = np.asarray(weights, dtype=float).ravel()
    order = np.argsort(th)
    th = th[order]
    w = w[order]
    z_hat = float(np.sum(w))
    if z_hat <= 0.0:
        raise ValueError("Soma dos pesos deve ser positiva.")
    if a <= float(th[0]):
        return 0.0
    if a >= float(th[-1]):
        return 1.0
    return float(np.sum(w[th <= float(a)])) / z_hat


def quantile_from_rule(rho: float, nodes: np.ndarray, weights: np.ndarray) -> float:
    if not (0.0 < rho < 1.0):
        raise ValueError("rho deve pertencer a (0,1).")
    th = np.asarray(nodes, dtype=float).ravel()
    w = np.asarray(weights, dtype=float).ravel()
    order = np.argsort(th)
    th = th[order]
    w = w[order]
    cw = np.cumsum(w)
    z_hat = float(cw[-1])
    if z_hat <= 0.0:
        raise ValueError("Soma dos pesos deve ser positiva.")
    target = float(rho) * z_hat
    idx = int(np.searchsorted(cw, target, side="left"))
    idx = min(max(idx, 0), len(th) - 1)
    return float(th[idx])


def cdf_sup_and_w1(
    case: KumaraswamyCase,
    nodes: np.ndarray,
    weights: np.ndarray,
    *,
    grid_size: int = GRID_SIZE,
) -> tuple[float, float]:
    grid = np.linspace(0.0, 1.0, int(grid_size), dtype=float)
    cdf_ref = np.array([float(cdf(float(t), case)) for t in grid], dtype=float)
    cdf_hat = np.array([cdf_from_rule(float(t), nodes, weights) for t in grid], dtype=float)
    err = np.abs(cdf_hat - cdf_ref)
    cdf_sup = float(np.max(err))
    w1 = float(np.trapezoid(err, grid))
    return cdf_sup, w1


def compute_qoi_metrics(
    case: KumaraswamyCase,
    nodes: np.ndarray,
    weights: np.ndarray,
    ref: AnalyticQoI,
    *,
    quantile_levels: list[float],
    grid_size: int = GRID_SIZE,
) -> dict[str, float]:
    mu_hat, var_hat = moments_from_rule(nodes, weights)
    out: dict[str, float] = {
        "rel_mean": _rel(mu_hat, ref.mean),
        "rel_var": _rel(var_hat, ref.variance),
    }
    rel_qs: list[float] = []
    for p in quantile_levels:
        key = f"q{int(round(100 * p))}"
        q_ref = ref.quantiles[key]
        q_hat = quantile_from_rule(p, nodes, weights)
        out[f"abs_{key}"] = abs(q_hat - q_ref)
        rel_q = _rel(q_hat, q_ref)
        out[f"rel_{key}"] = rel_q
        rel_qs.append(rel_q)
    out["rel_quant_mean"] = float(np.mean(rel_qs))
    cdf_sup, w1 = cdf_sup_and_w1(case, nodes, weights, grid_size=grid_size)
    out["CDF_sup"] = cdf_sup
    out["W1"] = w1
    return out
