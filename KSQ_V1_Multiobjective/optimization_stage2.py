from __future__ import annotations

import _bootstrap  # noqa: F401

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from cases import DOMAIN_HI, KumaraswamyCase, LAMBDA_M, LAMBDA_Q, LAMBDA_R, LAMBDA_V
from objective import (
    compute_multiobjective_F,
    compute_objective_terms,
    eta_from_weights,
    weighted_contributions,
    weights_from_eta,
)
from spectral_neumann import spectral_moment_residuals

EPS = 1e-30


@dataclass(frozen=True)
class Stage2Result:
    nodes: np.ndarray
    weights: np.ndarray
    score_F: float
    relative_spectral_error: float
    terms: dict[str, float]
    contributions: dict[str, float]
    converged: bool
    stop_message: str
    nfev: int
    nit: int


def stage2_multiobjective_refine(
    case: KumaraswamyCase,
    theta_init: np.ndarray,
    weights_init: np.ndarray,
    *,
    targets: np.ndarray,
    vecs: np.ndarray,
    dx: float,
    z_ref: float,
    lambda_r: float = LAMBDA_R,
    lambda_m: float = LAMBDA_M,
    lambda_v: float = LAMBDA_V,
    lambda_q: float = LAMBDA_Q,
    maxiter: int = 500,
    ftol: float = 1e-8,
) -> Stage2Result:
    th0 = np.clip(np.asarray(theta_init, dtype=float).ravel(), 0.0, DOMAIN_HI)
    eta0 = eta_from_weights(weights_init, z_ref)
    x0 = np.concatenate([th0, eta0])
    m_norm = max(float(np.linalg.norm(np.asarray(targets, dtype=float))), EPS)

    def objective(x: np.ndarray) -> float:
        th = np.clip(x[: len(th0)], 0.0, DOMAIN_HI)
        eta = x[len(th0) :]
        w = weights_from_eta(eta, z_ref)
        raw_r = spectral_moment_residuals(th, w, targets, vecs, dx)
        e_r = float(np.linalg.norm(raw_r)) / m_norm
        terms = compute_objective_terms(case, th, w, E_r=e_r)
        return compute_multiobjective_F(
            terms,
            lambda_r=lambda_r,
            lambda_m=lambda_m,
            lambda_v=lambda_v,
            lambda_q=lambda_q,
        )

    bounds = [(0.0, DOMAIN_HI)] * len(th0) + [(None, None)] * len(th0)
    result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds, options={"maxiter": maxiter, "ftol": ftol})
    x_opt = np.asarray(result.x, dtype=float)
    th = np.clip(x_opt[: len(th0)], 0.0, DOMAIN_HI)
    w = weights_from_eta(x_opt[len(th0) :], z_ref)
    raw_r = spectral_moment_residuals(th, w, targets, vecs, dx)
    e_r = float(np.linalg.norm(raw_r)) / m_norm
    terms = compute_objective_terms(case, th, w, E_r=e_r)
    contributions = weighted_contributions(
        terms,
        lambda_r=lambda_r,
        lambda_m=lambda_m,
        lambda_v=lambda_v,
        lambda_q=lambda_q,
    )
    return Stage2Result(
        nodes=th,
        weights=w,
        score_F=float(sum(contributions.values())),
        relative_spectral_error=float(e_r),
        terms=terms,
        contributions=contributions,
        converged=bool(result.success),
        stop_message=str(getattr(result, "message", "") or ""),
        nfev=int(getattr(result, "nfev", 0)),
        nit=int(getattr(result, "nit", 0)),
    )
