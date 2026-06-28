from __future__ import annotations

import _bootstrap  # noqa: F401

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, lsq_linear

from cases import DOMAIN_HI
from multistart_base import MultistartConfig
from multistart_config import DEFAULT_MULTISTART_CONFIG
from multistart_configs.ms_v2_stable_qmc import MultistartConfigV2
from multistart_initial_guess import generate_initial_guess
from spectral_neumann import phi_matrix

EPS = 1e-30
RESIDUAL_NORMALIZE_EPS = 1e-30


def _reg_coeff(n_pts: int, strategy_name: str) -> float:
    base = 1e-3 if n_pts <= 8 else 5e-4
    if "regularizacao" in strategy_name:
        base *= 2.0
    return float(base)


@dataclass(frozen=True)
class Stage1Result:
    n_nodes: int
    nodes: np.ndarray
    weights: np.ndarray
    residual_l2: float
    strategy: str
    nfev: int
    success: bool


def stage1_solve_with_strategy(
    n_pts: int,
    vecs: np.ndarray,
    dx: float,
    targets: np.ndarray,
    z_ref: float,
    strategy_name: str,
    *,
    rng: np.random.Generator | None = None,
    case_a: float | None = None,
    case_b: float | None = None,
    hereditary_prior_nodes: np.ndarray | None = None,
    hereditary_prior_weights: np.ndarray | None = None,
    lobatto_use_exact_endpoints: bool = True,
    multistart_seed: int | None = None,
) -> tuple[np.ndarray | None, np.ndarray | None, float, dict]:
    if rng is None:
        rng = np.random.default_rng(0)
    x0, w0 = generate_initial_guess(
        strategy_name,
        n_pts,
        z_ref,
        rng,
        case_a=case_a,
        case_b=case_b,
        hereditary_prior_nodes=hereditary_prior_nodes,
        hereditary_prior_weights=hereditary_prior_weights,
        lobatto_use_exact_endpoints=lobatto_use_exact_endpoints,
        multistart_seed=multistart_seed,
    )
    p0 = np.sqrt(np.maximum(w0, 1e-30))
    initial = np.concatenate([x0, p0])
    targets_arr = np.asarray(targets, dtype=float)
    reg_coeff = _reg_coeff(n_pts, strategy_name)

    def equations(v: np.ndarray) -> np.ndarray:
        xv = v[:n_pts]
        pv = v[n_pts : 2 * n_pts]
        wv = pv**2
        a = phi_matrix(xv, vecs, dx)[: 2 * n_pts, :]
        raw = a @ wv - targets_arr
        spectral_res = raw / np.maximum(np.abs(targets_arr), RESIDUAL_NORMALIZE_EPS)
        residuals = spectral_res.tolist()
        reg = reg_coeff * ((float(np.sum(wv)) - z_ref) / max(z_ref, EPS)) ** 2
        residuals.append(reg)
        return np.asarray(residuals, dtype=float)

    lo = np.concatenate([np.zeros(n_pts), np.zeros(n_pts)])
    hi = np.concatenate([np.full(n_pts, DOMAIN_HI), np.full(n_pts, np.inf)])
    try:
        result = least_squares(equations, initial, bounds=(lo, hi), max_nfev=1600, xtol=1e-12, ftol=1e-12)
    except Exception as exc:
        return None, None, float("inf"), {"success": False, "message": str(exc), "nfev": 0}
    xv = result.x[:n_pts]
    wv = result.x[n_pts : 2 * n_pts] ** 2
    a = phi_matrix(xv, vecs, dx)[: 2 * n_pts, :]
    raw_r = a @ wv - targets_arr
    resn = float(np.linalg.norm(raw_r))
    ok = bool(result.success) or resn < max(1e-9, 5e-4 * max(float(np.linalg.norm(targets_arr)), EPS))
    if not ok:
        return None, None, float("inf"), {"success": False, "nfev": int(getattr(result, "nfev", 0))}
    return xv, wv, resn, {"success": bool(result.success), "nfev": int(getattr(result, "nfev", 0))}


def _fallback(n_pts: int, vecs: np.ndarray, dx: float, targets: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    nodes = np.linspace(0.05 * DOMAIN_HI, 0.95 * DOMAIN_HI, n_pts)
    a = phi_matrix(nodes, vecs, dx)[: 2 * n_pts, :]
    res = lsq_linear(a, np.asarray(targets, dtype=float), bounds=(0.0, np.inf), lsmr_tol="auto")
    w = np.maximum(np.asarray(res.x, dtype=float), 0.0)
    rr = a @ w - np.asarray(targets, dtype=float)
    return nodes, w, float(np.linalg.norm(rr))


def stage1_solve(
    n_pts: int,
    vecs: np.ndarray,
    dx: float,
    targets: np.ndarray,
    z_ref: float,
    *,
    random_uniform_trials: int = 20,
    parallel: bool = True,
    max_workers: int | None = None,
    multistart_config: MultistartConfig | MultistartConfigV2 | None = None,
    hereditary_prior: tuple[np.ndarray, np.ndarray] | None = None,
    case_a: float | None = None,
    case_b: float | None = None,
) -> Stage1Result:
    from pathlib import Path

    from multistart import build_candidate_jobs, evaluate_candidates_parallel, evaluate_candidates_sequential, select_best_candidate

    cfg = DEFAULT_MULTISTART_CONFIG if multistart_config is None else multistart_config
    jobs = build_candidate_jobs(
        cfg, n_pts, random_uniform_trials, hereditary_prior=hereditary_prior, case_a=case_a, case_b=case_b
    )

    def _eval_one(job):
        rng = np.random.default_rng(int(job.seed))
        return stage1_solve_with_strategy(
            n_pts,
            vecs,
            dx,
            targets,
            z_ref,
            job.solve_strategy,
            rng=rng,
            case_a=job.case_a,
            case_b=job.case_b,
            hereditary_prior_nodes=None
            if job.hereditary_prior_nodes is None
            else np.asarray(job.hereditary_prior_nodes, dtype=float),
            hereditary_prior_weights=None
            if job.hereditary_prior_weights is None
            else np.asarray(job.hereditary_prior_weights, dtype=float),
            lobatto_use_exact_endpoints=job.lobatto_use_exact_endpoints,
            multistart_seed=int(job.seed),
        )

    if parallel and len(jobs) > 1:
        outcomes = evaluate_candidates_parallel(
            jobs,
            project_dir=str(Path(__file__).resolve().parent),
            n_pts=n_pts,
            dx=dx,
            z_ref=z_ref,
            vecs=vecs,
            targets=targets,
            max_workers=max_workers,
        )
    else:
        outcomes = evaluate_candidates_sequential(jobs, _eval_one)

    best = select_best_candidate(outcomes)
    if best is not None:
        return Stage1Result(
            n_nodes=n_pts,
            nodes=np.asarray(best.nodes, dtype=float),
            weights=np.asarray(best.weights, dtype=float),
            residual_l2=float(best.residual_l2),
            strategy=best.label,
            nfev=best.nfev,
            success=True,
        )
    best_nodes, best_weights, best_res = _fallback(n_pts, vecs, dx, targets)
    return Stage1Result(
        n_nodes=n_pts,
        nodes=np.asarray(best_nodes, dtype=float),
        weights=np.asarray(best_weights, dtype=float),
        residual_l2=float(best_res),
        strategy="fallback_lsq_linear",
        nfev=0,
        success=True,
    )
