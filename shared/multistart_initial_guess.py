"""Geracao de chutes iniciais para multi-start KSQ V1 (V1 + V2)."""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import legder, leggauss, legroots
from scipy.special import roots_jacobi
from scipy.stats import qmc

from cases import DOMAIN_HI

EPS = 1e-30
BOUNDARY_EPS = 1e-12

HEREDITARY_SIGMA_NODES = 0.01
HEREDITARY_SIGMA_WEIGHTS = 0.01
HEREDITARY_NEW_WEIGHT_FRAC = 0.01


def _normalize_weights(w: np.ndarray, z_ref: float) -> np.ndarray:
    w = np.maximum(np.asarray(w, dtype=float), 1e-16)
    w *= z_ref / max(float(np.sum(w)), EPS)
    return w


def _legendre_lobatto_nodes_01(n: int, *, use_exact_endpoints: bool) -> np.ndarray:
    if n == 1:
        return np.array([0.5 * DOMAIN_HI])
    if n == 2:
        if use_exact_endpoints:
            return np.array([0.0, DOMAIN_HI])
        return np.array([BOUNDARY_EPS, DOMAIN_HI - BOUNDARY_EPS])
    coeffs = [0.0] * (n - 1) + [1.0]
    interior = legroots(legder(coeffs))
    nodes = np.concatenate([[-1.0], interior, [1.0]])
    x = 0.5 * (nodes + 1.0) * DOMAIN_HI
    if use_exact_endpoints:
        x[0] = 0.0
        x[-1] = DOMAIN_HI
    else:
        x = np.clip(x, BOUNDARY_EPS, DOMAIN_HI - BOUNDARY_EPS)
    return np.sort(x)


def _qmc_unit_samples(strategy: str, n: int, seed: int) -> np.ndarray:
    if strategy == "halton_sequence":
        sampler = qmc.Halton(d=1, scramble=True, seed=int(seed))
    elif strategy == "latin_hypercube":
        sampler = qmc.LatinHypercube(d=1, seed=int(seed))
    else:
        raise ValueError(f"estrategia QMC desconhecida: {strategy}")
    return np.asarray(sampler.random(n), dtype=float).ravel()


def generate_initial_guess(
    strategy_name: str,
    n: int,
    z_ref: float,
    rng: np.random.Generator,
    *,
    case_a: float | None = None,
    case_b: float | None = None,
    hereditary_prior_nodes: np.ndarray | None = None,
    hereditary_prior_weights: np.ndarray | None = None,
    lobatto_use_exact_endpoints: bool = True,
    multistart_seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Chute inicial (nos, pesos). Estrategias V1 e V2."""
    if strategy_name == "hereditary_N_minus_1":
        if hereditary_prior_nodes is None or hereditary_prior_weights is None:
            return generate_initial_guess(
                "gauss_legendre", n, z_ref, rng, multistart_seed=multistart_seed
            )
        return _hereditary_guess(
            hereditary_prior_nodes, hereditary_prior_weights, n, z_ref, rng
        )

    if strategy_name in ("halton_sequence", "latin_hypercube"):
        qmc_seed = int(multistart_seed) if multistart_seed is not None else int(rng.integers(0, 2**31 - 1))
        u = _qmc_unit_samples(strategy_name, n, qmc_seed)
        x = np.sort(np.clip(u * DOMAIN_HI, BOUNDARY_EPS, DOMAIN_HI - BOUNDARY_EPS))
        w = np.ones(n) * (z_ref / n) * (1.0 + 0.05 * rng.standard_normal(n))
        return np.sort(x), _normalize_weights(w, z_ref)

    if strategy_name == "chebyshev":
        x = np.cos((2 * np.arange(1, n + 1) - 1) * np.pi / (2 * n)) * DOMAIN_HI / 2 + DOMAIN_HI / 2
        if n % 2 == 1:
            x[n // 2] = DOMAIN_HI / 2
        w = np.ones(n) * (z_ref / n) * (1.0 + 0.05 * rng.standard_normal(n))
    elif strategy_name == "gauss_legendre":
        gl_nodes, gl_weights = leggauss(n)
        x = 0.5 * (gl_nodes + 1.0) * DOMAIN_HI
        w = 0.5 * DOMAIN_HI * gl_weights
        w *= z_ref / max(float(np.sum(w)), EPS)
        w *= 1.0 + 0.05 * rng.standard_normal(n)
    elif strategy_name == "gauss_jacobi_approx":
        alpha = max(float(case_a or 2.0), 0.51)
        beta = max(float(case_b or 2.0), 0.51)
        nodes_j, weights_j = roots_jacobi(n, alpha - 1.0, beta - 1.0)
        x = 0.5 * (nodes_j + 1.0) * DOMAIN_HI
        w = 0.5 * DOMAIN_HI * weights_j
        w *= z_ref / max(float(np.sum(w)), EPS)
        w *= 1.0 + 0.04 * rng.standard_normal(n)
    elif strategy_name == "gauss_lobatto":
        x = _legendre_lobatto_nodes_01(n, use_exact_endpoints=lobatto_use_exact_endpoints)
        w = np.ones(n) * (z_ref / n) * (1.0 + 0.04 * rng.standard_normal(n))
    elif strategy_name == "uniforme" or strategy_name.startswith("uniforme_extra"):
        x = np.linspace(0.1 * DOMAIN_HI, 0.9 * DOMAIN_HI, n) + 0.01 * DOMAIN_HI * rng.standard_normal(n)
        x = np.sort(np.clip(x, 0.05 * DOMAIN_HI, 0.95 * DOMAIN_HI))
        w = np.ones(n) * (z_ref / n) * (1.0 + 0.1 * rng.standard_normal(n))
    elif strategy_name == "regularizacao_minima":
        x = np.linspace(0.1 * DOMAIN_HI, 0.9 * DOMAIN_HI, n)
        w = np.ones(n) * (z_ref / n)
    else:
        cheb = np.cos((2 * np.arange(1, n + 1) - 1) * np.pi / (2 * n)) * DOMAIN_HI / 2 + DOMAIN_HI / 2
        uni = np.linspace(0.1 * DOMAIN_HI, 0.9 * DOMAIN_HI, n)
        x = np.sort(0.7 * cheb + 0.3 * uni)
        w = np.ones(n) * (z_ref / n) * (1.0 + 0.04 * rng.standard_normal(n))

    return np.sort(x), _normalize_weights(w, z_ref)


def _hereditary_guess(
    prior_nodes: np.ndarray,
    prior_weights: np.ndarray,
    n: int,
    z_ref: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    pn = np.sort(np.asarray(prior_nodes, dtype=float).ravel())
    pw = np.maximum(np.asarray(prior_weights, dtype=float).ravel(), 1e-16)
    if len(pn) != n - 1 or len(pw) != n - 1:
        return generate_initial_guess("gauss_legendre", n, z_ref, rng)

    gaps = np.diff(pn)
    idx = int(np.argmax(gaps))
    new_node = 0.5 * (pn[idx] + pn[idx + 1])
    new_node += HEREDITARY_SIGMA_NODES * DOMAIN_HI * float(rng.standard_normal())
    new_node = float(np.clip(new_node, BOUNDARY_EPS, DOMAIN_HI - BOUNDARY_EPS))

    x = np.sort(np.concatenate([pn, [new_node]]))
    small_w = HEREDITARY_NEW_WEIGHT_FRAC * z_ref
    scale = max(z_ref - small_w, EPS) / max(float(np.sum(pw)), EPS)
    w = np.concatenate([pw * scale, [small_w]])
    w *= 1.0 + HEREDITARY_SIGMA_WEIGHTS * rng.standard_normal(n)
    return x, _normalize_weights(w, z_ref)
