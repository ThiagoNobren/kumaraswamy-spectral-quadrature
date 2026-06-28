from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh

from cases import DOMAIN_HI, KumaraswamyCase
from kumaraswamy import pdf

DEFAULT_NX_SPECTRUM = 401
EPS = 1e-30


def create_neumann_matrix(nx: int, dx: float) -> np.ndarray:
    mat = np.zeros((nx, nx), dtype=float)
    for i in range(nx):
        mat[i, i] = 2.0
    for i in range(nx - 1):
        mat[i, i + 1] = mat[i + 1, i] = -1.0
    mat[0, 0] = 1.0
    mat[-1, -1] = 1.0
    return mat * (1.0 / dx**2)


def _normalize_columns_uniform_l2(vecs: np.ndarray, dx: float) -> np.ndarray:
    out = np.array(vecs, dtype=float, copy=True)
    for j in range(out.shape[1]):
        v = out[:, j]
        nrm = float(np.sqrt(np.sum(v * v) * dx))
        if nrm > 0:
            out[:, j] = v / nrm
    return out


@lru_cache(maxsize=32)
def _build_neumann_spectrum_cached(nx: int, ldom: float) -> tuple[float, np.ndarray, np.ndarray]:
    if nx < 4:
        raise ValueError("nx deve ser pelo menos 4.")
    dx = float(ldom) / float(nx - 1)
    mat = create_neumann_matrix(nx, dx)
    lam, vecs = eigh(mat)
    vecs = _normalize_columns_uniform_l2(vecs, dx)
    vecs[:, 0] = 1.0
    return dx, lam, vecs


def build_neumann_spectrum(nx: int = DEFAULT_NX_SPECTRUM, *, ldom: float = DOMAIN_HI) -> tuple[float, np.ndarray, np.ndarray]:
    return _build_neumann_spectrum_cached(int(nx), float(ldom))


_SPLINE_CACHE: dict[tuple[int, int, int, float], tuple[np.ndarray, list[CubicSpline]]] = {}


def _natural_splines(vecs: np.ndarray, dx: float) -> tuple[np.ndarray, list[CubicSpline]]:
    key = (
        int(vecs.__array_interface__["data"][0]),
        int(vecs.shape[0]),
        int(vecs.shape[1]),
        round(float(dx), 15),
    )
    cached = _SPLINE_CACHE.get(key)
    if cached is not None:
        return cached
    nx = int(vecs.shape[0])
    x_grid = np.linspace(0.0, (nx - 1) * float(dx), nx)
    m = int(vecs.shape[1])
    splines = [
        CubicSpline(x_grid, vecs[:, j], bc_type="natural", extrapolate=False) for j in range(m)
    ]
    cached = (x_grid, splines)
    _SPLINE_CACHE[key] = cached
    if len(_SPLINE_CACHE) > 32:
        _SPLINE_CACHE.pop(next(iter(_SPLINE_CACHE)))
    return cached


def phi_matrix(nodes: np.ndarray, vecs: np.ndarray, dx: float) -> np.ndarray:
    """Matriz G[i,j] = phi_i(theta_j), shape (m_modos, n_nos)."""
    x_grid, splines = _natural_splines(vecs, dx)
    th = np.clip(np.asarray(nodes, dtype=float).ravel(), float(x_grid[0]), float(x_grid[-1]))
    m = int(vecs.shape[1])
    values = np.column_stack([splines[k](th) for k in range(m)])
    return np.asarray(values.T, dtype=float)


def _phi_at_point(x: float, vecs: np.ndarray, dx: float, k: int) -> float:
    x_grid, splines = _natural_splines(vecs, dx)
    x_clip = float(np.clip(x, float(x_grid[0]), float(x_grid[-1])))
    return float(splines[k](x_clip))


def spectral_targets(case: KumaraswamyCase, vecs: np.ndarray, dx: float) -> np.ndarray:
    """Alvos espectrais t_k = int phi_k(x) f(x) dx, k=0..m-1. t_0 = 1 (PDF normalizada)."""
    m = int(vecs.shape[1])
    targets = np.zeros(m, dtype=float)
    targets[0] = 1.0
    quad_points = [1e-12, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0 - 1e-12]

    for i in range(1, m):

        def integrand(tt: float, idx: int = i) -> float:
            return _phi_at_point(tt, vecs, dx, idx) * float(pdf(tt, case))

        val, _ = quad(integrand, 0.0, DOMAIN_HI, limit=500, points=quad_points, epsabs=1e-12, epsrel=1e-12)
        targets[i] = val
    return targets


def spectral_moment_residuals(
    nodes: np.ndarray,
    weights: np.ndarray,
    targets: np.ndarray,
    vecs: np.ndarray,
    dx: float,
) -> np.ndarray:
    g_mat = phi_matrix(nodes, vecs, dx)
    m = int(len(targets))
    if g_mat.shape[0] > m:
        g_mat = g_mat[:m, :]
    return np.asarray(g_mat @ np.asarray(weights, dtype=float) - np.asarray(targets, dtype=float), dtype=float)


def setup_spectral_bundle(
    case: KumaraswamyCase,
    n_nodes: int,
    *,
    nx_spectrum: int = DEFAULT_NX_SPECTRUM,
) -> tuple[float, float, np.ndarray, np.ndarray, np.ndarray]:
    z_ref = 1.0
    dx, _lam, vecs_full = build_neumann_spectrum(nx_spectrum, ldom=DOMAIN_HI)
    m = 2 * int(n_nodes)
    vecs = vecs_full[:, :m]
    targets = spectral_targets(case, vecs, dx)
    return z_ref, dx, vecs, targets, vecs_full
