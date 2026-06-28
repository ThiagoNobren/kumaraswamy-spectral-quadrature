from __future__ import annotations

import numpy as np
from scipy.special import beta as beta_fn

from cases import DOMAIN_HI, DOMAIN_LO, KumaraswamyCase

EPS = 1e-30


def pdf(x: float | np.ndarray, case: KumaraswamyCase) -> float | np.ndarray:
    """PDF Kumaraswamy: f(x;a,b) = a*b*x^(a-1)*(1-x^a)^(b-1) em (0,1)."""
    a, b = float(case.a), float(case.b)
    x_arr = np.asarray(x, dtype=float)
    out = np.zeros_like(x_arr, dtype=float)
    mask = (x_arr > DOMAIN_LO) & (x_arr < DOMAIN_HI)
    xm = x_arr[mask]
    out[mask] = a * b * np.power(xm, a - 1.0) * np.power(1.0 - np.power(xm, a), b - 1.0)
    if np.ndim(x) == 0:
        return float(out)
    return out


def cdf(x: float | np.ndarray, case: KumaraswamyCase) -> float | np.ndarray:
    """CDF analitica: F(x;a,b) = 1 - (1 - x^a)^b."""
    a, b = float(case.a), float(case.b)
    x_arr = np.asarray(x, dtype=float)
    out = np.zeros_like(x_arr, dtype=float)
    out[x_arr <= DOMAIN_LO] = 0.0
    out[x_arr >= DOMAIN_HI] = 1.0
    mask = (x_arr > DOMAIN_LO) & (x_arr < DOMAIN_HI)
    xm = x_arr[mask]
    out[mask] = 1.0 - np.power(1.0 - np.power(xm, a), b)
    if np.ndim(x) == 0:
        return float(out)
    return out


def quantile(p: float, case: KumaraswamyCase) -> float:
    """Quantil analitico: Q(p;a,b) = [1 - (1-p)^(1/b)]^(1/a)."""
    if not (0.0 < float(p) < 1.0):
        raise ValueError("p deve pertencer a (0,1).")
    a, b = float(case.a), float(case.b)
    return float(np.power(1.0 - np.power(1.0 - float(p), 1.0 / b), 1.0 / a))


def raw_moment(r: int | float, case: KumaraswamyCase) -> float:
    """Momento bruto m_r = E[X^r] = b * B(1 + r/a, b)."""
    a, b = float(case.a), float(case.b)
    return float(b * beta_fn(1.0 + float(r) / a, b))
