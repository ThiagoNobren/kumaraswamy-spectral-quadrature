from __future__ import annotations

from dataclasses import dataclass

from cases import QUANTILE_LEVELS, KumaraswamyCase
from kumaraswamy import cdf, pdf, quantile, raw_moment


@dataclass(frozen=True)
class AnalyticQoI:
    mean: float
    variance: float
    m2: float
    quantiles: dict[str, float]

    def quantile(self, p: float) -> float:
        key = _quant_key(p)
        return float(self.quantiles[key])


def _quant_key(p: float) -> str:
    pct = int(round(100.0 * float(p)))
    return f"q{pct}"


def analytic_qoi(case: KumaraswamyCase, *, quantile_levels: list[float] | None = None) -> AnalyticQoI:
    """QoI analiticas: media, variancia e quantis da Kumaraswamy pura."""
    levels = QUANTILE_LEVELS if quantile_levels is None else quantile_levels
    mu = raw_moment(1, case)
    m2 = raw_moment(2, case)
    var = float(m2 - mu * mu)
    qmap = {_quant_key(p): quantile(p, case) for p in levels}
    return AnalyticQoI(mean=float(mu), variance=float(var), m2=float(m2), quantiles=qmap)


def analytic_cdf_on_grid(case: KumaraswamyCase, grid_size: int) -> tuple[list[float], list[float]]:
    import numpy as np

    grid = np.linspace(0.0, 1.0, int(grid_size), dtype=float)
    vals = [float(cdf(float(t), case)) for t in grid]
    return [float(t) for t in grid], vals
