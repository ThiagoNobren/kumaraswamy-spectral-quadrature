"""Parametros congelados da campanha do capitulo (K4R/K4S).

NAO EDITAR durante a campanha completa de simulacoes.
"""
from __future__ import annotations

# K4Q-mini — malha espectral
NX_SPECTRUM: int = 801

# K4P — faixa operacional segura
N_VALUES: list[int] = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

SOLVER: str = "KSQ_V1_Multiobjective"
STAGE1_CONFIG: str = "MS_V1_REFERENCE"
INTERPOLATION: str = "CubicSpline(bc_type=natural)"
K4R_REF: str = "docs/auditorias/auditoria_k4r_baseline_freeze.md"
