"""MS_V1_REFERENCE — configuracao congelada do multi-start usada em K1–K4.

IMUTAVEL: nao editar valores apos congelamento. Nova experimentacao em ms_v2_stable_qmc.py.
"""
from __future__ import annotations

from types import MappingProxyType

from multistart_base import MultistartConfig

_STRATEGY_ID = MappingProxyType(
    {
        "gauss_legendre": 0,
        "chebyshev": 1,
        "hibrido": 2,
        "uniforme": 3,
        "regularizacao_minima": 4,
    }
)

MS_V1_REFERENCE = MultistartConfig(
    name="MS_V1_REFERENCE",
    spectral_strategies=(
        "gauss_legendre",
        "chebyshev",
        "hibrido",
        "uniforme",
        "regularizacao_minima",
    ),
    uniform_extra_enabled=True,
    uniform_extra_trials=20,
    uniform_extra_solve_strategy="uniforme",
    strategy_id=_STRATEGY_ID,
    base_seed_fixed=1234,
    base_seed_uniform_extra=91000,
    fixed_strategy_trial=0,
    parallel_default=True,
)
