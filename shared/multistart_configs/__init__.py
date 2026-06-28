"""Configuracoes de multi-start KSQ V1."""
from __future__ import annotations

from multistart_configs.ms_v1_reference import MS_V1_REFERENCE
from multistart_configs.ms_v2_stable_qmc import MS_V2_STABLE_QMC, MultistartConfigV2

DEFAULT_MULTISTART_CONFIG = MS_V1_REFERENCE

__all__ = [
    "DEFAULT_MULTISTART_CONFIG",
    "MS_V1_REFERENCE",
    "MS_V2_STABLE_QMC",
    "MultistartConfigV2",
]
