from __future__ import annotations

from multistart_base import MultistartConfig
from multistart_configs.ms_v1_reference import MS_V1_REFERENCE

# Alias retrocompativel — aponta para MS_V1_REFERENCE congelado.
KSQ_MS_V1 = MS_V1_REFERENCE
DEFAULT_MULTISTART_CONFIG = MS_V1_REFERENCE

__all__ = [
    "DEFAULT_MULTISTART_CONFIG",
    "KSQ_MS_V1",
    "MS_V1_REFERENCE",
    "MultistartConfig",
]
