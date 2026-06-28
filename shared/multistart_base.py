from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True)
class MultistartConfig:
    """Configuracao V1: estrategias fixas + uniforme_extra (baseline K1–K4)."""

    name: str
    spectral_strategies: tuple[str, ...]
    uniform_extra_enabled: bool
    uniform_extra_trials: int
    uniform_extra_solve_strategy: str
    strategy_id: MappingProxyType
    base_seed_fixed: int
    base_seed_uniform_extra: int
    fixed_strategy_trial: int
    parallel_default: bool

    def candidate_seed(self, n_pts: int, strategy_name: str, trial: int | None = None) -> int:
        t = self.fixed_strategy_trial if trial is None else int(trial)
        sid = self.strategy_id[strategy_name]
        return int(self.base_seed_fixed + n_pts + 1000 * t + sid)

    def uniform_extra_seed(self, n_pts: int, trial: int) -> int:
        return int(self.base_seed_uniform_extra + 1000 * n_pts + int(trial))

    def uniform_extra_label(self, trial: int) -> str:
        return f"uniforme_extra_{int(trial)}"
