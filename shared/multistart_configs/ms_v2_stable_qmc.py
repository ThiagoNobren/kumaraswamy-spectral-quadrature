"""MS_V2_STABLE_QMC — configuracao experimental (nao substitui MS_V1_REFERENCE)."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

_V2_STRATEGY_ID = MappingProxyType(
    {
        "gauss_legendre": 0,
        "chebyshev": 1,
        "gauss_jacobi_approx": 2,
        "gauss_lobatto": 3,
        "halton_sequence": 20,
        "latin_hypercube": 30,
        "hereditary_N_minus_1": 40,
    }
)


@dataclass(frozen=True)
class MultistartConfigV2:
    """Multi-start experimental: 4 deterministicos + 15 QMC/LHS + 1 hereditario."""

    name: str
    deterministic_strategies: tuple[str, ...]
    qmc_trials: tuple[tuple[str, int], ...]
    hereditary_strategy: str
    hereditary_fallback: str
    strategy_id: MappingProxyType
    base_seed_fixed: int
    base_seed_qmc: int
    fixed_strategy_trial: int
    parallel_default: bool
    lobatto_use_exact_endpoints: bool

    def candidate_seed(self, n_pts: int, strategy_name: str, trial: int | None = None) -> int:
        t = self.fixed_strategy_trial if trial is None else int(trial)
        sid = self.strategy_id[strategy_name]
        return int(self.base_seed_fixed + n_pts + 1000 * t + sid)

    def qmc_seed(self, n_pts: int, strategy_name: str, trial: int) -> int:
        sid = self.strategy_id[strategy_name]
        return int(self.base_seed_qmc + 1000 * n_pts + 100 * sid + int(trial))

    def qmc_label(self, strategy_name: str, trial: int) -> str:
        short = {"halton_sequence": "halton", "latin_hypercube": "lhs"}[strategy_name]
        return f"{short}_{int(trial)}"

    @property
    def total_candidates(self) -> int:
        qmc_count = sum(c for _, c in self.qmc_trials)
        hereditary = 1
        return len(self.deterministic_strategies) + qmc_count + hereditary


MS_V2_STABLE_QMC = MultistartConfigV2(
    name="MS_V2_STABLE_QMC",
    deterministic_strategies=(
        "gauss_legendre",
        "chebyshev",
        "gauss_jacobi_approx",
        "gauss_lobatto",
    ),
    qmc_trials=(
        ("halton_sequence", 8),
        ("latin_hypercube", 7),
    ),
    hereditary_strategy="hereditary_N_minus_1",
    hereditary_fallback="gauss_legendre",
    strategy_id=_V2_STRATEGY_ID,
    base_seed_fixed=2234,
    base_seed_qmc=92000,
    fixed_strategy_trial=0,
    parallel_default=True,
    lobatto_use_exact_endpoints=False,
)
