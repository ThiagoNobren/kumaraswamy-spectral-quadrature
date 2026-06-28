"""Construcao de jobs de multi-start (V1 reference + V2 experimental)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from multistart_base import MultistartConfig
from multistart_configs.ms_v2_stable_qmc import MultistartConfigV2


@dataclass(frozen=True)
class CandidateJob:
    order: int
    label: str
    solve_strategy: str
    seed: int
    case_a: float | None = None
    case_b: float | None = None
    hereditary_prior_nodes: tuple[float, ...] | None = None
    hereditary_prior_weights: tuple[float, ...] | None = None
    lobatto_use_exact_endpoints: bool = True


def is_v2_config(config: Any) -> bool:
    return isinstance(config, MultistartConfigV2)


def build_candidate_jobs(
    config: MultistartConfig | MultistartConfigV2,
    n_pts: int,
    random_uniform_trials: int | None = None,
    *,
    hereditary_prior: tuple[np.ndarray, np.ndarray] | None = None,
    case_a: float | None = None,
    case_b: float | None = None,
) -> list[CandidateJob]:
    if is_v2_config(config):
        jobs = _build_v2_jobs(config, n_pts, hereditary_prior=hereditary_prior)
    else:
        jobs = _build_v1_jobs(config, n_pts, random_uniform_trials)
    if case_a is None and case_b is None:
        return jobs
    return [
        CandidateJob(
            order=j.order,
            label=j.label,
            solve_strategy=j.solve_strategy,
            seed=j.seed,
            case_a=case_a,
            case_b=case_b,
            hereditary_prior_nodes=j.hereditary_prior_nodes,
            hereditary_prior_weights=j.hereditary_prior_weights,
            lobatto_use_exact_endpoints=j.lobatto_use_exact_endpoints,
        )
        for j in jobs
    ]


def _build_v1_jobs(
    cfg: MultistartConfig,
    n_pts: int,
    random_uniform_trials: int | None,
) -> list[CandidateJob]:
    trials = cfg.uniform_extra_trials if random_uniform_trials is None else int(random_uniform_trials)
    jobs: list[CandidateJob] = []
    order = 0
    for strat in cfg.spectral_strategies:
        jobs.append(
            CandidateJob(
                order=order,
                label=strat,
                solve_strategy=strat,
                seed=cfg.candidate_seed(n_pts, strat, cfg.fixed_strategy_trial),
            )
        )
        order += 1
    if cfg.uniform_extra_enabled:
        for k in range(trials):
            jobs.append(
                CandidateJob(
                    order=order,
                    label=cfg.uniform_extra_label(k),
                    solve_strategy=cfg.uniform_extra_solve_strategy,
                    seed=cfg.uniform_extra_seed(n_pts, k),
                )
            )
            order += 1
    return jobs


def _build_v2_jobs(
    cfg: MultistartConfigV2,
    n_pts: int,
    *,
    hereditary_prior: tuple[np.ndarray, np.ndarray] | None,
) -> list[CandidateJob]:
    jobs: list[CandidateJob] = []
    order = 0
    for strat in cfg.deterministic_strategies:
        jobs.append(
            CandidateJob(
                order=order,
                label=strat,
                solve_strategy=strat,
                seed=cfg.candidate_seed(n_pts, strat, cfg.fixed_strategy_trial),
                lobatto_use_exact_endpoints=cfg.lobatto_use_exact_endpoints,
            )
        )
        order += 1

    for strat, count in cfg.qmc_trials:
        for k in range(int(count)):
            jobs.append(
                CandidateJob(
                    order=order,
                    label=cfg.qmc_label(strat, k),
                    solve_strategy=strat,
                    seed=cfg.qmc_seed(n_pts, strat, k),
                )
            )
            order += 1

    h_label = cfg.hereditary_strategy
    h_strat = cfg.hereditary_strategy
    if hereditary_prior is None:
        h_strat = cfg.hereditary_fallback
        h_label = f"{cfg.hereditary_strategy}_fallback"
    prior_nodes = None
    prior_weights = None
    if hereditary_prior is not None:
        prior_nodes = tuple(float(v) for v in np.asarray(hereditary_prior[0], dtype=float).ravel())
        prior_weights = tuple(float(v) for v in np.asarray(hereditary_prior[1], dtype=float).ravel())

    jobs.append(
        CandidateJob(
            order=order,
            label=h_label,
            solve_strategy=h_strat,
            seed=cfg.candidate_seed(n_pts, cfg.hereditary_strategy, cfg.fixed_strategy_trial),
            hereditary_prior_nodes=prior_nodes,
            hereditary_prior_weights=prior_weights,
            lobatto_use_exact_endpoints=cfg.lobatto_use_exact_endpoints,
        )
    )
    return jobs
