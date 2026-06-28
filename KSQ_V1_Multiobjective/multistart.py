from __future__ import annotations

import _bootstrap  # noqa: F401

import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

from multistart_base import MultistartConfig
from multistart_config import MS_V1_REFERENCE
from multistart_configs.ms_v2_stable_qmc import MultistartConfigV2
from multistart_jobs import CandidateJob, build_candidate_jobs

_ACTIVE_CONFIG = MS_V1_REFERENCE


@dataclass(frozen=True)
class CandidateOutcome:
    order: int
    label: str
    nodes: np.ndarray | None
    weights: np.ndarray | None
    residual_l2: float
    nfev: int
    success: bool


def get_active_multistart_config() -> MultistartConfig | MultistartConfigV2:
    return _ACTIVE_CONFIG


def set_active_multistart_config(config: MultistartConfig | MultistartConfigV2) -> None:
    global _ACTIVE_CONFIG
    _ACTIVE_CONFIG = config


def select_best_candidate(outcomes: list[CandidateOutcome]) -> CandidateOutcome | None:
    best: CandidateOutcome | None = None
    for outcome in sorted(outcomes, key=lambda o: o.order):
        if outcome.nodes is None or not np.isfinite(outcome.residual_l2):
            continue
        if best is None or outcome.residual_l2 < best.residual_l2:
            best = outcome
    return best


def _outcome_from_raw(job: CandidateJob, nodes, weights, resn: float, meta: dict) -> CandidateOutcome:
    return CandidateOutcome(
        order=job.order,
        label=job.label,
        nodes=None if nodes is None else np.asarray(nodes, dtype=float),
        weights=None if weights is None else np.asarray(weights, dtype=float),
        residual_l2=float(resn),
        nfev=int(meta.get("nfev", 0)),
        success=bool(meta.get("success", False)),
    )


def evaluate_candidates_sequential(
    jobs: list[CandidateJob],
    evaluate_one: Callable[[CandidateJob], tuple],
) -> list[CandidateOutcome]:
    outcomes: list[CandidateOutcome] = []
    for job in jobs:
        nodes, weights, resn, meta = evaluate_one(job)
        outcomes.append(_outcome_from_raw(job, nodes, weights, resn, meta))
    return outcomes


def _job_to_payload(job: CandidateJob) -> dict[str, Any]:
    return {
        "order": job.order,
        "label": job.label,
        "solve_strategy": job.solve_strategy,
        "seed": job.seed,
        "case_a": job.case_a,
        "case_b": job.case_b,
        "hereditary_prior_nodes": job.hereditary_prior_nodes,
        "hereditary_prior_weights": job.hereditary_prior_weights,
        "lobatto_use_exact_endpoints": job.lobatto_use_exact_endpoints,
    }


def _worker_evaluate(payload: dict[str, Any]) -> dict[str, Any]:
    project_dir = Path(payload["project_dir"])
    shared_dir = project_dir.parent / "shared"
    for p in (shared_dir, project_dir):
        s = str(p)
        if s not in sys.path:
            sys.path.insert(0, s)

    from optimization_stage1 import stage1_solve_with_strategy

    job = payload["job"]
    n_pts = int(payload["n_pts"])
    dx = float(payload["dx"])
    z_ref = float(payload["z_ref"])
    vecs = np.asarray(payload["vecs"], dtype=float)
    targets = np.asarray(payload["targets"], dtype=float)

    prior_nodes = job.get("hereditary_prior_nodes")
    prior_weights = job.get("hereditary_prior_weights")
    rng = np.random.default_rng(int(job["seed"]))
    nodes, weights, resn, meta = stage1_solve_with_strategy(
        n_pts,
        vecs,
        dx,
        targets,
        z_ref,
        job["solve_strategy"],
        rng=rng,
        case_a=job.get("case_a"),
        case_b=job.get("case_b"),
        hereditary_prior_nodes=None if prior_nodes is None else np.asarray(prior_nodes, dtype=float),
        hereditary_prior_weights=None if prior_weights is None else np.asarray(prior_weights, dtype=float),
        lobatto_use_exact_endpoints=bool(job.get("lobatto_use_exact_endpoints", True)),
        multistart_seed=int(job["seed"]),
    )
    return {
        "order": int(job["order"]),
        "label": job["label"],
        "nodes": None if nodes is None else np.asarray(nodes, dtype=float),
        "weights": None if weights is None else np.asarray(weights, dtype=float),
        "residual_l2": float(resn),
        "nfev": int(meta.get("nfev", 0)),
        "success": bool(meta.get("success", False)),
    }


def evaluate_candidates_parallel(
    jobs: list[CandidateJob],
    *,
    project_dir: str,
    n_pts: int,
    dx: float,
    z_ref: float,
    vecs: np.ndarray,
    targets: np.ndarray,
    max_workers: int | None = None,
) -> list[CandidateOutcome]:
    if not jobs:
        return []
    workers = max_workers if max_workers is not None else min(len(jobs), os.cpu_count() or 1)
    workers = max(1, min(workers, len(jobs)))
    payloads = [
        {
            "project_dir": project_dir,
            "job": _job_to_payload(j),
            "n_pts": n_pts,
            "dx": dx,
            "z_ref": z_ref,
            "vecs": np.asarray(vecs, dtype=float),
            "targets": np.asarray(targets, dtype=float),
        }
        for j in jobs
    ]
    raw_results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_worker_evaluate, p) for p in payloads]
        for fut in as_completed(futures):
            raw_results.append(fut.result())
    outcomes: list[CandidateOutcome] = []
    for r in raw_results:
        outcomes.append(
            CandidateOutcome(
                order=int(r["order"]),
                label=str(r["label"]),
                nodes=r["nodes"],
                weights=r["weights"],
                residual_l2=float(r["residual_l2"]),
                nfev=int(r["nfev"]),
                success=bool(r["success"]),
            )
        )
    return outcomes


__all__ = [
    "CandidateJob",
    "CandidateOutcome",
    "build_candidate_jobs",
    "evaluate_candidates_parallel",
    "evaluate_candidates_sequential",
    "get_active_multistart_config",
    "select_best_candidate",
    "set_active_multistart_config",
]
