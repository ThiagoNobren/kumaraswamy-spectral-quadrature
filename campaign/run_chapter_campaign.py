"""Campanha completa do capitulo — baseline K4R congelada (K4S).

Entrypoint oficial da campanha numerica. Nao usar run_all_cases.py para o capitulo.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

_CAMPAIGN = Path(__file__).resolve().parent
_ROOT = _CAMPAIGN.parent
_MULTI = _ROOT / "KSQ_V1_Multiobjective"
_SHARED = _ROOT / "shared"
for p in (_CAMPAIGN, _SHARED, _MULTI):
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)

import _bootstrap  # noqa: F401

from analytic_qoi import analytic_qoi
from cases import QUANTILE_LEVELS, all_cases
from chapter_campaign_config import (
    INTERPOLATION,
    K4R_REF,
    N_VALUES as CAMPAIGN_N_VALUES,
    NX_SPECTRUM as CAMPAIGN_NX,
    SOLVER,
    STAGE1_CONFIG,
)
from multistart_config import MS_V1_REFERENCE
from optimization_stage1 import stage1_solve
from optimization_stage2 import stage2_multiobjective_refine
from qoi_metrics import compute_qoi_metrics
from reporting import run_with_tee
from spectral_neumann import setup_spectral_bundle

REPORT_DIR = _ROOT / "relatorio_cli"
OUTPUT_CSV = _ROOT / "outputs" / "csv" / "chapter_campaign_multi_results.csv"

RANDOM_UNIFORM_TRIALS = 20


@dataclass(frozen=True)
class CampaignRow:
    case_id: str
    label: str
    n: int
    nx_spectrum: int
    rel_mean: float
    rel_var: float
    rel_q25: float
    rel_q50: float
    rel_q75: float
    abs_q25: float
    abs_q50: float
    abs_q75: float
    CDF_sup: float
    W1: float
    residual_l2: float
    F: float
    elapsed_s: float
    strategy_stage1: str
    stage1_config: str
    solver: str


def run_one(case, n_nodes: int, *, nx_spectrum: int, parallel: bool) -> CampaignRow:
    ref = analytic_qoi(case)
    z_ref, dx, vecs, targets, _ = setup_spectral_bundle(
        case, n_nodes, nx_spectrum=nx_spectrum,
    )
    t0 = time.perf_counter()
    stage1 = stage1_solve(
        n_nodes, vecs, dx, targets, z_ref,
        random_uniform_trials=RANDOM_UNIFORM_TRIALS,
        parallel=parallel,
        multistart_config=MS_V1_REFERENCE,
        case_a=case.a,
        case_b=case.b,
    )
    stage2 = stage2_multiobjective_refine(
        case,
        stage1.nodes,
        stage1.weights,
        targets=targets,
        vecs=vecs,
        dx=dx,
        z_ref=z_ref,
    )
    elapsed = time.perf_counter() - t0
    metrics = compute_qoi_metrics(case, stage2.nodes, stage2.weights, ref, quantile_levels=QUANTILE_LEVELS)
    return CampaignRow(
        case_id=case.case_id,
        label=case.label,
        n=n_nodes,
        nx_spectrum=nx_spectrum,
        rel_mean=metrics["rel_mean"],
        rel_var=metrics["rel_var"],
        rel_q25=metrics["rel_q25"],
        rel_q50=metrics["rel_q50"],
        rel_q75=metrics["rel_q75"],
        abs_q25=metrics["abs_q25"],
        abs_q50=metrics["abs_q50"],
        abs_q75=metrics["abs_q75"],
        CDF_sup=metrics["CDF_sup"],
        W1=metrics["W1"],
        residual_l2=stage1.residual_l2,
        F=stage2.score_F,
        elapsed_s=elapsed,
        strategy_stage1=stage1.strategy,
        stage1_config=STAGE1_CONFIG,
        solver=SOLVER,
    )


def _write_csv(rows: list[CampaignRow]) -> Path:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    return OUTPUT_CSV


def cli_main(
    *,
    n_values: list[int] | None = None,
    nx_spectrum: int | None = None,
    parallel: bool = True,
) -> list[CampaignRow]:
    n_list = CAMPAIGN_N_VALUES if n_values is None else n_values
    nx = CAMPAIGN_NX if nx_spectrum is None else nx_spectrum
    print(f"Campanha capitulo — {SOLVER} | {STAGE1_CONFIG} | {INTERPOLATION}")
    print(f"Baseline: {K4R_REF}")
    print(f"nx_spectrum={nx}  N={n_list[0]}..{n_list[-1]}  ({len(n_list)} ordens)")
    rows: list[CampaignRow] = []
    for case in all_cases():
        print(f"\n=== {case.label} ({case.case_id}) ===")
        for n_nodes in n_list:
            row = run_one(case, n_nodes, nx_spectrum=nx, parallel=parallel)
            rows.append(row)
            print(
                f"  N={n_nodes:2d}  nx={row.nx_spectrum}  rel_mean={row.rel_mean:.3e}  "
                f"rel_var={row.rel_var:.3e}  rel_q50={row.rel_q50:.3e}  "
                f"F={row.F:.3e}  ||r||={row.residual_l2:.3e}  "
                f"stg1={row.strategy_stage1}  t={row.elapsed_s:.2f}s"
            )
    path = _write_csv(rows)
    print(f"\nCSV gravado: {path}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Campanha capitulo — baseline K4R/K4S")
    parser.add_argument("--serial", action="store_true", help="desliga paralelismo Stage1")
    parser.add_argument(
        "--n", type=int, nargs="*", default=None,
        help="sobrescreve N (uso de auditoria apenas; defeito = congelado K4R)",
    )
    parser.add_argument(
        "--nx", type=int, default=None,
        help="sobrescreve nx (uso de auditoria apenas; defeito = 801)",
    )
    args = parser.parse_args()
    run_with_tee(
        lambda: cli_main(
            n_values=args.n,
            nx_spectrum=args.nx,
            parallel=not args.serial,
        ),
        REPORT_DIR,
        "chapter_campaign_multi",
    )


if __name__ == "__main__":
    main()
