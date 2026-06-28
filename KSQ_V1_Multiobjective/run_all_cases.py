from __future__ import annotations

import _bootstrap  # noqa: F401

import argparse
import csv
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from analytic_qoi import analytic_qoi
from cases import N_VALUES, QUANTILE_LEVELS, all_cases
from optimization_stage1 import stage1_solve
from optimization_stage2 import stage2_multiobjective_refine
from qoi_metrics import compute_qoi_metrics
from reporting import run_with_tee
from spectral_neumann import setup_spectral_bundle

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
REPORT_DIR = SCRIPT_DIR / "relatorio_cli"
OUTPUT_CSV = ROOT_DIR / "outputs" / "csv" / "k2_multi_results.csv"

RANDOM_UNIFORM_TRIALS = 20


@dataclass(frozen=True)
class MultiRunRow:
    case_id: str
    label: str
    n: int
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


def run_one(case, n_nodes: int, *, parallel: bool = True) -> MultiRunRow:
    ref = analytic_qoi(case)
    z_ref, dx, vecs, targets, _ = setup_spectral_bundle(case, n_nodes)
    t0 = time.perf_counter()
    stage1 = stage1_solve(
        n_nodes, vecs, dx, targets, z_ref, random_uniform_trials=RANDOM_UNIFORM_TRIALS, parallel=parallel
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
    return MultiRunRow(
        case_id=case.case_id,
        label=case.label,
        n=n_nodes,
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
    )


def _write_csv(rows: list[MultiRunRow]) -> Path:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", encoding="utf-8", newline="\n") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    return OUTPUT_CSV


def cli_main(*, n_values: list[int] | None = None, parallel: bool = True) -> list[MultiRunRow]:
    n_list = N_VALUES if n_values is None else n_values
    print("KSQ V1 Multiobjective — todos os casos")
    rows: list[MultiRunRow] = []
    for case in all_cases():
        print(f"\n=== {case.label} ({case.case_id}) ===")
        for n_nodes in n_list:
            row = run_one(case, n_nodes, parallel=parallel)
            rows.append(row)
            print(
                f"  N={n_nodes:2d}  rel_mean={row.rel_mean:.3e}  rel_var={row.rel_var:.3e}  "
                f"rel_q50={row.rel_q50:.3e}  F={row.F:.3e}  ||r||={row.residual_l2:.3e}  "
                f"stg1={row.strategy_stage1}  t={row.elapsed_s:.2f}s"
            )
    path = _write_csv(rows)
    print(f"\nCSV gravado: {path}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="KSQ V1 Multiobjective — run all cases")
    parser.add_argument("--serial", action="store_true")
    parser.add_argument("--n", type=int, nargs="*", default=None)
    args = parser.parse_args()
    run_with_tee(
        lambda: cli_main(n_values=args.n, parallel=not args.serial),
        REPORT_DIR,
        "saida_terminal_ksq_v1_multiobjective",
    )


if __name__ == "__main__":
    main()
