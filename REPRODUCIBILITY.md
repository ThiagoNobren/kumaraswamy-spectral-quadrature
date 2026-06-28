# Reproducibility Guide

This document describes how to reproduce the numerical results and figure presented in the Springer chapter, starting from the frozen baseline **`v1.0-article-baseline`**.

## Requirements

- **Python** 3.8 or newer
- **OS**: Windows, Linux, or macOS (use `--serial` on Windows for the campaign)
- Dependencies listed in `requirements.txt`

## Installation

From the repository root:

```bash
python -m venv .venv
```

Activate the virtual environment:

```bash
# Windows
.venv\Scripts\activate

# Linux / macOS
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Official campaign

The chapter results were produced with the dedicated campaign runner (not `run_all_cases.py`).

```bash
python campaign/run_chapter_campaign.py --serial
```

### Frozen parameters

Defined in `campaign/chapter_campaign_config.py`:

| Parameter | Value |
|-----------|-------|
| `NX_SPECTRUM` | 801 |
| `N_VALUES` | 4 … 14 |
| `SOLVER` | `KSQ_V1_Multiobjective` |
| `STAGE1_CONFIG` | `MS_V1_REFERENCE` |
| `INTERPOLATION` | `CubicSpline(bc_type=natural)` |

### Expected outputs

| Artefact | Path |
|----------|------|
| Results CSV | `outputs/csv/chapter_campaign_multi_results.csv` |
| Terminal log | `relatorio_cli/chapter_campaign_multi_YYYYMMDD_HHMMSS.md` |

The reference CSV contains **33 data rows** (3 Kumaraswamy cases × 11 quadrature orders) plus a header.

Key columns: `rel_mean`, `rel_var`, `rel_q50`, `n`, `nx_spectrum`, `stage1_config`, `solver`.

### Runtime

A full serial campaign typically requires **1–2 hours** depending on hardware. Order `N = 14` is the most expensive per case.

## Chapter figure

Generate the Springer figure from the official CSV:

```bash
python figures/make_relative_errors_kumaraswamy.py
```

### Expected outputs

| Artefact | Path |
|----------|------|
| Vector figure | `outputs/figures/relative_errors_kumaraswamy.pdf` |
| Raster figure | `outputs/figures/relative_errors_kumaraswamy.png` |

The script reads only `outputs/csv/chapter_campaign_multi_results.csv` and plots **Mean**, **Variance**, and **Q_0.50** relative errors versus quadrature order `N` for the three Kumaraswamy cases.

## Validating against reference artefacts

After a campaign run, compare your CSV with the bundled reference:

```bash
# Example: row count (33 data rows)
python -c "import csv; print(sum(1 for _ in csv.DictReader(open('outputs/csv/chapter_campaign_multi_results.csv'))))"
```

Reference checksums are not enforced by this package; byte-identical reproduction depends on platform floating-point behaviour. Relative errors and campaign structure should match at reported precision.

## Baseline freeze records

For provenance of the frozen configuration:

- **K4R** — `docs/auditorias/auditoria_k4r_baseline_freeze.md`
- **K4S** — `docs/auditorias/auditoria_k4s_campaign_freeze.md`
- **Official campaign log** — `relatorio_cli/chapter_campaign_multi_20260626_001524.md`

## What is not included

This package intentionally excludes exploratory audits, the Simple solver, diagnostic CSV archives, and intermediate campaign attempts. Those belong to the private research repository.
