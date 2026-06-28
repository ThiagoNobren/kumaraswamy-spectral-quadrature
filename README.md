# Kumaraswamy Spectral Quadrature — Article Baseline

This repository provides the **reference implementation** used to produce the computational results of the submitted Springer book chapter on spectral Neumann quadrature for Kumaraswamy distributions.

The present release corresponds **exactly** to the frozen scientific baseline employed at submission time (baseline tag `v1.0-article-baseline`, development commit `c021af5`). Its purpose is to enable independent reproduction of the numerical experiments reported in the chapter. It does **not** represent the full research environment, the complete development history, or ongoing exploratory work carried out in the private research repository.

## Scientific scope

The chapter reports results from the **KSQ_V1_Multiobjective** solver with:

| Parameter | Value |
|-----------|-------|
| Stage 1 multistart | `MS_V1_REFERENCE` |
| Eigenfunction interpolation | natural cubic spline (`CubicSpline`, `bc_type="natural"`) |
| Spectral mesh | `nx = 801` |
| Quadrature orders | `N = 4, …, 14` |
| Cases | Kum(2,8), Kum(2,2.41), Kum(4,1.5) |

Quantities of interest: mean, variance, and median (`Q_0.50`) relative errors, plus the chapter figure of errors versus `N`.

## Project structure

```
.
├── campaign/                  # Official chapter campaign entrypoint
├── figures/                   # Figure generation script
├── KSQ_V1_Multiobjective/     # Multi-objective spectral solver
├── shared/                    # Kumaraswamy model, spectral Neumann, QoI metrics
├── docs/auditorias/           # Baseline freeze documentation (K4R, K4S)
├── outputs/
│   ├── csv/                   # Official campaign results (reference)
│   └── figures/               # Chapter figure (PDF + PNG)
└── relatorio_cli/             # Official campaign terminal log
```

## Quick start

See **[REPRODUCIBILITY.md](REPRODUCIBILITY.md)** for full installation and reproduction steps.

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt

# Official campaign (serial mode recommended on Windows)
python campaign/run_chapter_campaign.py --serial

# Chapter figure from the official CSV
python figures/make_relative_errors_kumaraswamy.py
```

## Reference outputs

Pre-computed reference artefacts shipped with this package:

- `outputs/csv/chapter_campaign_multi_results.csv`
- `outputs/figures/relative_errors_kumaraswamy.pdf`
- `outputs/figures/relative_errors_kumaraswamy.png`
- `relatorio_cli/chapter_campaign_multi_20260626_001524.md`

## Baseline documentation

- `docs/auditorias/auditoria_k4r_baseline_freeze.md` — scientific baseline freeze (K4R)
- `docs/auditorias/auditoria_k4s_campaign_freeze.md` — campaign infrastructure freeze (K4S)

## Citation

This repository contains the reference implementation accompanying the submitted book chapter reporting the generalized spectral quadrature methodology and the numerical experiments reproduced herein.

Until the chapter is officially published, the full bibliographic reference is not yet available. In the interim, the chapter may be cited as a **submitted book chapter**. After publication, this section will be updated with the definitive bibliographic reference, the DOI, and the official BibTeX entry of the published work.

```bibtex
@unpublished{NascimentoCostaSolheid2026,
  author  = {Nascimento, Antonio Thiago Nobre and Costa, Eliardo Guimar{\~a}es da and Solheid, Bruno dos Santos},
  title   = {Generalized Spectral Quadrature for Statistical Computation},
  year    = {2026},
  note    = {Submitted book chapter.}
}
```

## License

This software is distributed under the terms of the **GNU General Public License, version 3**. See [LICENSE](LICENSE) for the full text.

## Scope of this release

This public package contains **only** the material required to reproduce the chapter results. It does not include the broader research history, the Simple solver, exploratory audits, or diagnostic datasets from the private development repository.
