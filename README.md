# Kumaraswamy Spectral Quadrature — Article Baseline

This repository provides the **reference implementation** used to produce the computational results reported in the chapter *A Spectral Quadrature for Computing Moments and Quantiles of Continuous Probability Distributions*, submitted to the Springer Nature contributed volume *Statistical Modeling and Analysis in the Age of AI* (LISA 2020 Global Network).

The repository corresponds to the frozen computational baseline used to generate all numerical experiments reported in the submitted chapter. Its purpose is to enable the independent reproduction of the computational results presented in the chapter. It does **not** represent the complete research environment, the full development history, or ongoing exploratory work carried out in the private research repository.

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

This repository contains the reference implementation accompanying the submitted chapter *A Spectral Quadrature for Computing Moments and Quantiles of Continuous Probability Distributions* and provides all computational resources required to reproduce the numerical experiments reported therein.

Since the chapter is currently under editorial review, the definitive bibliographic reference is not yet available. The BibTeX entry below may be used to cite the submitted manuscript. After publication, this section will be updated with the final bibliographic information, including the DOI and the official citation provided by Springer Nature.

```bibtex
@unpublished{NascimentoCostaSolheid2026,
  author = {Nascimento, Antonio Thiago Nobre and Costa, Eliardo Guimar{\~a}es da and Solheid, Bruno dos Santos},
  title  = {A Spectral Quadrature for Computing Moments and Quantiles of Continuous Probability Distributions},
  year   = {2026},
  note   = {Submitted chapter to the Springer Nature contributed volume Statistical Modeling and Analysis in the Age of AI (LISA 2020 Global Network).}
}
```

## License

This software is distributed under the terms of the **GNU General Public License, version 3**. See [LICENSE](LICENSE) for the full text.

## Scope of this release

This public package contains **only** the material required to reproduce the chapter results. It does not include the broader research history, the Simple solver, exploratory audits, or diagnostic datasets from the private development repository.
