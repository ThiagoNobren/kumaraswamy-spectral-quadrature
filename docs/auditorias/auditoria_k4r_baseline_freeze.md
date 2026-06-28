# K4R — Baseline Freeze (gate campanha completa)

**Data:** 2026-06-25 21:41  
**Branch:** `feature/natural-cubic-eigenfunctions`  
**Commit:** `c8fde528eba95fff1a39998d28711cd20810932d`  
**CSV:** `outputs/csv/k4r_baseline_freeze.csv`

## Baseline congelada

| componente | valor |
| --- | --- |
| Solver | `KSQ_V1_Multiobjective` |
| Stage1 | `MS_V1_REFERENCE` |
| Interpolação φ | `CubicSpline(..., bc_type="natural")` |
| Malha espectral | **nx = 801** (K4Q-mini) |
| Ordens | **N = 4..14** (K4P) |
| Casos | Kum(2,8), Kum(2,2.41), Kum(4,1.5) |

**Sem merge** em `main`. **Sem alteracao** de algoritmo nesta auditoria.

## Verificacoes

| categoria | item | esperado | observado | status |
| --- | --- | --- | --- | --- |
| reprodutibilidade | git_hash | (registado) | c8fde528eba9 | PASS |
| reprodutibilidade | branch | feature/natural-cubic-eigenfunctions | feature/natural-cubic-eigenfunctions | PASS |
| reprodutibilidade | data_utc | (ISO) | 2026-06-26 00:41 UTC | PASS |
| solver | entrypoint | KSQ_V1_Multiobjective | KSQ_V1_Multiobjective | PASS |
| stage1 | multistart ativo | MS_V1_REFERENCE | MS_V1_REFERENCE | PASS |
| stage1 | DEFAULT_MULTISTART_CONFIG | MS_V1_REFERENCE | MS_V1_REFERENCE | PASS |
| interpolação | CubicSpline bc_type | bc_type="natural" | presente | PASS |
| interpolação | phi_matrix | spline natural via _natural_splines / delegação | OK | PASS |
| interpolação | _phi_at_point | spline natural via _natural_splines / delegação | OK | PASS |
| interpolação | spectral_targets | spline natural via _natural_splines / delegação | OK | PASS |
| interpolação | spectral_moment_residuals | spline natural via _natural_splines / delegação | OK | PASS |
| interpolação | linear em producao | nenhum | nenhum | PASS |
| numerico | nx campanha | 801 | congelado K4Q-mini (DEFAULT_NX_SPECTRUM=401 em codigo) | PASS |
| numerico | N campanha | 4..14 | congelado K4P (cases.N_VALUES=4..20 no modulo) | PASS |
| casos | Kumaraswamy | Kum(2,8); Kum(2,2.41); Kum(4,1.5) | Kum(2,8); Kum(2,2.41); Kum(4,1.5) | PASS |
| legado | parametros legados ativos | nenhum em conflito com baseline congelada | DEFAULT_NX_SPECTRUM=401 (campanha usa 801 via nx_spectrum); cases.N_VALUES ate 20 (campanha limita a 14); MS_V2_STABLE_QMC presente mas inativo (OK); run_all_cases.py nao passa nx_spectrum=801 (campanha deve faze-lo) | PASS |

## Perguntas

| # | resposta |
| --- | --- |
| Q1 | SIM |
| Q2 | SIM |
| Q3 | SIM |
| Q4 | SIM |
| Q5 | SIM |
| Q6 | SIM |
| Q7 | PARCIAL — constantes legadas no modulo; campanha sobrepoe via spec K4R |
| Q8 | SIM |

### Notas Q7 / Q8

- `DEFAULT_NX_SPECTRUM` no codigo permanece **401**; a campanha usa **801** via `nx_spectrum=801` em `setup_spectral_bundle` (spec K4Q-mini, nao alterada neste gate).
- `cases.N_VALUES` cobre 4..20; a campanha executa apenas **4..14** (spec K4P).
- `MS_V2_STABLE_QMC` existe mas **nao** esta ativo; Stage1 permanece `MS_V1_REFERENCE`.
- Interpolacao **linear** restrita a scripts de auditoria (`tests/audit_k4j_*`, `tests/audit_k4l_*`).

## Autorizacao

Campanha completa: **AUTORIZADA** — baseline congelada neste commit.
