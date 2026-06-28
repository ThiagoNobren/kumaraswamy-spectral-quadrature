# K4S — Campaign Freeze

**Data:** 2026-06-25 21:50  
**Branch:** `feature/natural-cubic-eigenfunctions`  
**Commit:** `c8fde528eba95fff1a39998d28711cd20810932d`  
**CSV:** `outputs/csv/k4s_campaign_freeze.csv`

## Infraestrutura congelada

| item | valor |
| --- | --- |
| Entrypoint | `campaign/run_chapter_campaign.py` |
| Config | `campaign/chapter_campaign_config.py` |
| nx | **801** (explicito) |
| N | **4..14** (explicito) |
| Output CSV | `outputs/csv/chapter_campaign_multi_results.csv` |

`run_all_cases.py` permanece legado (defaults 401 / N=4..20); **nao** usar para o capitulo.

## Verificacoes

| categoria | item | esperado | observado | status |
| --- | --- | --- | --- | --- |
| reprodutibilidade | commit | (hash) | c8fde528eba9 | PASS |
| reprodutibilidade | branch | feature/natural-cubic-eigenfunctions | feature/natural-cubic-eigenfunctions | PASS |
| antes | run_all_cases dependia defaults | sim (legado) | sim | PASS |
| config | NX_SPECTRUM | 801 | 801 | PASS |
| config | N_VALUES | 4..14 | [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14] | PASS |
| campanha | nx=801 explicito | setup_spectral_bundle(..., nx_spectrum=...) | sim | PASS |
| campanha | N=4..14 explicito | CAMPAIGN_N_VALUES de chapter_campaign_config | sim | PASS |
| campanha | DEFAULT_NX_SPECTRUM no runner | ausente | ausente | PASS |
| campanha | cases.N_VALUES no runner | ausente | ausente | PASS |
| legado | run_all_cases ainda usa defaults | sim (nao usado pela campanha) | N_VALUES default=True; nx default=True | PASS |
| baseline | nucleo K4R inalterado | sem diff em solver/stage1/spline | sem diff | PASS |
| campanha | entrypoint oficial | campaign/run_chapter_campaign.py | presente | PASS |

## Perguntas

| # | resposta |
| --- | --- |
| Q1 | SIM |
| Q2 | SIM |
| Q3 | SIM |
| Q4 | NAO (campanha); SIM apenas em run_all_cases legado |
| Q5 | NAO (campanha); SIM apenas em run_all_cases legado |
| Q6 | SIM |
| Q7 | SIM |
| Q8 | SIM |

## Estado

Campanha: **CONGELADA E PRONTA**.

Sem merge em `main`. Sem simulacoes nesta auditoria.
