# Dengue exposure dynamics and vaccination scenarios in Brazil

Code for *Reconstructing Dengue Exposure Dynamics and Assessing the Public Health
Impact of Vaccination Scenarios Using Mathematical Models in Brazil*
(Ko, Santos, Huang, Salje, Cummings, Ko).

The analysis has two stages:

1. **Catalytic model (R + Stan).** A Bayesian age-structured catalytic model is
   fitted to SINAN dengue case notifications (27 federative units, 17 age groups,
   2000-2025) to reconstruct state- and age-specific immune profiles.
2. **SIR-type model (MATLAB).** The reconstructed immune profiles and forces of
   infection initialize an age-structured transmission model with vaccination,
   used to project the burden of disease in 2026-2030 under several vaccine
   allocation and supply scenarios. MATLAB also produces all figures and tables.

## Repository layout

```
R/                    catalytic model pipeline (run from the project root)
  00_common.R           shared settings and file names
  01_download_sinan.R   download SINAN line lists          -> data/raw/sinan/
  02_aggregate_cases.R  confirmed cases, year x age x state -> data/processed/
  03_aggregate_population.R  population, year x age x state -> data/processed/
  04_fit_catalytic_model.R   fit the Stan model            -> results_big/catalytic_fit.rds
  05_export_posterior.R      export draws for MATLAB       -> results_big/*.mat
stan/dengue_catalytic.stan   catalytic model
data/processed/       aggregated case and population counts (tracked)
data/raw/             raw downloads (not tracked)
results/              small results shared between steps (tracked)
results_big/          large posterior files (not tracked, see below)
matlab/
  M01 ... M12         one script per group of figures/tables (see mapping below)
  run_all.m           runs M01-M12 in order
  functions/          model and plotting functions; project_config.m holds all paths
  data/               life table, birth rate, vaccine distribution, IBGE state shapefile
  output/             figures and simulation files (not tracked), tables (tracked)
```

## Requirements

The analysis was run on Windows.

- **R** (>= 4.3) with `data.table`, `dplyr`, `tidyr`, `purrr`, `readr`, `posterior`,
  `R.matlab`, `brpop`, and `cmdstanr` with a CmdStan installation (>= 2.33,
  `cmdstanr::install_cmdstan()`).
- **MATLAB** R2024b (the version used for the paper) with the Statistics and
  Machine Learning, Optimization and Mapping toolboxes (Parallel Computing
  Toolbox is optional for M09).
- **Hardware.** The MATLAB scripts run on a standard desktop computer; the
  allocation optimization (M09) takes several hours. The catalytic model fit
  (`R/04`) was run on a 128-core workstation and would take a very long time on
  a personal computer; its results can be downloaded instead (see
  [Large files](#large-files-results_big)).

**Installation.** Install R and MATLAB, then the R packages above and CmdStan
(`cmdstanr::install_cmdstan()`, which compiles CmdStan). This typically takes
15-30 minutes. No further installation is needed; clone or download this
repository and run the scripts as described below.

## Data

**SINAN case notifications.** The yearly dengue line lists (several GB) are not
included. `R/01_download_sinan.R` downloads them from the OpenDataSUS mirror

```
https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SINAN/Dengue/csv/DENGBRyy.csv.zip   (yy = 00..25)
```

into `data/raw/sinan/`. `R/02_aggregate_cases.R` keeps cases confirmed by
laboratory or clinical-epidemiological criteria and writes
`data/processed/dengue_cases_year_age5_state_2000_2025.csv`.

**Population.** `R/03_aggregate_population.R` retrieves municipal population by age
(DataSUS estimates, via the `brpop` package), aggregates it to states and
extrapolates linearly after the last observed year. Output:
`data/processed/pop_state_age5_2000_2025.csv`.

The two processed files are tracked, so the analysis can start at `R/04` without
downloading SINAN.

**MATLAB inputs** (`matlab/data/`): IBGE complete life table (`data_lifeexp.mat`,
`data_survrate.mat`), World Bank crude birth rate 1960-2022 (`data_birthrate.mat`),
dengue vaccine doses distributed by the Ministry of Health in 2024/25
(`vaccine_doses_by_state_2024_2025.csv`) and the IBGE 2023 state boundaries
(`shapefile/BR_UF_2023.*`).

## Large files (`results_big/`)

These files are needed to reproduce the catalytic-model figures and tables
(M01-M05) but are too large for git. They are created by the R pipeline under
exactly these names, and MATLAB reads the same names:

| File | Created by | Content |
|---|---|---|
| `catalytic_fit.rds` | `R/04_fit_catalytic_model.R` | cmdstanr fit object with all draws |
| `catalytic_ppc_draws.mat` | `R/05_export_posterior.R` | immune profiles (`op_seroposall_*`), total infections (`inf_total`), cumulative predicted cases (`cum_cases_age`) |
| `catalytic_param_draws.mat` | `R/05_export_posterior.R` | `lambda_draws`, `lambda_pre_draws`, `kappa_draws`, `omega_draws`, `delta_draws`, `rho_draws` (kappa relative to ages 0-4) |

The SIR-type model only needs the small `results/sir_inputs.mat` written by M05
(posterior medians of the 2026 immune profile, age multipliers and projected
reproduction numbers), so M06-M12 run without `results_big/`.

**Download.** Because of their size, the `results_big/` files used for the
paper are not included in this repository. Fitting the catalytic model is
computationally demanding: we ran it on a 128-core high-performance
workstation, where the run time was manageable, but reproducing the fit on a
personal computer would take a very long time. The three files above can
instead be downloaded from
[this shared folder](https://1drv.ms/f/c/ca3ee35878ee3566/IgBZE_nzOrrDTa8a2EluSa2uASsZPYxNOIs1UlBNLIhjvhE?e=xYVD2A).
Place them in `results_big/` without renaming them; the MATLAB scripts can
then be run directly, skipping `R/04` and `R/05`. If the link is unavailable,
please contact Youngsuk Ko (yko@berkeley.edu).

## How to run

R, from the project root (open `denv_cat_sir_2025-2026.Rproj`):

```r
source("R/01_download_sinan.R")         # optional, data/processed is tracked
source("R/02_aggregate_cases.R")        # optional
source("R/03_aggregate_population.R")   # optional
source("R/04_fit_catalytic_model.R")    # several hours
source("R/05_export_posterior.R")
```

MATLAB, from `matlab/`:

```matlab
run_all        % all figures and tables, reusing results/opt_allocation_*.mat
run_all(true)  % also rerun the allocation optimization (M09, slow)
```

Scripts can also be run one at a time in numeric order. Figures are written to
`matlab/output/figures/`, tables to `matlab/output/tables/`.

The optimized allocations (`results/opt_allocation_shen.mat`,
`results/opt_allocation_daniels.mat`) come from numerical optimization and may
differ slightly between runs; they are tracked so that M10-M12 reproduce the paper.

## Figures and tables

| Paper | Script | Output (`matlab/output/`) |
|---|---|---|
| Fig. 1 (map, age distribution) | M01 | `figures/fig1_*` |
| Fig. 2A-F | M02 | `figures/fig2*` |
| Fig. 3A-F (incl. inset map) | M06 | `figures/fig3*` |
| Fig. 4 | M08 | `figures/fig4_bod_averted_baseline.jpg` |
| Fig. 5 | M11 | `figures/fig5_mixed_strategy.jpg` |
| Table 1 | M12 | `tables/table1_bod_averted_shen.csv` (Daniels version also written) |
| Fig. S1, S2, S4-S8 | M02 | `figures/figS1*`, `figS2*`, `figS4*` ... `figS8*` |
| Fig. S3 | M05 | `figures/figS3_R0_projection.jpg` |
| Fig. S9 | M04 | `figures/figS9_serosurvey_*` |
| Fig. S11, Table S9 | M07 | `figures/figS11*`, `tables/tableS9*` |
| Fig. S12, Table S12 | M06 | `figures/figS12*`, `tables/tableS12*` |
| Fig. S13, Table S10 | M08 | `figures/figS13_*`, `tables/tableS10*` |
| Fig. S14-S18 | M10 | `figures/figS14*` ... `figS18*` |
| Tables S1-S7 | M03 | `tables/tableS1*` ... `tableS7*` |

Fig. 1 (workflow), Fig. S10 (model diagram) and Tables S8 and S11 (parameter
values) were prepared by hand; the parameter values are set in
`matlab/functions/sir_base_parameters.m`. Multi-panel figures are exported panel
by panel and were assembled, with panel letters and region shading of the state
panel maps, outside MATLAB.

## License

The code is released under the MIT License (see `LICENSE`).
