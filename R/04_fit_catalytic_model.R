# Fit the age-structured catalytic model (stan/dengue_catalytic.stan) to the
# SINAN case counts, 27 states x 26 years x 17 age groups.
#
# Input : data/processed/dengue_cases_year_age5_state_2000_2025.csv
#         data/processed/pop_state_age5_2000_2025.csv
# Output: results_big/catalytic_fit.rds   (cmdstanr fit object, all draws)
#
# Runtime is long (hours on a workstation). Afterwards run
# R/05_export_posterior.R to create the files read by MATLAB.

source("R/00_common.R")

library(dplyr)
library(readr)
library(cmdstanr)

# -- Sampler settings (Methods) -------------------------------------
N_CHAINS          <- 4
ITER_WARMUP       <- 500
ITER_SAMPLING     <- 500
MAX_TREEDEPTH     <- 12
SEED              <- 42
THREADS_PER_CHAIN <- 27   # reduce_sum parallelizes the likelihood over states

years <- YEAR_MIN:YEAR_MAX
T     <- length(years)
A5    <- length(AGE_LABELS)
A1    <- 81                                   # single-year cohorts 0..80+
band_of_age <- c(rep(1:16, each = 5), 17)
stopifnot(length(band_of_age) == A1)

# -- 1. Data --------------------------------------------------------
case_summary <- read_csv(FILE_CASES, col_types = cols(uf = "c", year = "i",
                                                      age_group = "c", cases = "i"))
pop_state    <- read_csv(FILE_POP, col_types = cols(uf = "c", year = "i",
                                                    age_group = "c", pop = "d"))

states <- sort(unique(case_summary$uf))
S      <- length(states)
stopifnot(S == 27)

idx_cases <- cbind(match(case_summary$age_group, AGE_LABELS),
                   match(case_summary$year, years),
                   match(case_summary$uf, states))
stopifnot(!anyNA(idx_cases))
cases_arr <- array(0L, dim = c(A5, T, S))
cases_arr[idx_cases] <- case_summary$cases

pop_filt <- pop_state |> filter(year >= YEAR_MIN, year <= YEAR_MAX)
idx_pop  <- cbind(match(pop_filt$age_group, AGE_LABELS),
                  match(pop_filt$year, years),
                  match(pop_filt$uf, states))
stopifnot(!anyNA(idx_pop))
pop_arr <- array(1.0, dim = c(A5, T, S))
pop_arr[idx_pop] <- pmax(pop_filt$pop, 1.0)

cat("Total cases", YEAR_MIN, "-", YEAR_MAX, ":", sum(cases_arr), "\n")

# -- 2. State weights: inverse square root of total cases (mean 1) --
state_cases <- apply(cases_arr, 3, sum)
pop_weight  <- mean(sqrt(state_cases)) / sqrt(state_cases)

# -- 3. Overdispersion phi_s, method of moments on annual total cases --
annual_totals <- apply(cases_arr, c(2, 3), sum)          # [T x S]
mean_annual   <- colMeans(annual_totals)
var_annual    <- apply(annual_totals, 2, var)
if (any(var_annual <= mean_annual)) {
  stop("Annual case counts are not overdispersed in: ",
       paste(states[var_annual <= mean_annual], collapse = ", "))
}
phi_state <- mean_annual^2 / (var_annual - mean_annual)

print(data.frame(state = states, weight = round(pop_weight, 3), phi = round(phi_state, 3)))

stan_data <- list(
  S           = S,
  A5          = A5,
  A1          = A1,
  T           = T,
  band_of_age = as.array(band_of_age),
  cases       = cases_arr,
  pop         = pop_arr,
  pop_weight  = as.array(pop_weight),
  phi_fixed   = as.array(phi_state)
)

# -- 4. Compile and sample ------------------------------------------
mod <- cmdstan_model(
  FILE_STAN_MODEL,
  cpp_options = list(stan_threads = TRUE,
                     CXXFLAGS = "-O3 -march=native -mtune=native")
)

init_fn <- function() list(
  lambda          = matrix(0.02, S, T),
  lambda_pre      = rep(0.005, S),
  kappa_age       = rep(1.0, A5),
  omega_age       = rep(0.5, A5),
  delta_primary   = 0.05,
  delta_secondary = 0.5,
  delta_tert_frac = 0.5,
  delta_quat_frac = 0.5,
  rho             = matrix(0.1, S, T)
)

output_dir <- file.path(DIR_BIG, "cmdstan_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fit <- mod$sample(
  data              = stan_data,
  seed              = SEED,
  chains            = N_CHAINS,
  parallel_chains   = N_CHAINS,
  threads_per_chain = THREADS_PER_CHAIN,
  iter_warmup       = ITER_WARMUP,
  iter_sampling     = ITER_SAMPLING,
  refresh           = 20,
  init              = init_fn,
  max_treedepth     = MAX_TREEDEPTH,
  output_dir        = output_dir
)

fit$save_object(FILE_FIT)
cat("Saved:", FILE_FIT, "\n")

# -- 5. Quick look --------------------------------------------------
print(fit$diagnostic_summary())
print(fit$summary("kappa_age") |> mutate(age = AGE_LABELS))
print(fit$summary("omega_age") |> mutate(age = AGE_LABELS))
print(fit$summary("delta"))
