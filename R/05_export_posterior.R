# Export the posterior of the catalytic model for MATLAB and write the MCMC
# convergence diagnostics reported in Supplementary S3-1.
#
# Input : results_big/catalytic_fit.rds          (R/04_fit_catalytic_model.R)
# Output: results_big/catalytic_ppc_draws.mat    immune profiles, predicted cases,
#                                                total infections
#         results_big/catalytic_param_draws.mat  model parameters
#         results/mcmc_sampler_diagnostics.csv   divergences, treedepth, E-BFMI
#         results/mcmc_convergence_diagnostics.csv  R-hat, bulk/tail ESS

source("R/00_common.R")

library(dplyr)
library(readr)
library(cmdstanr)
library(posterior)
library(R.matlab)

fit <- readRDS(FILE_FIT)

years  <- YEAR_MIN:YEAR_MAX
T      <- length(years)
A5     <- length(AGE_LABELS)
states <- sort(unique(read_csv(FILE_CASES, col_types = cols(uf = "c"))$uf))
S      <- length(states)

# -- 1. MCMC diagnostics (S3-1) -------------------------------------
dir.create(DIR_RESULTS, showWarnings = FALSE)

sampler_diag <- as.data.frame(fit$diagnostic_summary()) |>
  mutate(chain = seq_len(n()), .before = 1)
write_csv(sampler_diag, file.path(DIR_RESULTS, "mcmc_sampler_diagnostics.csv"))
print(sampler_diag)

param_vars <- c("lambda", "lambda_pre", "kappa_age", "omega_age", "delta", "rho")
conv_diag  <- fit$summary(variables = param_vars, "rhat", "ess_bulk", "ess_tail")
write_csv(conv_diag, file.path(DIR_RESULTS, "mcmc_convergence_diagnostics.csv"))

cat("\nmax R-hat:", round(max(conv_diag$rhat), 3),
    " min bulk ESS:", round(min(conv_diag$ess_bulk)),
    " min tail ESS:", round(min(conv_diag$ess_tail)), "\n")
cat("Parameters with R-hat > 1.05 or ESS < 400:\n")
print(conv_diag |> filter(rhat > 1.05 | ess_bulk < 400 | ess_tail < 400), n = Inf)

# -- 2. Draw extraction helpers -------------------------------------
# For array[A5, T, S] Stan variables the draws_matrix columns run with the
# first index fastest, so filling array(dim = c(n_draw, A5, T, S)) keeps the
# Stan index order.
n_draw <- nrow(fit$draws(variables = "delta[1]", format = "draws_matrix"))

extract_3d <- function(varname) {
  m <- fit$draws(variables = varname, format = "draws_matrix")
  stopifnot(ncol(m) == A5 * T * S)
  array(as.numeric(m), dim = c(n_draw, A5, T, S))
}
extract_2d <- function(varname) {
  m <- fit$draws(variables = varname, format = "draws_matrix")
  stopifnot(ncol(m) == A5 * S)
  array(as.numeric(m), dim = c(n_draw, A5, S))
}

# MATLAB layout: [n_draw, A5, T, S] -> [n_draw, S, T, A5] and [n_draw, A5, S] -> [n_draw, S, A5]
to_matlab_3d <- function(x) aperm(x, c(1, 4, 3, 2))
to_matlab_2d <- function(x) aperm(x, c(1, 3, 2))

# -- 3. Immune profiles, predicted cases, total infections ---------
cases_pred <- extract_3d("cases_pred")

cum_cases_age <- array(0, dim = c(n_draw, S, A5))   # cumulative 2000-2025 predicted cases
for (s in 1:S) {
  for (b in 1:A5) {
    cum_cases_age[, s, b] <- rowSums(cases_pred[, b, , s])
  }
}
rm(cases_pred)

writeMat(FILE_PPC_MAT,
         op_seroposall_0     = to_matlab_3d(extract_3d("seropos_0")),     # [n_draw x S x T x A5]
         op_seroposall_1     = to_matlab_3d(extract_3d("seropos_1")),
         op_seroposall_2     = to_matlab_3d(extract_3d("seropos_2")),
         op_seroposall_3     = to_matlab_3d(extract_3d("seropos_3")),
         op_seroposall_4     = to_matlab_3d(extract_3d("seropos_4")),
         op_seroposall_pre_0 = to_matlab_2d(extract_2d("seropos_pre_0")), # [n_draw x S x A5]
         op_seroposall_pre_1 = to_matlab_2d(extract_2d("seropos_pre_1")),
         op_seroposall_pre_2 = to_matlab_2d(extract_2d("seropos_pre_2")),
         op_seroposall_pre_3 = to_matlab_2d(extract_2d("seropos_pre_3")),
         op_seroposall_pre_4 = to_matlab_2d(extract_2d("seropos_pre_4")),
         inf_total           = to_matlab_3d(extract_3d("inf_total")),     # [n_draw x S x T x A5]
         cum_cases_age       = cum_cases_age)                             # [n_draw x S x A5]
cat("Saved:", FILE_PPC_MAT, "\n")

# -- 4. Parameters, normalized to the 0-4 age group ----------------
# Only the products kappa_a * lambda enter the likelihood, so kappa is
# reported relative to kappa[1] (0-4 years) and lambda, lambda_pre are
# rescaled by kappa[1] to keep the products unchanged.
lambda_draws     <- array(as.numeric(fit$draws(variables = "lambda", format = "draws_matrix")),
                          dim = c(n_draw, S, T))
rho_draws        <- array(as.numeric(fit$draws(variables = "rho", format = "draws_matrix")),
                          dim = c(n_draw, S, T))
kappa_draws      <- as.matrix(fit$draws(variables = "kappa_age",  format = "draws_matrix"))
omega_draws      <- as.matrix(fit$draws(variables = "omega_age",  format = "draws_matrix"))
delta_draws      <- as.matrix(fit$draws(variables = "delta",      format = "draws_matrix"))
lambda_pre_draws <- as.matrix(fit$draws(variables = "lambda_pre", format = "draws_matrix"))

kappa1           <- kappa_draws[, 1]
kappa_draws      <- sweep(kappa_draws, 1, kappa1, FUN = "/")
lambda_pre_draws <- sweep(lambda_pre_draws, 1, kappa1, FUN = "*")
lambda_draws     <- array(as.numeric(lambda_draws) * rep(kappa1, times = S * T),
                          dim = c(n_draw, S, T))

writeMat(FILE_PARAM_MAT,
         lambda_draws     = lambda_draws,      # [n_draw x S x T]
         kappa_draws      = kappa_draws,       # [n_draw x A5], first column == 1
         omega_draws      = omega_draws,       # [n_draw x A5]
         delta_draws      = delta_draws,       # [n_draw x 4]
         rho_draws        = rho_draws,         # [n_draw x S x T]
         lambda_pre_draws = lambda_pre_draws)  # [n_draw x S]
cat("Saved:", FILE_PARAM_MAT, "\n")
