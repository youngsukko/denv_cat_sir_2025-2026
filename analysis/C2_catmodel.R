# ══════════════════════════════════════════════════════════════════
# Brazil DENV Catalytic Model — v4 fixed phi + full PPC extraction
# - phi[s,b]: 주×연령 moment estimator로 고정값 입력
# - kappa_age: truncated normal RW (lower=0)
# - omega_age: truncated normal RW ([0,1])
# - generated quantities: cases_pred (PPC), seropos_0~4 (post-2000),
#   seropos_pre_0~4 (pre-2000 baseline)
# ══════════════════════════════════════════════════════════════════

library(tidyverse)
library(cmdstanr)
library(R.matlab)

#% ── 0. 기본 설정 ────────────────────────────────────────────────
YEAR_MIN <- 2000
YEAR_MAX <- 2025
years    <- YEAR_MIN:YEAR_MAX
T        <- length(years)

age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                "35-39","40-44","45-49","50-54","55-59","60-64",
                "65-69","70-74","75-79","80+")
A5 <- length(age_labels)
A1 <- 81
band_of_age <- c(rep(1:16, each = 5), 17)
stopifnot(length(band_of_age) == A1)

#% ── 1. 데이터 로드 ───────────────────────────────────────────────
case_summary <- readRDS("data/processed/dengue_year_age5_state_2000_2025.rds")
pop_full     <- readRDS("data/processed/pop_state_age5_2000_2035.rds")

states <- sort(unique(case_summary$uf))
S      <- length(states)
stopifnot(S == 27, "ES" %in% states)
cat("States (", S, "):", paste(states, collapse = ", "), "\n")

cases_arr <- array(0L, dim = c(A5, T, S))
for (k in seq_len(nrow(case_summary))) {
  b <- match(as.character(case_summary$age_group[k]), age_labels)
  t <- match(case_summary$year[k], years)
  s <- match(case_summary$uf[k], states)
  if (!anyNA(c(b, t, s))) cases_arr[b, t, s] <- case_summary$cases[k]
}

pop_filt <- pop_full %>% filter(year >= YEAR_MIN, year <= YEAR_MAX)
pop_arr  <- array(1.0, dim = c(A5, T, S))
for (k in seq_len(nrow(pop_filt))) {
  b <- match(as.character(pop_filt$age_grp[k]), age_labels)
  t <- match(pop_filt$year[k], years)
  s <- match(pop_filt$uf[k], states)
  if (!anyNA(c(b, t, s))) pop_arr[b, t, s] <- max(pop_filt$pop[k], 1.0)
}

cat("Total cases", years[1], "-", years[T], ":", sum(cases_arr), "\n")
cat("Total pop in", years[T], ":", sum(pop_arr[, T, ]), "\n\n")

#% ── 2. 케이스 수 가중치: sqrt 기반, cap 없음 ──────────────────────
state_cases <- apply(cases_arr, 3, sum)
pop_weight  <- mean(sqrt(state_cases)) / sqrt(state_cases)
cat("Case weights (mean=1):\n")
print(round(sort(pop_weight), 3))

#% ── 3. phi[s,b]: 주×연령 moment estimator ─────────────────────────
phi_matrix <- case_summary %>%
  group_by(uf, age_group) %>%
  summarise(
    mean_cases = mean(cases),
    var_cases  = var(cases),
    .groups = "drop"
  ) %>%
  mutate(
    phi_approx = ifelse(
      var_cases > mean_cases,
      mean_cases^2 / (var_cases - mean_cases),
      0.05
    ),
    age_group = factor(age_group, levels = age_labels),
    uf        = factor(uf, levels = states)
  ) %>%
  arrange(uf, age_group) %>%
  pivot_wider(
    id_cols     = uf,
    names_from  = age_group,
    values_from = phi_approx
  ) %>%
  select(-uf) %>%
  as.matrix()  # S×A5

cat("\nphi_matrix (주×연령, 일부):\n")
print(round(phi_matrix[1:5, 1:5], 3))

stan_data <- list(
  S           = S,
  A5          = A5,
  A1          = A1,
  T           = T,
  band_of_age = as.array(band_of_age),
  cases       = cases_arr,
  pop         = pop_arr,
  pop_weight  = as.array(pop_weight),
  phi_fixed   = phi_matrix   # S×A5 고정값
)

#% ── 4. Stan 코드 ─────────────────────────────────────────────────
stan_code <- "functions {
  real partial_sum(array[] int slice_s,
                   int start, int end,
                   int T, int A1, int A5,
                   array[] int band_of_age,
                   array[,,] int cases,
                   array[,,] real pop,
                   matrix lambda,
                   vector lambda_pre,
                   matrix rho,
                   array[,] real phi_fixed,
                   vector kappa_full,
                   vector omega_age,
                   vector delta,
                   vector pop_weight) {
    real target_acc = 0.0;
    for (idx in 1:size(slice_s)) {
      int s = slice_s[idx];

      array[A1, T + 1] real Lam;
      Lam[1, 1] = 0.0;
      for (i in 2:A1)
        Lam[i, 1] = Lam[i-1, 1] + lambda_pre[s] * kappa_full[i-1];
      for (t in 1:T) {
        Lam[1, t+1] = kappa_full[1] * lambda[s, t];
        for (i in 2:(A1-1))
          Lam[i, t+1] = Lam[i-1, t] + kappa_full[i] * lambda[s, t];
        Lam[A1, t+1] = Lam[A1, t] + kappa_full[A1] * lambda[s, t];
      }

      for (t in 1:T) {
        array[A5] real band_sum;
        array[A5] int  band_n;
        for (b in 1:A5) { band_sum[b] = 0.0; band_n[b] = 0; }
        vector[A1] dh = kappa_full * lambda[s, t];

        for (i in 1:A1) {
          real Lb;
          if (i == 1)       Lb = 0.0;
          else if (i == A1) Lb = Lam[A1, t];
          else              Lb = Lam[i-1, t];

          real p = -expm1(-Lb);
          real q = 1 - p;
          real p2 = p*p; real p3 = p2*p;
          real q2 = q*q; real q3 = q2*q; real q4 = q3*q;

          real app =
            delta[1] *      q4       * (-expm1(-4 * dh[i])) +
            delta[2] * (4*p *q3)     * (-expm1(-3 * dh[i])) +
            delta[3] * (6*p2*q2)     * (-expm1(-2 * dh[i])) +
            delta[4] * (4*p3*q)      * (-expm1(-1 * dh[i]));

          int b = band_of_age[i];
          band_sum[b] += app;
          band_n[b]   += 1;
        }

        for (b in 1:A5) {
          real mu = fmax(rho[s,t] * omega_age[b] * (band_sum[b]/band_n[b]) * pop[b,t,s], 1e-8);
          target_acc += pop_weight[s] * neg_binomial_2_lpmf(cases[b,t,s] | mu, phi_fixed[s,b]);
        }
      }
    }
    return target_acc;
  }
}

data {
  int<lower=1> S;
  int<lower=1> A5;
  int<lower=1> A1;
  int<lower=1> T;
  array[A1] int<lower=1, upper=A5>  band_of_age;
  array[A5, T, S] int<lower=0>      cases;
  array[A5, T, S] real<lower=0>     pop;
  vector<lower=0>[S]                pop_weight;
  array[S, A5] real<lower=0>        phi_fixed;  // 주×연령 고정값
}

transformed data {
  array[S] int seq_S;
  for (s in 1:S) seq_S[s] = s;
}

parameters {
  matrix<lower=0>[S, T]   lambda;
  vector<lower=0>[S]      lambda_pre;

  vector<lower=0>[A5]              kappa_age;
  vector<lower=0, upper=1>[A5]     omega_age;

  real<lower=0, upper=1> delta_primary;
  real<lower=0, upper=1> delta_secondary;
  real<lower=0, upper=1> delta_tert_frac;
  real<lower=0, upper=1> delta_quat_frac;

  matrix<lower=0, upper=1>[S, T]   rho;
}

transformed parameters {
  vector[A1] kappa_full;
  for (i in 1:A1) kappa_full[i] = kappa_age[band_of_age[i]];

  vector[4] delta;
  delta[1] = delta_primary;
  delta[2] = delta_secondary;
  delta[3] = delta_primary * delta_tert_frac;
  delta[4] = delta[3]      * delta_quat_frac;
}

model {
  to_vector(lambda) ~ exponential(4/mean(kappa_age));
  lambda_pre        ~ exponential(40/mean(kappa_age));
  to_vector(rho) ~ beta(1, 1.5);

  // kappa_age ~ normal(1, 4);
  target += S * normal_lpdf(kappa_age[1]    | 1, 1);
  for (i in 2:A5) target += S * normal_lpdf(kappa_age[i]    | kappa_age[i-1], 2);
  
  target += S * beta_lpdf(omega_age[1]    | 1.5, 1);
  for (i in 2:A5) target += S * normal_lpdf(omega_age[i]    | omega_age[i-1], 2);

  target += S * beta_lpdf(delta_primary   | 4, 16);
  target += S * beta_lpdf(delta_secondary | 2, 2);

  target += reduce_sum(partial_sum, seq_S, 1,
                       T, A1, A5, band_of_age, cases, pop,
                       lambda, lambda_pre, rho, phi_fixed,
                       kappa_full, omega_age, delta,
                       pop_weight);
}

generated quantities {
  array[A5, T, S] int   cases_pred;
  array[A5, T, S] real  seropos_0;
  array[A5, T, S] real  seropos_1;
  array[A5, T, S] real  seropos_2;
  array[A5, T, S] real  seropos_3;
  array[A5, T, S] real  seropos_4;
  array[A5, S]    real  seropos_pre_0;
  array[A5, S]    real  seropos_pre_1;
  array[A5, S]    real  seropos_pre_2;
  array[A5, S]    real  seropos_pre_3;
  array[A5, S]    real  seropos_pre_4;

  for (s in 1:S) {
    array[A1, T+1] real Lam;
    Lam[1,1] = 0.0;
    for (i in 2:A1) Lam[i,1] = Lam[i-1,1] + lambda_pre[s] * kappa_full[i-1];
    for (t in 1:T) {
      Lam[1,t+1] = kappa_full[1] * lambda[s,t];
      for (i in 2:(A1-1)) Lam[i,t+1] = Lam[i-1,t] + kappa_full[i] * lambda[s,t];
      Lam[A1,t+1] = Lam[A1,t] + kappa_full[A1] * lambda[s,t];
    }

    // ── pre-2000 baseline (t=0 시점, Lam[,1]) ──────────────────
    {
      array[A5] real bx0; array[A5] real bx1; array[A5] real bx2;
      array[A5] real bx3; array[A5] real bx4; array[A5] int nx;
      for (b in 1:A5) { bx0[b]=0.0; bx1[b]=0.0; bx2[b]=0.0; bx3[b]=0.0; bx4[b]=0.0; nx[b]=0; }

      for (i in 1:A1) {
        real pe = -expm1(-Lam[i,1]); real qe = 1-pe;
        real pe2=pe*pe; real pe3=pe2*pe;
        real qe2=qe*qe; real qe3=qe2*qe; real qe4=qe3*qe;

        int b = band_of_age[i];
        bx0[b] += qe4;
        bx1[b] += 4*pe*qe3;
        bx2[b] += 6*pe2*qe2;
        bx3[b] += 4*pe3*qe;
        bx4[b] += pe2*pe2;
        nx[b]  += 1;
      }
      for (b in 1:A5) {
        seropos_pre_0[b,s] = bx0[b]/nx[b];
        seropos_pre_1[b,s] = bx1[b]/nx[b];
        seropos_pre_2[b,s] = bx2[b]/nx[b];
        seropos_pre_3[b,s] = bx3[b]/nx[b];
        seropos_pre_4[b,s] = bx4[b]/nx[b];
      }
    }

    // ── 연도별 (post-2000): cases_pred + seropos_0~4 ────────────
    for (t in 1:T) {
      array[A5] real bapp; array[A5] int napp;
      array[A5] real bx0; array[A5] real bx1; array[A5] real bx2;
      array[A5] real bx3; array[A5] real bx4; array[A5] int nx;
      for (b in 1:A5) {
        bapp[b]=0.0; napp[b]=0;
        bx0[b]=0.0; bx1[b]=0.0; bx2[b]=0.0; bx3[b]=0.0; bx4[b]=0.0; nx[b]=0;
      }
      vector[A1] dh = kappa_full * lambda[s,t];

      for (i in 1:A1) {
        real Lb;
        if (i==1)       Lb = 0.0;
        else if (i==A1) Lb = Lam[A1,t];
        else            Lb = Lam[i-1,t];

        real p = -expm1(-Lb); real q = 1-p;
        real p2=p*p; real p3=p2*p;
        real q2=q*q; real q3=q2*q; real q4=q3*q;

        real app =
          delta[1]*q4       *(-expm1(-4*dh[i])) +
          delta[2]*(4*p*q3) *(-expm1(-3*dh[i])) +
          delta[3]*(6*p2*q2)*(-expm1(-2*dh[i])) +
          delta[4]*(4*p3*q) *(-expm1(-1*dh[i]));

        int b = band_of_age[i];
        bapp[b] += app; napp[b] += 1;

        real pe = -expm1(-Lam[i,t+1]); real qe=1-pe;
        real pe2=pe*pe; real pe3=pe2*pe;
        real qe2=qe*qe; real qe3=qe2*qe; real qe4=qe3*qe;

        bx0[b] += qe4;
        bx1[b] += 4*pe*qe3;
        bx2[b] += 6*pe2*qe2;
        bx3[b] += 4*pe3*qe;
        bx4[b] += pe2*pe2;
        nx[b]  += 1;
      }

      for (b in 1:A5) {
        real mu = fmax(rho[s,t]*omega_age[b]*(bapp[b]/napp[b])*pop[b,t,s], 1e-8);
        cases_pred[b,t,s] = neg_binomial_2_rng(mu, phi_fixed[s,b]);
        seropos_0[b,t,s] = bx0[b]/nx[b];
        seropos_1[b,t,s] = bx1[b]/nx[b];
        seropos_2[b,t,s] = bx2[b]/nx[b];
        seropos_3[b,t,s] = bx3[b]/nx[b];
        seropos_4[b,t,s] = bx4[b]/nx[b];
      }
    }
  }
}
"

dir.create("models", showWarnings = FALSE)
writeLines(stan_code, "models/dengue_catalytic_v4_fixedphi_ppc.stan")

mod <- cmdstan_model(
  "models/dengue_catalytic_v4_fixedphi_ppc.stan",
  force_recompile = TRUE,
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

dir.create("results/temp_samples", recursive = TRUE, showWarnings = FALSE)

#% ── 5. 샘플링 ────────────────────────────────────────────────────
fit <- mod$sample(
  data              = stan_data,
  seed              = 42,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = 27,
  iter_warmup       = 500,
  iter_sampling     = 500,
  refresh           = 20,
  init              = init_fn,
  max_treedepth     = 12,
  output_dir        = "results/temp_samples"
)

fit$save_object("results/fit_fincand1_ppc.rds")

cat("\nDiagnostics:\n")
print(fit$diagnostic_summary())

cat("\nkappa_age:\n")
print(fit$summary("kappa_age") %>% mutate(age = age_labels))

cat("\nomega_age:\n")
print(fit$summary("omega_age") %>% mutate(age = age_labels))

cat("\ndelta:\n")
print(fit$summary("delta"))

cat("\nlambda_pre:\n")
print(n = 27, fit$summary("lambda_pre") %>% mutate(state = states))

cat("\nmean rho by state:\n")
print(fit$summary("rho") %>%
        mutate(state = rep(states, each = T), year = rep(years, S)) %>%
        group_by(state) %>%
        summarise(mean_rho = mean(mean)))

# ══════════════════════════════════════════════════════════════════
#  PPC / Immune profile draws 추출 + 저장 (.rds + .mat)
# ══════════════════════════════════════════════════════════════════

n_draw <- fit$draws(variables = "cases_pred[1,1,1]", format = "draws_matrix") %>% nrow()
cat("\nn_draw =", n_draw, "\n")

# Stan array[A5,T,S] -> draws_matrix 컬럼 순서는 "var[1,1,1]","var[2,1,1]",...
# 즉 b(A5)가 가장 빠르게 변함(column-major) → array(dim=c(n_draw,A5,T,S))로
# 그대로 채우면 Stan 인덱스 순서(b,t,s)와 일치
extract_array <- function(fit, varname, A5, T, S, n_draw) {
  m <- fit$draws(variables = varname, format = "draws_matrix")  # [n_draw x A5*T*S]
  stopifnot(ncol(m) == A5*T*S)
  array(as.numeric(m), dim = c(n_draw, A5, T, S))
}

extract_array_2d <- function(fit, varname, A5, S, n_draw) {
  m <- fit$draws(variables = varname, format = "draws_matrix")  # [n_draw x A5*S]
  stopifnot(ncol(m) == A5*S)
  array(as.numeric(m), dim = c(n_draw, A5, S))
}

cat("Extracting cases_pred...\n")
cases_pred <- extract_array(fit, "cases_pred", A5, T, S, n_draw)   # [n_draw,A5,T,S]

cat("Extracting seropos_0..4 (post-2000)...\n")
seropos_0 <- extract_array(fit, "seropos_0", A5, T, S, n_draw)
seropos_1 <- extract_array(fit, "seropos_1", A5, T, S, n_draw)
seropos_2 <- extract_array(fit, "seropos_2", A5, T, S, n_draw)
seropos_3 <- extract_array(fit, "seropos_3", A5, T, S, n_draw)
seropos_4 <- extract_array(fit, "seropos_4", A5, T, S, n_draw)

cat("Extracting seropos_pre_0..4 (pre-2000 baseline)...\n")
seropos_pre_0 <- extract_array_2d(fit, "seropos_pre_0", A5, S, n_draw)
seropos_pre_1 <- extract_array_2d(fit, "seropos_pre_1", A5, S, n_draw)
seropos_pre_2 <- extract_array_2d(fit, "seropos_pre_2", A5, S, n_draw)
seropos_pre_3 <- extract_array_2d(fit, "seropos_pre_3", A5, S, n_draw)
seropos_pre_4 <- extract_array_2d(fit, "seropos_pre_4", A5, S, n_draw)

#% ── 누적(2000-2025) 케이스 연령분포 (모델 PPC 기반) [n_draw x S x A5] ──
cum_cases_age <- array(0, dim = c(n_draw, S, A5))
for (s in 1:S) {
  for (b in 1:A5) {
    cum_cases_age[, s, b] <- rowSums(cases_pred[, b, , s])  # T에 대해 합산
  }
}
cum_case_age_prop_model <- array(0, dim = c(n_draw, S, A5))
for (s in 1:S) {
  rs <- rowSums(cum_cases_age[, s, ])
  cum_case_age_prop_model[, s, ] <- cum_cases_age[, s, ] / rs
}

#% ── 주별 연도별 incidence (PPC 기반, 연령 합산) [n_draw x T x S] ──
pop_bystate_byyear <- apply(pop_arr, c(2,3), sum)   # [T x S]

model_inc <- array(0, dim = c(n_draw, T, S))
for (s in 1:S) {
  cases_total_t <- apply(cases_pred[, , , s], c(1,3), sum)  # [n_draw x T]
  for (t in 1:T) {
    model_inc[, t, s] <- cases_total_t[, t] / pop_bystate_byyear[t, s] * 100
  }
}

#% ── 저장 ────────────────────────────────────────────────────────
dir.create("results", showWarnings = FALSE)

saveRDS(list(
  cases_pred              = cases_pred,
  seropos_0 = seropos_0, seropos_1 = seropos_1, seropos_2 = seropos_2,
  seropos_3 = seropos_3, seropos_4 = seropos_4,
  seropos_pre_0 = seropos_pre_0, seropos_pre_1 = seropos_pre_1,
  seropos_pre_2 = seropos_pre_2, seropos_pre_3 = seropos_pre_3,
  seropos_pre_4 = seropos_pre_4,
  cum_cases_age           = cum_cases_age,
  cum_case_age_prop_model = cum_case_age_prop_model,
  model_inc               = model_inc,
  states = states, years = years, age_labels = age_labels
), "results/catalytic_ppc_draws.rds", compress = "xz")

# MATLAB용: [n_draw, A5, T, S] -> [n_draw, S, T, A5] (기존 MATLAB 변수 구조와 호환)
to_matlab_order <- function(arr_b_t_s) {
  aperm(arr_b_t_s, c(1,4,3,2))   # [n_draw,A5,T,S] -> [n_draw,S,T,A5]
}
to_matlab_order_2d <- function(arr_b_s) {
  aperm(arr_b_s, c(1,3,2))       # [n_draw,A5,S] -> [n_draw,S,A5]
}

writeMat("results/catalytic_ppc_draws.mat",
         op_seroposall_0 = to_matlab_order(seropos_0),
         op_seroposall_1 = to_matlab_order(seropos_1),
         op_seroposall_2 = to_matlab_order(seropos_2),
         op_seroposall_3 = to_matlab_order(seropos_3),
         op_seroposall_4 = to_matlab_order(seropos_4),
         op_seroposall_pre_0 = to_matlab_order_2d(seropos_pre_0),
         op_seroposall_pre_1 = to_matlab_order_2d(seropos_pre_1),
         op_seroposall_pre_2 = to_matlab_order_2d(seropos_pre_2),
         op_seroposall_pre_3 = to_matlab_order_2d(seropos_pre_3),
         op_seroposall_pre_4 = to_matlab_order_2d(seropos_pre_4),
         model_inc               = model_inc,                 # [n_draw x T x S]
         cum_cases_age           = cum_cases_age,              # [n_draw x S x A5]
         cum_case_age_prop_model = cum_case_age_prop_model     # [n_draw x S x A5]
)

cat("\nSaved: results/catalytic_ppc_draws.rds (R backup, full draws)\n")
cat("Saved: results/catalytic_ppc_draws.mat (for MATLAB plotting)\n")
cat("Done.\n")


S <- 27; T <- 26; A5 <- 17

lambda_m <- fit$draws(variables = "lambda", format = "draws_matrix")  # [n_draw x S*T], s가 빠르게 변함
n_draw   <- nrow(lambda_m)
lambda_draws <- array(as.numeric(lambda_m), dim = c(n_draw, S, T))    # [n_draw,S,T]

kappa_draws <- as.matrix(fit$draws(variables = "kappa_age", format = "draws_matrix"))  # [n_draw, A5]
omega_draws <- as.matrix(fit$draws(variables = "omega_age", format = "draws_matrix"))  # [n_draw, A5]
delta_draws <- as.matrix(fit$draws(variables = "delta",     format = "draws_matrix"))  # [n_draw, 4]

rho_m <- fit$draws(variables = "rho", format = "draws_matrix")   # [n_draw x S*T], s가 빠르게 변함
rho_draws <- array(as.numeric(rho_m), dim = c(n_draw, S, T))     # [n_draw,S,T]

lambda_pre_draws <- as.matrix(fit$draws(variables = "lambda_pre", format = "draws_matrix"))  # [n_draw,S]

saveRDS(list(
  lambda_draws = lambda_draws, kappa_draws = kappa_draws,
  omega_draws = omega_draws, delta_draws = delta_draws,
  rho_draws = rho_draws, lambda_pre_draws = lambda_pre_draws,
  n_draw = n_draw, S = S, T = T, A5 = A5
), "results/catalytic_lambda_kappa_draws.rds", compress = "xz")

writeMat("results/catalytic_lambda_kappa_draws.mat",
         lambda_draws = lambda_draws,        # [n_draw x S x T]
         kappa_draws  = kappa_draws,         # [n_draw x A5]
         omega_draws  = omega_draws,         # [n_draw x A5]
         delta_draws  = delta_draws,         # [n_draw x 4]
         rho_draws    = rho_draws,           # [n_draw x S x T]
         lambda_pre_draws = lambda_pre_draws # [n_draw x S]
)

cat("Saved: results/catalytic_lambda_kappa_draws.mat (with omega/delta/rho/lambda_pre)\n")