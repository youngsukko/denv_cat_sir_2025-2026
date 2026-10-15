library(tidyverse)
library(cmdstanr)

fit    <- readRDS("results/dengue_catalytic_v4_fixedphi_fit_E.rds")
years  <- 2000:2025; T <- length(years)

dat <- readRDS("data/processed/dengue_year_age5_state_2000_2025.rds")
pop <- readRDS("data/processed/pop_state_age5_2000_2035.rds") %>%
  filter(year %in% years) %>%
  rename(age_group = age_grp) %>%
  select(uf, year, age_group, pop)

states <- sort(unique(dat$uf))
S <- length(states); A5 <- 17

age_labels <- sort(unique(dat$age_group)) %>% as.character()

# age_group 레벨이 dat과 pop에서 동일한지 체크
pop_age_labels <- sort(unique(pop$age_group)) %>% as.character()
if (!setequal(age_labels, pop_age_labels)) {
  warning("dat$age_group 와 pop$age_group 의 값이 다릅니다. 매핑 확인 필요!\n",
          "dat: ", paste(age_labels, collapse=", "), "\n",
          "pop: ", paste(pop_age_labels, collapse=", "))
}
stopifnot(length(age_labels) == A5)

# ── cases_arr / pop_arr 구성 ────────────────────────────────────
full_grid <- expand_grid(uf = states, year = years, age_group = age_labels)

cases_long <- full_grid %>%
  left_join(dat,  by = c("uf", "year", "age_group")) %>%
  mutate(cases = replace_na(cases, 0))

pop_long <- full_grid %>%
  left_join(pop,  by = c("uf", "year", "age_group")) %>%
  mutate(pop = replace_na(pop, 0))

cases_arr <- array(0, dim = c(A5, T, S))
pop_arr   <- array(0, dim = c(A5, T, S))

for (s in seq_len(S)) {
  for (t in seq_len(T)) {
    cs <- cases_long %>%
      filter(uf == states[s], year == years[t]) %>%
      arrange(match(age_group, age_labels))
    ps <- pop_long %>%
      filter(uf == states[s], year == years[t]) %>%
      arrange(match(age_group, age_labels))
    cases_arr[, t, s] <- cs$cases
    pop_arr[, t, s]   <- ps$pop
  }
}

# ── draws 추출 ──────────────────────────────────────────────────
kappa_draws      <- fit$draws("kappa_age",  format = "draws_matrix")
lambda_draws     <- fit$draws("lambda",     format = "draws_matrix")
delta_draws      <- fit$draws("delta",      format = "draws_matrix")
omega_draws      <- fit$draws("omega_age",  format = "draws_matrix")
rho_draws        <- fit$draws("rho",        format = "draws_matrix")
lambda_pre_draws <- fit$draws("lambda_pre", format = "draws_matrix")

# ── kappa 정규화 ─────────────────────────────────────────────────
mean_kappa     <- rowMeans(kappa_draws)
kappa_adj      <- kappa_draws      / mean_kappa
lambda_adj     <- lambda_draws     * mean_kappa
lambda_pre_adj <- lambda_pre_draws * mean_kappa

# ── 컬럼명 ──────────────────────────────────────────────────────
kappa_names      <- paste0("kappa_age_", 1:A5)
lambda_names     <- as.vector(outer(1:S, 1:T, function(s,t)
  paste0("lambda_s", s, "t", t)))
delta_names      <- paste0("delta_", 1:4)
omega_names      <- paste0("omega_age_", 1:A5)
rho_names        <- as.vector(outer(1:S, 1:T, function(s,t)
  paste0("rho_s", s, "t", t)))
lambda_pre_names <- paste0("lambda_pre_s", 1:S)

# ── 순서대로 cbind → CSV ─────────────────────────────────────────
out <- cbind(
  kappa_adj,
  lambda_adj,
  delta_draws,
  omega_draws,
  rho_draws,
  lambda_pre_adj
)
colnames(out) <- c(kappa_names, lambda_names, delta_names,
                   omega_names, rho_names, lambda_pre_names)

dir.create("results", showWarnings = FALSE)
write.csv(out, "results/mcmc_samples.csv", row.names = FALSE)

cat("Saved: results/mcmc_samples.csv\n")
cat("  Rows (draws):", nrow(out), "\n")
cat("  Cols:", ncol(out), "\n")
cat("  kappa_age:   cols 1 ~", A5, "\n")
cat("  lambda:      cols", A5+1, "~", A5+S*T, "\n")
cat("  delta:       cols", A5+S*T+1, "~", A5+S*T+4, "\n")
cat("  omega_age:   cols", A5+S*T+5, "~", A5+S*T+4+A5, "\n")
cat("  rho:         cols", A5+S*T+5+A5, "~", A5+S*T+4+A5+S*T, "\n")
cat("  lambda_pre:  cols", A5+S*T+5+A5+S*T, "~", ncol(out), "\n")

# ── cases_arr / pop_arr → CSV (행 = (s-1)*T+t, 열 = 연령밴드) ─────
cases_export <- matrix(0, S*T, A5)
pop_export   <- matrix(0, S*T, A5)
for (s in 1:S) for (t in 1:T) {
  cases_export[(s-1)*T+t, ] <- cases_arr[, t, s]
  pop_export  [(s-1)*T+t, ] <- pop_arr  [, t, s]
}
colnames(cases_export) <- age_labels
colnames(pop_export)   <- age_labels

write.csv(cases_export, "results/cases_age_state_year.csv", row.names = FALSE)
write.csv(pop_export,   "results/pop_age_state_year.csv",   row.names = FALSE)

# ── incidence per 100 ───────────────────────────────────────────
incidence_export <- ifelse(pop_export > 0,
                           cases_export / pop_export * 100,
                           NA_real_)
colnames(incidence_export) <- age_labels
write.csv(incidence_export, "results/incidence_per100.csv", row.names = FALSE)

# ── state 순서 매핑 테이블 (MATLAB reg_order 대조용) ──────────────
write.csv(tibble(idx = seq_along(states), uf = states),
          "results/state_order.csv", row.names = FALSE)

cat("Done\n")