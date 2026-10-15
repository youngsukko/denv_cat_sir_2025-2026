# ══════════════════════════════════════════════════════════════════
# Brazil DENV Catalytic Model — results visualization (v2 fit)
#
# Run AFTER C2_dengue_catalytic_stan.R has finished and saved
# results/dengue_catalytic_v2_fit.rds. This script doesn't refit
# anything -- it reloads the saved fit + the raw case data (for the
# "observed" overlays) and produces:
#
#   1. ppc_year_by_state.png  -- 27-panel, age-summed cases by year
#                                 (observed points vs model ribbon)
#   2. ppc_age_by_state.png   -- 27-panel, cumulative cases by age band
#                                 (observed bars vs model pointrange)
#   3. kappa_age.png          -- national age-FOI multiplier
#   4. omega_age.png          -- national age-reporting multiplier
#   5. lambda_by_state.png    -- 27-panel FOI lambda(s,t), pre-2000
#                                 lambda_pre shown as dashed reference
#   6. rho_by_state.png       -- 27-panel reporting ratio rho(s,t)
#   7. delta.png              -- infection-order apparent-illness bars
#
# Note: mu_cases isn't saved in the v2 (reduce_sum) model -- it's
# local inside partial_sum -- so PPC plots use cases_pred (noisy
# posterior-predictive draws) instead, which is the standard way to
# do this anyway.
# ══════════════════════════════════════════════════════════════════

library(tidyverse)
library(cmdstanr)

# ══════════════════════════════════════════════════════════════════
# 0. Reconstruct the deterministic setup (must match the fitting script)
# ══════════════════════════════════════════════════════════════════

YEAR_MIN <- 2000
YEAR_MAX <- 2025
years    <- YEAR_MIN:YEAR_MAX
T        <- length(years)

age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                "35-39","40-44","45-49","50-54","55-59","60-64",
                "65-69","70-74","75-79","80+")
A5 <- length(age_labels)

delta_labels <- c("Primary", "Secondary", "Tertiary", "Quaternary")

# ══════════════════════════════════════════════════════════════════
# 1. Reload raw case data, just for the "observed" overlays
# ══════════════════════════════════════════════════════════════════

case_summary <- readRDS("data/processed/dengue_year_age5_state_2000_2025.rds")
states <- sort(unique(case_summary$uf))
S      <- length(states)
stopifnot(S == 27)

cases_arr <- array(0L, dim = c(A5, T, S))
for (k in seq_len(nrow(case_summary))) {
  b <- match(as.character(case_summary$age_group[k]), age_labels)
  t <- match(case_summary$year[k], years)
  s <- match(case_summary$uf[k], states)
  if (!anyNA(c(b, t, s))) cases_arr[b, t, s] <- case_summary$cases[k]
}

# ══════════════════════════════════════════════════════════════════
# 2. Load the fit
# ══════════════════════════════════════════════════════════════════

# fit <- readRDS("results/dengue_catalytic_v2_fit.rds")
fit <- readRDS("results/dengue_catalytic_v2_fit_plain_default.rds")

# parse_idx(): splits "var[1,2,3]"-style names into columns i1, i2, i3
parse_idx <- function(df) {
  idx_str   <- str_match(df$variable, "\\[([0-9,]+)\\]")[, 2]
  idx_split <- str_split(idx_str, ",")
  max_n     <- max(lengths(idx_split))
  idx_mat   <- do.call(rbind, lapply(idx_split, function(x) {
    x <- as.integer(x)
    length(x) <- max_n
    x
  }))
  colnames(idx_mat) <- paste0("i", seq_len(max_n))
  bind_cols(df, as_tibble(idx_mat))
}

dir.create("figures", showWarnings = FALSE)
theme_min <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

# ══════════════════════════════════════════════════════════════════
# 3. Case PPC by year, per state (age-summed)
# ══════════════════════════════════════════════════════════════════

cases_pred_s <- fit$summary("cases_pred") %>%
  parse_idx() %>%
  rename(b = i1, t = i2, s = i3) %>%
  mutate(year = years[t], state = states[s],
         age_grp = factor(age_labels[b], levels = age_labels))

ppc_year <- cases_pred_s %>%
  group_by(state, year) %>%
  summarise(pred_mean = sum(mean), pred_q5 = sum(q5), pred_q95 = sum(q95), .groups = "drop")

obs_year <- map_dfr(seq_len(S), function(s) {
  tibble(state = states[s], year = years, obs = apply(cases_arr[, , s], 2, sum))
})

p_ppc_year <- ggplot() +
  geom_ribbon(data = ppc_year, aes(x = year, ymin = pred_q5, ymax = pred_q95),
              fill = "#4575B4", alpha = 0.25) +
  geom_line(data = ppc_year, aes(x = year, y = pred_mean), color = "#4575B4", linewidth = 0.5) +
  geom_point(data = obs_year, aes(x = year, y = obs), color = "black", size = 0.6) +
  facet_wrap(~state, scales = "free_y", ncol = 5) +
  labs(title = "주별 연도별 케이스 (점 = 관측, 리본/선 = 모델)", x = NULL, y = "Cases") +
  theme_min +
  theme(strip.text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 6))
print(p_ppc_year)
ggsave("figures/ppc_year_by_state.png", p_ppc_year, width = 16, height = 18, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 4. Case PPC by age band, per state (cumulative across all years)
# ══════════════════════════════════════════════════════════════════

ppc_age <- cases_pred_s %>%
  group_by(state, age_grp) %>%
  summarise(pred_mean = sum(mean), pred_q5 = sum(q5), pred_q95 = sum(q95), .groups = "drop")

obs_age <- map_dfr(seq_len(S), function(s) {
  tibble(state = states[s], age_grp = factor(age_labels, levels = age_labels),
         obs = apply(cases_arr[, , s], 1, sum))
})

p_ppc_age <- ggplot() +
  geom_col(data = obs_age, aes(x = age_grp, y = obs), fill = "#4575B4", alpha = 0.35) +
  geom_pointrange(data = ppc_age,
                  aes(x = age_grp, y = pred_mean, ymin = pred_q5, ymax = pred_q95),
                  color = "#08306B", size = 0.2, linewidth = 0.4) +
  facet_wrap(~state, scales = "free_y", ncol = 5) +
  labs(title = "주별 연령대별 누적 케이스 (막대 = 관측, 점 = 모델)", x = NULL, y = "Cumulative cases") +
  theme_min +
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5))
print(p_ppc_age)
ggsave("figures/ppc_age_by_state.png", p_ppc_age, width = 18, height = 18, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 5. National age multipliers: kappa_age, omega_age
# ══════════════════════════════════════════════════════════════════

kappa_s <- fit$summary("kappa_age") %>% parse_idx() %>%
  mutate(age_grp = factor(age_labels[i1], levels = age_labels))
omega_s <- fit$summary("omega_age") %>% parse_idx() %>%
  mutate(age_grp = factor(age_labels[i1], levels = age_labels))

p_kappa <- ggplot(kappa_s, aes(x = age_grp, y = mean, ymin = q5, ymax = q95)) +
  geom_pointrange(color = "#3182BD") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(title = "연령별 FOI 보정인자 (kappa_age)", x = NULL, y = "kappa_age") +
  theme_min +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_kappa)
ggsave("figures/kappa_age.png", p_kappa, width = 7, height = 4, dpi = 300)

p_omega <- ggplot(omega_s, aes(x = age_grp, y = mean, ymin = q5, ymax = q95)) +
  geom_pointrange(color = "#E6550D") +
  labs(title = "연령별 보고율(현성감염+의료이용) 보정인자 (omega_age)", x = NULL, y = "omega_age") +
  theme_min +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_omega)
ggsave("figures/omega_age.png", p_omega, width = 7, height = 4, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 6. FOI by state: lambda(s,t), pre-2000 lambda_pre as reference
# ══════════════════════════════════════════════════════════════════

lambda_s <- fit$summary("lambda") %>% parse_idx() %>%
  rename(s = i1, t = i2) %>%
  mutate(year = years[t], state = states[s])

lambda_pre_s <- fit$summary("lambda_pre") %>% parse_idx() %>%
  rename(s = i1) %>%
  mutate(state = states[s])

p_lambda <- ggplot(lambda_s, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#D73027", alpha = 0.25) +
  geom_line(color = "#D73027", linewidth = 0.5) +
  geom_hline(data = lambda_pre_s, aes(yintercept = mean), linetype = "dashed", color = "grey30") +
  facet_wrap(~state, scales = "free_y", ncol = 5) +
  labs(title = "주별 FOI (lambda; 점선 = pre-2000 lambda_pre)", x = NULL, y = "lambda") +
  theme_min +
  theme(strip.text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 6))
print(p_lambda)
ggsave("figures/lambda_by_state.png", p_lambda, width = 16, height = 18, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 7. Reporting ratio by state: rho(s,t)
# ══════════════════════════════════════════════════════════════════

rho_s <- fit$summary("rho") %>% parse_idx() %>%
  rename(s = i1, t = i2) %>%
  mutate(year = years[t], state = states[s])

p_rho <- ggplot(rho_s, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#4DAC26", alpha = 0.25) +
  geom_line(color = "#4DAC26", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~state, ncol = 5) +
  labs(title = "주별 신고율 (rho)", x = NULL, y = "rho") +
  theme_min +
  theme(strip.text = element_text(size = 7), axis.text.x = element_text(angle = 90, size = 6))
print(p_rho)
ggsave("figures/rho_by_state.png", p_rho, width = 16, height = 18, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 8. Infection-order apparent-illness ratio (delta) -- bar chart
# ══════════════════════════════════════════════════════════════════

delta_s <- fit$summary("delta") %>% parse_idx() %>%
  mutate(order = factor(delta_labels[i1], levels = delta_labels))

p_delta <- ggplot(delta_s, aes(x = order, y = mean)) +
  geom_col(fill = "#756BB1", width = 0.6) +
  geom_errorbar(aes(ymin = q5, ymax = q95), width = 0.2, color = "#3F007D") +
  labs(title = "감염 순서별 증상발현율 (delta: 1~4차)", x = NULL, y = "delta") +
  theme_min
print(p_delta)
ggsave("figures/delta.png", p_delta, width = 5, height = 4, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 9. Immune profile by state: 2000 & 2025
#
# pi_exposed[b, t, s] = 1 - mean_i_in_band(q_end^4)  where q = 1-p.
# We need p (per-serotype cumulative infection prob) to reconstruct the
# full Binomial(4, p) distribution. Since pi_exposed = 1 - (1-p)^4
# under the band-average approximation:
#   p_approx = 1 - (1 - pi_exposed)^(1/4)
#
# This approximation is exact for a homogeneous band and accurate to
# within a few pp for realistic within-band heterogeneity. We use the
# posterior mean of pi_exposed for the point estimates (the CI would
# require propagating through the nonlinear transform, which inflates
# slightly -- omitted here for readability).
#
# Two panels per state:
#   left  = 2000 baseline  (t = 1,  pi_exposed[b, 1, s])
#   right = 2025 end       (t = T,  pi_exposed[b, T, s])
# ══════════════════════════════════════════════════════════════════

pi_s <- fit$summary("pi_exposed") %>%
  parse_idx() %>%
  rename(b = i1, t = i2, s = i3) %>%
  mutate(age_grp = factor(age_labels[b], levels = age_labels),
         state   = states[s],
         year    = years[t])

# helper: given a tibble with one row per age band (pi_mean col),
# return long tibble with Binomial(4,p) probabilities for k=0..4
expand_binom4 <- function(df_pi) {
  k_levels <- c("0회 (미감염)", "1차", "2차", "3차", "4차")
  map_dfr(seq_len(nrow(df_pi)), function(i) {
    p  <- 1 - (1 - pmin(pmax(df_pi$mean[i], 0), 1))^0.25
    q  <- 1 - p
    pk <- c(
      q^4,
      4   * p   * q^3,
      6   * p^2 * q^2,
      4   * p^3 * q,
      p^4
    )
    tibble(
      age_grp  = df_pi$age_grp[i],
      state    = df_pi$state[i],
      k        = factor(k_levels, levels = k_levels),
      prob     = pk
    )
  })
}

immune_2000 <- pi_s %>% filter(t == 1)  %>% expand_binom4()
immune_2025 <- pi_s %>% filter(t == T)  %>% expand_binom4()

k_colors <- c(
  "0회 (미감염)" = "#CCCCCC",
  "1차"          = "#74ADD1",
  "2차"          = "#4575B4",
  "3차"          = "#F46D43",
  "4차"          = "#A50026"
)

make_immune_plot <- function(df, year_label) {
  ggplot(df, aes(x = age_grp, y = prob, color = k, group = k)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 0.8) +
    scale_color_manual(values = k_colors, name = "감염 횟수") +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    facet_wrap(~state, ncol = 5) +
    labs(title = paste0("주별 연령별 면역 프로파일 — ", year_label,
                        " (Binomial(4,p) 분포)"),
         x = NULL, y = "비율") +
    theme_min +
    theme(strip.text      = element_text(size = 7),
          axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          legend.position = "bottom")
}

p_immune_2000 <- make_immune_plot(immune_2000, "2000년 기준")
p_immune_2025 <- make_immune_plot(immune_2025, "2025년 기준")

print(p_immune_2000)
ggsave("figures/immune_profile_2000.png", p_immune_2000, width = 18, height = 18, dpi = 300)

print(p_immune_2025)
ggsave("figures/immune_profile_2025.png", p_immune_2025, width = 18, height = 18, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 10. Population-weighted seroprevalence by state: 2000 vs 2025
#
# 인구수로 가중평균한 주별 항체양성률 비교.
# 정확한 posterior 분포를 얻으려면 draws가 필요하지만,
# 각 밴드의 분위수를 인구수로 가중합산하는 근사를 사용.
# (정확한 분위수는 아니지만 시각화 목적으로 합리적인 근사)
# ══════════════════════════════════════════════════════════════════

# 인구 데이터 로드 (가중치용)
pop_full_viz <- readRDS("data/processed/pop_state_age5_2000_2035.rds")
pop_arr_viz  <- array(1.0, dim = c(A5, T, S))
pop_filt_viz <- pop_full_viz %>% filter(year >= YEAR_MIN, year <= YEAR_MAX)
for (k in seq_len(nrow(pop_filt_viz))) {
  b <- match(as.character(pop_filt_viz$age_grp[k]), age_labels)
  t <- match(pop_filt_viz$year[k], years)
  s <- match(pop_filt_viz$uf[k], states)
  if (!anyNA(c(b, t, s))) pop_arr_viz[b, t, s] <- max(pop_filt_viz$pop[k], 1.0)
}

# pi_exposed 분위수 추출 (t=1, t=T만)
pi_q <- fit$summary(
  "pi_exposed",
  mean = mean,
  ~quantile(.x, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
) %>%
  rename(q5 = `5%`, q25 = `25%`, q50 = `50%`, q75 = `75%`, q95 = `95%`) %>%
  parse_idx() %>%
  rename(b = i1, t = i2, s = i3) %>%
  filter(t %in% c(1, T)) %>%
  mutate(
    year_label = factor(ifelse(t == 1, "2000년", "2025년"), levels = c("2000년","2025년")),
    state      = states[s],
    pop_wt     = map2_dbl(b, t, ~pop_arr_viz[.x, .y, s[cur_group_id()]])
  )

# 인구 가중합산 (밴드별 분위수의 인구 가중 근사)
sero_box <- pi_q %>%
  rowwise() %>%
  mutate(pop_wt = pop_arr_viz[b, t, s]) %>%
  ungroup() %>%
  group_by(state, year_label) %>%
  summarise(
    pop_total = sum(pop_wt),
    w_mean = sum(mean * pop_wt) / pop_total,
    w_q5   = sum(q5   * pop_wt) / pop_total,
    w_q25  = sum(q25  * pop_wt) / pop_total,
    w_q50  = sum(q50  * pop_wt) / pop_total,
    w_q75  = sum(q75  * pop_wt) / pop_total,
    w_q95  = sum(q95  * pop_wt) / pop_total,
    .groups = "drop"
  ) %>%
  mutate(state = factor(state, levels = sort(unique(state))))

year_colors <- c("2000년" = "#74ADD1", "2025년" = "#D73027")

p_sero_state <- ggplot(sero_box,
                       aes(x = state, fill = year_label, color = year_label,
                           ymin = w_q5, lower = w_q25, middle = w_q50,
                           upper = w_q75, ymax = w_q95)) +
  geom_boxplot(stat = "identity", position = position_dodge(width = 0.7),
               width = 0.6, alpha = 0.7) +
  scale_fill_manual(values  = year_colors, name = NULL) +
  scale_color_manual(values = year_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(title = "주별 전체인구 가중 항체양성률 (2000년 vs 2025년)",
       x = NULL, y = "항체양성률") +
  theme_min +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "top")

print(p_sero_state)
ggsave("figures/seroprevalence_by_state.png", p_sero_state, width = 16, height = 6, dpi = 300)

cat("\nDone. Figures saved under figures/:\n")
cat("  ppc_year_by_state.png, ppc_age_by_state.png,\n")
cat("  kappa_age.png, omega_age.png,\n")
cat("  lambda_by_state.png, rho_by_state.png,\n")
cat("  delta.png\n")
cat("  immune_profile_2000.png, immune_profile_2025.png\n")
cat("  seroprevalence_by_state.png\n")