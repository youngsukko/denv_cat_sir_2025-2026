# ══════════════════════════════════════════════════════════════════
# Brazil DENV — FOI panel map (kappa-normalized)
# ══════════════════════════════════════════════════════════════════

library(tidyverse)
library(cmdstanr)
library(patchwork)

YEAR_MIN <- 2000; YEAR_MAX <- 2025
years    <- YEAR_MIN:YEAR_MAX; T <- length(years)
A5       <- 17

# ── 1. 로드 ──────────────────────────────────────────────────────

# fit          <- readRDS("results/dengue_catalytic_v2_fit_plain_default.rds")
fit          <- readRDS("results/dengue_catalytic_v4_fixedphi_fit.rds")
case_summary <- readRDS("data/processed/dengue_year_age5_state_2000_2025.rds")
pop_full     <- readRDS("data/processed/pop_state_age5_2000_2035.rds")

states <- sort(unique(case_summary$uf))
S      <- length(states)
stopifnot(S == 27)

# 풀네임 룩업
state_names <- c(
  "AC"="Acre",            "AL"="Alagoas",           "AP"="Amapá",
  "AM"="Amazonas",        "BA"="Bahia",              "CE"="Ceará",
  "DF"="Distrito Federal","ES"="Espírito Santo",     "GO"="Goiás",
  "MA"="Maranhão",        "MT"="Mato Grosso",        "MS"="Mato Grosso do Sul",
  "MG"="Minas Gerais",    "PA"="Pará",               "PB"="Paraíba",
  "PR"="Paraná",          "PE"="Pernambuco",         "PI"="Piauí",
  "RJ"="Rio de Janeiro",  "RN"="Rio Grande do Norte","RS"="Rio Grande do Sul",
  "RO"="Rondônia",        "RR"="Roraima",            "SC"="Santa Catarina",
  "SP"="São Paulo",       "SE"="Sergipe",            "TO"="Tocantins"
)

# 인시던스
obs_year <- case_summary %>%
  group_by(uf, year) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

pop_year <- pop_full %>%
  filter(year >= YEAR_MIN, year <= YEAR_MAX) %>%
  group_by(uf, year) %>%
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

incidence <- obs_year %>%
  left_join(pop_year, by = c("uf","year")) %>%
  mutate(inc_per100k = cases / pop * 1e2)

# 인시던스 전국 최대값으로 통일
inc_max_global <- max(incidence$inc_per100k, na.rm = TRUE) * 1.05

# ── 2. kappa 정규화 ──────────────────────────────────────────────

kappa_draws      <- fit$draws("kappa_age",  format = "draws_matrix")
lambda_draws     <- fit$draws("lambda",     format = "draws_matrix")

mean_kappa   <- rowMeans(kappa_draws)
lambda_adj   <- lambda_draws * mean_kappa

lambda_adj_s <- tibble(
  variable = colnames(lambda_draws),
  q5  = apply(lambda_adj, 2, quantile, 0.05),
  q50 = apply(lambda_adj, 2, quantile, 0.50),
  q95 = apply(lambda_adj, 2, quantile, 0.95)
) %>%
  mutate(
    s     = as.integer(str_extract(variable, "(?<=\\[)\\d+")),
    t     = as.integer(str_extract(variable, "(?<=,)\\d+(?=\\])")),
    state = states[s],
    year  = years[t]
  )

# ── 3. 패널 레이아웃 ─────────────────────────────────────────────

panel_layout <- tribble(
  ~state, ~row, ~col,
  "RR",      1,    3,   # 1→3
  "AP",      1,    5,   # 3→5
  "AM",      2,    3,   # 1→3
  "PA",      2,    5,   # 3→5
  "AC",      3,    1,
  "RO",      3,    3,
  "TO",      3,    5,
  "MA",      2,    7,   # 5→7
  "CE",      2,    9,   # 7→9
  "PI",      3,    7,
  "PB",      3,    9,
  "RN",      3,   11,
  "MT",      4,    3,
  "GO",      4,    5,
  "BA",      4,    7,
  "PE",      4,    9,
  "AL",      4,   11,
  "MS",      5,    3,
  "DF",      5,    5,
  "MG",      5,    7,
  "SE",      5,    9,
  "ES",      6,    9,
  "SP",      6,    5,
  "RJ",      6,    7,
  "PR",      7,    5,
  "SC",      7,    7,
  "RS",      8,    5
) # panel_layout tribble 끝

max_col_idx <- max(panel_layout$col)   # 11 (빈 열 포함)
max_row_idx <- max(panel_layout$row)   # 8

# ── 4. 패널 생성 함수 ────────────────────────────────────────────

FOI_MAX <- 0.5
FOI_COL <- "#C0392B"

make_panel <- function(st) {
  foi     <- lambda_adj_s %>% filter(state == st)
  inc     <- incidence    %>% filter(uf    == st)
  inc_max <- inc_max_global
  scale_f <- FOI_MAX / inc_max
  
  ggplot() +
    geom_ribbon(data = foi,
                aes(x = year, ymin = q5, ymax = q95),
                fill = "#E8998D", alpha = 0.65) +
    geom_line(data = foi,
              aes(x = year, y = q50),
              color = FOI_COL, linewidth = 0.55) +
    geom_line(data = inc,
              aes(x = year, y = inc_per100k * scale_f),
              linetype = "dotted", color = "black", linewidth = 0.45) +
    scale_y_continuous(
      limits   = c(0, FOI_MAX),
      breaks   = c(0, 0.15, 0.30, 0.45),
      labels   = c("0","0.15","0.30","0.45"),
      sec.axis = sec_axis(
        transform = ~. / scale_f,
        breaks    = pretty(c(0, inc_max), n = 3)
      )
    ) +
    scale_x_continuous(
      limits = c(YEAR_MIN, YEAR_MAX),
      breaks = c(2000, 2010, 2020),
      labels = c("2000","2010","2020")
    ) +
    ggtitle(state_names[st]) +
    theme_bw(base_size = 7) +
    theme(
      panel.grid         = element_blank(),
      axis.title         = element_blank(),
      axis.text.x        = element_text(size = 7),
      axis.text.y.left   = element_text(color = FOI_COL, size = 7),
      axis.ticks.y.left  = element_line(color = FOI_COL),
      axis.text.y.right  = element_text(color = "black", size = 7),
      axis.ticks.y.right = element_line(color = "black"),
      plot.title         = element_text(size = 7, face = "bold",
                                        hjust = 0.5, margin = margin(b=1)),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5)   # 기존 margin(2,2,2,2) → 0.5로 축소
    )
}

panels <- setNames(lapply(panel_layout$state, make_panel), panel_layout$state)

# ── 5. 조립 ──────────────────────────────────────────────────────

blank <- ggplot() + theme_void()

grid_plots <- vector("list", max_row_idx * max_col_idx)
for (r in seq_len(max_row_idx)) {
  for (c in seq_len(max_col_idx)) {
    idx <- (r - 1) * max_col_idx + c
    st  <- panel_layout %>% filter(row == r, col == c) %>% pull(state)
    grid_plots[[idx]] <- if (length(st) == 1) panels[[st]] else blank
  }
}

fig <- wrap_plots(grid_plots, nrow = max_row_idx, ncol = max_col_idx) +
  plot_layout(widths = ifelse(seq_len(max_col_idx) %% 2 == 1, 1, 0.15)) &
  theme(plot.margin = margin(0.0, 0.0, 0.0, 0.0))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/foi_panel_map.png", fig,
       width = 7, height = 8, dpi = 500)
ggsave("figures/foi_panel_map.pdf", fig,
       width = 7, height = 8)

cat("Saved: figures/foi_panel_map.png\n")