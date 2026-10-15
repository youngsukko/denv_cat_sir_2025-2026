# ════════════════════════════════════════════════════════════════════
# SINAN Dengue — QC visuals, per 100,000 population
# ════════════════════════════════════════════════════════════════════
# Plot A: 27-panel facet (one per state), ANNUAL INCIDENCE by year
#         (cases summed over age groups / population summed over age
#         groups, matched year-by-year).
# Plot B: nationwide age distribution, AVERAGE ANNUAL INCIDENCE per
#         age group (cases summed 2000-2025 / population summed over
#         the same years -- i.e. average annual rate over the period,
#         not a single-year snapshot).
# ════════════════════════════════════════════════════════════════════

library(tidyverse)

case_summary <- read_csv("data/processed/dengue_year_age5_state_2000_2025.csv")
pop_state    <- read_csv("data/processed/pop_state_age5_2000_2035.csv") %>%
  rename(age_group = age_grp)

# ── Region / state ordering (same convention as the CHIKV C1 scripts) ──
regiao_map <- tibble(
  uf     = c("RO","AC","AM","RR","PA","AP","TO",
             "MA","PI","CE","RN","PB","PE","AL","SE","BA",
             "MG","ES","RJ","SP",
             "PR","SC","RS",
             "MS","MT","GO","DF"),
  regiao = c(rep("North",        7),
             rep("Northeast",    9),
             rep("Southeast",    4),
             rep("South",        3),
             rep("Central-West", 4))
)

regiao_order  <- c("North", "Northeast", "Central-West", "Southeast", "South")
regiao_colors <- c(
  "North"        = "#1565C0",
  "Northeast"    = "#E65100",
  "Central-West" = "#6A1B9A",
  "Southeast"    = "#AD1457",
  "South"        = "#2E7D32"
)

uf_order <- regiao_map %>%
  mutate(regiao = factor(regiao, levels = regiao_order)) %>%
  arrange(regiao, uf) %>%
  pull(uf)

age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                "35-39","40-44","45-49","50-54","55-59","60-64",
                "65-69","70-74","75-79","80+")

case_summary <- case_summary %>%
  left_join(regiao_map, by = "uf") %>%
  mutate(
    uf        = factor(uf, levels = uf_order),
    age_group = factor(age_group, levels = age_labels),
    regiao    = factor(regiao, levels = regiao_order)
  )

pop_state <- pop_state %>%
  mutate(
    uf        = factor(uf, levels = uf_order),
    age_group = factor(age_group, levels = age_labels)
  )

# ── Join cases + population ────────────────────────────────────────
# YEAR_MIN:YEAR_MAX from the case script only goes to 2025, well within
# pop_state's observed+extrapolated 2000-2035 range, so every case row
# should find a population match. If pop is NA anywhere below, that
# uf/year/age_group combo fell outside pop_state's coverage.
case_pop <- case_summary %>%
  left_join(
    pop_state %>% select(uf, year, age_group, pop),
    by = c("uf", "year", "age_group")
  )

n_missing_pop <- sum(is.na(case_pop$pop))
if (n_missing_pop > 0) {
  message(n_missing_pop, " case rows have no matching population row -- ",
          "incidence will be NA for those.")
}

# ══════════════════════════════════════════════════════════════════
# Plot A: annual incidence by year, faceted by state
# (cases and population both summed across age groups first, THEN
#  divided -- not an average of per-age-group incidences)
# ══════════════════════════════════════════════════════════════════

inc_year_state <- case_pop %>%
  group_by(uf, regiao, year) %>%
  summarise(cases = sum(cases), pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
  mutate(incidence = if_else(pop > 0, cases / pop * 1e5, NA_real_))

p_year_state <- ggplot(inc_year_state, aes(x = year, y = incidence, color = regiao)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ uf, scales = "free_y", ncol = 5) +
  scale_color_manual(values = regiao_colors, guide = "none") +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_y_continuous(labels = scales::comma_format(accuracy = 1)) +
  labs(x = NULL, y = "Cases per 100,000") +
  theme_bw(base_size = 10) +
  theme(
    strip.text   = element_text(face = "bold", size = 8),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 6)
  )

print(p_year_state)

# ══════════════════════════════════════════════════════════════════
# Plot B: nationwide age distribution, average annual incidence
# (cases and population summed across ALL years 2000-2025 first,
#  so the ratio is an average annual rate over the period)
# ══════════════════════════════════════════════════════════════════

inc_age_nat <- case_pop %>%
  group_by(age_group) %>%
  summarise(cases = sum(cases), pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
  mutate(incidence = cases / pop * 1e5)

p_age_nat <- ggplot(inc_age_nat, aes(x = age_group, y = incidence)) +
  geom_col(fill = "steelblue", width = 0.7) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Age group", y = "Avg. annual cases per 100,000 (2000-2025, nationwide)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(p_age_nat)

# ── Save ──────────────────────────────────────────────────────────
dir.create("figures", showWarnings = FALSE)
ggsave("figures/dengue_incidence_year_by_state.png", p_year_state, width = 12, height = 9, dpi = 400)
ggsave("figures/dengue_incidence_age_nationwide.png", p_age_nat,    width = 6,  height = 4, dpi = 400)