# ════════════════════════════════════════════════════════════════════
# SINAN Dengue — QC visuals from the year x age5 x state case summary
# ════════════════════════════════════════════════════════════════════
# Plot A: 27-panel facet (one per state), cases by year, summed over
#         age groups.
# Plot B: nationwide age distribution, summed over all years/states.
# ════════════════════════════════════════════════════════════════════

library(tidyverse)

case_summary <- read_csv("data/processed/dengue_year_age5_state_2000_2025.csv")

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

# ══════════════════════════════════════════════════════════════════
# Plot A: cases by year, faceted by state (summed across age groups)
# ══════════════════════════════════════════════════════════════════

cases_year_state <- case_summary %>%
  group_by(uf, regiao, year) %>%
  summarise(cases = sum(cases), .groups = "drop")

p_year_state <- ggplot(cases_year_state, aes(x = year, y = cases, color = regiao)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ uf, scales = "free_y", ncol = 5) +
  scale_color_manual(values = regiao_colors, guide = "none") +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5)) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = NULL, y = "Cases") +
  theme_bw(base_size = 10) +
  theme(
    strip.text   = element_text(face = "bold", size = 8),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 6)
  )

print(p_year_state)

# ══════════════════════════════════════════════════════════════════
# Plot B: nationwide age distribution (summed across all years/states)
# ══════════════════════════════════════════════════════════════════

cases_age_nat <- case_summary %>%
  group_by(age_group) %>%
  summarise(cases = sum(cases), .groups = "drop")

p_age_nat <- ggplot(cases_age_nat, aes(x = age_group, y = cases)) +
  geom_col(fill = "steelblue", width = 0.7) +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Age group", y = "Total cases (2000-2025, nationwide)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(p_age_nat)

# ── Save ──────────────────────────────────────────────────────────
dir.create("figures", showWarnings = FALSE)
ggsave("figures/dengue_year_by_state.png", p_year_state, width = 12, height = 9, dpi = 400)
ggsave("figures/dengue_age_nationwide.png", p_age_nat,    width = 6,  height = 4, dpi = 400)