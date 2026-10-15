# ════════════════════════════════════════════════════════════════════
# State x 5-year Age Group population, 2000-2035 (extrapolated)
# ════════════════════════════════════════════════════════════════════
# Aggregates municipality-level population (brpop::mun_pop_age) up to
# state level, then linearly extrapolates beyond the last year of real
# data through 2035 -- same two-point slope method used in
# C1_case_age_region.R (there: 2019->2021 extrapolated to 2022-2025),
# just generalized to whatever the actual last observed year is.
# ════════════════════════════════════════════════════════════════════

library(tidyverse)

YEAR_MIN        <- 2000
YEAR_TARGET_MAX <- 2035

# ── 1. Load raw municipality-level population ─────────────────────
# library(brpop)
# pop_all_raw <- mun_pop_age(source = "datasus", sex = "all")
# saveRDS(pop_all_raw, "data/pop/pop_all_raw_datasus.rds")
pop_all_raw <- readRDS("data/pop/pop_all_raw_datasus.rds")

# ── 2. UF mapping (IBGE 2-digit -> abbreviation), same as the case script ──
estado_name_map <- tibble(
  co_estado = c("11","12","13","14","15","16","17",
                "21","22","23","24","25","26","27","28","29",
                "31","32","33","35",
                "41","42","43",
                "50","51","52","53"),
  uf        = c("RO","AC","AM","RR","PA","AP","TO",
                "MA","PI","CE","RN","PB","PE","AL","SE","BA",
                "MG","ES","RJ","SP",
                "PR","SC","RS",
                "MS","MT","GO","DF")
)

# ── 3. Age group labels -> match the 5-year bins used in the case data ──
age_label_map <- c(
  "From 0 to 4 years"     = "0-4",
  "From 5 to 9 years"     = "5-9",
  "From 10 to 14 years"   = "10-14",
  "From 15 to 19 years"   = "15-19",
  "From 20 to 24 years"   = "20-24",
  "From 25 to 29 years"   = "25-29",
  "From 30 to 34 years"   = "30-34",
  "From 35 to 39 years"   = "35-39",
  "From 40 to 44 years"   = "40-44",
  "From 45 to 49 years"   = "45-49",
  "From 50 to 54 years"   = "50-54",
  "From 55 to 59 years"   = "55-59",
  "From 60 to 64 years"   = "60-64",
  "From 65 to 69 years"   = "65-69",
  "From 70 to 74 years"   = "70-74",
  "From 75 to 79 years"   = "75-79",
  "From 80 years or more" = "80+"
)
age_labels <- unname(age_label_map)

# ── 4. Aggregate municipality pop -> state x year x age group ──────
pop_state_year_age <- pop_all_raw %>%
  filter(age_group != "Total") %>%
  mutate(
    co_municipio = as.integer(code_muni),
    co_estado    = sprintf("%02d", co_municipio %/% 10000L),
    age_grp      = age_label_map[age_group]
  ) %>%
  filter(!is.na(age_grp)) %>%
  left_join(estado_name_map, by = "co_estado") %>%
  filter(!is.na(uf)) %>%
  group_by(uf, year, age_grp) %>%
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

# ── 5. Extrapolate beyond the last observed year ────────────────────
last_real_year <- max(pop_state_year_age$year)
anchor_year    <- last_real_year - 2L   # same 2-year lookback as C1's extrapolation

message("Real population data observed through ", last_real_year,
        " -- extrapolating ", last_real_year + 1L, " to ", YEAR_TARGET_MAX,
        " using the ", anchor_year, "-", last_real_year, " trend.")

pop_extrap <- pop_state_year_age %>%
  filter(year %in% c(anchor_year, last_real_year)) %>%
  group_by(uf, age_grp) %>%
  summarise(
    slope      = (pop[year == last_real_year] - pop[year == anchor_year]) /
      (last_real_year - anchor_year),
    pop_anchor = pop[year == last_real_year],
    .groups    = "drop"
  ) %>%
  crossing(year = (last_real_year + 1L):YEAR_TARGET_MAX) %>%
  mutate(pop = pmax(0L, as.integer(round(pop_anchor + slope * (year - last_real_year))))) %>%
  select(uf, year, age_grp, pop)

pop_state_full <- bind_rows(
  pop_state_year_age %>% mutate(source = "observed"),
  pop_extrap          %>% mutate(source = "extrapolated")
) %>%
  filter(year >= YEAR_MIN, year <= YEAR_TARGET_MAX) %>%
  mutate(age_grp = factor(age_grp, levels = age_labels)) %>%
  arrange(uf, year, age_grp)

# ── 6. Sanity check & save ─────────────────────────────────────────
pop_state_full %>%
  distinct(year, source) %>%
  arrange(year) %>%
  print(n = 50)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(pop_state_full, "data/processed/pop_state_age5_2000_2035.csv")
saveRDS(pop_state_full,   "data/processed/pop_state_age5_2000_2035.rds")