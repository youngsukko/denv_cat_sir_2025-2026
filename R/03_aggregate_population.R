# State x five-year age group population, 2000-2025.
#
# Municipality-level population by age (DataSUS estimates, via the brpop
# package) is aggregated to states. Years after the last observed year are
# filled by linear extrapolation using the slope over the last two years.
#
# Input : brpop::mun_pop_age(), cached in data/raw/population/
# Output: data/processed/pop_state_age5_2000_2025.csv

source("R/00_common.R")

library(dplyr)
library(tidyr)
library(readr)

# ── 1. Municipality-level population (downloaded once, then cached) ──
cache_file <- file.path(DIR_RAW_POP, "pop_mun_age_datasus.rds")
if (!file.exists(cache_file)) {
  dir.create(DIR_RAW_POP, showWarnings = FALSE, recursive = TRUE)
  pop_all_raw <- brpop::mun_pop_age(source = "datasus", sex = "all")
  saveRDS(pop_all_raw, cache_file)
}
pop_all_raw <- readRDS(cache_file)

# ── 2. brpop age labels -> five-year bands used for the case data ──
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
stopifnot(identical(unname(age_label_map), AGE_LABELS))

# ── 3. Aggregate municipalities -> state x year x age group ────────
pop_state <- pop_all_raw |>
  filter(age_group != "Total") |>
  mutate(
    co_uf         = sprintf("%02d", as.integer(code_muni) %/% 10000L),
    age_group_new = unname(age_label_map[as.character(age_group)])
  ) |>
  filter(!is.na(age_group_new)) |>
  inner_join(UF_CODES, by = "co_uf") |>
  group_by(uf, year, age_group = age_group_new) |>
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

# ── 4. Linear extrapolation beyond the last observed year ─────────
last_real_year <- max(pop_state$year)
if (last_real_year < YEAR_MAX) {
  anchor_year <- last_real_year - 2L
  message("Observed population through ", last_real_year, "; extrapolating ",
          last_real_year + 1L, "-", YEAR_MAX, " from the ", anchor_year, "-",
          last_real_year, " trend.")

  pop_extrap <- pop_state |>
    filter(year %in% c(anchor_year, last_real_year)) |>
    group_by(uf, age_group) |>
    summarise(
      slope      = (pop[year == last_real_year] - pop[year == anchor_year]) /
                   (last_real_year - anchor_year),
      pop_anchor = pop[year == last_real_year],
      .groups    = "drop"
    ) |>
    crossing(year = (last_real_year + 1L):YEAR_MAX) |>
    mutate(pop = pmax(0L, as.integer(round(pop_anchor + slope * (year - last_real_year))))) |>
    select(uf, year, age_group, pop)

  pop_state <- bind_rows(pop_state, pop_extrap)
}

pop_state <- pop_state |>
  filter(year >= YEAR_MIN, year <= YEAR_MAX) |>
  mutate(age_group = factor(age_group, levels = AGE_LABELS)) |>
  arrange(uf, year, age_group) |>
  mutate(age_group = as.character(age_group)) |>
  select(uf, year, age_group, pop)

stopifnot(nrow(pop_state) == 27 * (YEAR_MAX - YEAR_MIN + 1) * length(AGE_LABELS))

dir.create(DIR_PROCESSED, showWarnings = FALSE, recursive = TRUE)
write_csv(pop_state, FILE_POP)
cat("Saved:", FILE_POP, "\n")
