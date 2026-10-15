# Aggregate SINAN dengue notifications into confirmed case counts by
# year x five-year age group x state, 2000-2025.
#
# Two file layouts exist:
#   2000-2006: NU_IDADE / CON_CLASSI / CON_CRITER, SG_UF_NOT = state abbreviation
#   2007-2025: NU_IDADE_N / CLASSI_FIN / CRITERIO, SG_UF_NOT = IBGE numeric code
#
# Input : data/raw/sinan/DENGBRyy.csv  (R/01_download_sinan.R)
# Output: data/processed/dengue_cases_year_age5_state_2000_2025.csv

source("R/00_common.R")

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# -- Age parsing ----------------------------------------------------
# New layout: 4-digit code, first digit is the unit
# (4 = years, 3 = months, 2 = days, 1 = hours)
parse_age_new <- function(x) {
  x     <- suppressWarnings(as.integer(x))
  unit  <- x %/% 1000L
  value <- x %% 1000L
  case_when(
    unit == 4L              ~ as.numeric(value),
    unit %in% c(1L, 2L, 3L) ~ 0,
    TRUE                    ~ NA_real_
  )
}

# Old layout: unit letter + 3 digits (A = years, M = months, D = days, H = hours)
parse_age_old <- function(x) {
  unit  <- substr(x, 1, 1)
  value <- suppressWarnings(as.integer(substr(x, 2, 4)))
  case_when(
    unit == "A"                ~ as.numeric(value),
    unit %in% c("M", "D", "H") ~ 0,
    TRUE                       ~ NA_real_
  )
}

# -- Read one yearly file into a common layout ---------------------
read_one_csv <- function(f) {
  message("Reading ", basename(f), " ...")
  all_cols <- names(fread(f, nrows = 0))
  is_old   <- ("NU_IDADE" %in% all_cols) && !("NU_IDADE_N" %in% all_cols)

  if (is_old) {
    need <- intersect(all_cols, c("DT_NOTIFIC", "NU_IDADE", "SG_UF_NOT",
                                  "SG_UF", "ID_MN_RESI", "CON_CLASSI", "CON_CRITER"))
    fread(f, select = need, colClasses = "character") |>
      as_tibble() |>
      rename(CLASSI_FIN = CON_CLASSI, CRITERIO = CON_CRITER) |>
      mutate(NU_IDADE_N = as.character(as.integer(parse_age_old(NU_IDADE)))) |>
      select(-NU_IDADE)
  } else {
    need <- intersect(all_cols, c("DT_NOTIFIC", "NU_IDADE_N", "SG_UF_NOT",
                                  "SG_UF", "ID_MN_RESI", "CLASSI_FIN", "CRITERIO"))
    fread(f, select = need, colClasses = "character") |>
      as_tibble()
  }
}

csv_files <- list.files(DIR_RAW_SINAN, pattern = "\\.csv$",
                        full.names = TRUE, ignore.case = TRUE)
if (length(csv_files) == 0) stop("No CSV files in ", DIR_RAW_SINAN, "; run R/01_download_sinan.R first.")

df_raw <- map(csv_files, read_one_csv) |>
  bind_rows() |>
  mutate(across(everything(), as.character))
cat("Rows read:", nrow(df_raw), "\n")

# -- Keep confirmed cases (laboratory or clinical-epidemiological) --
confirmed_codes <- c("1", "2", "3", "4", "10", "11", "12")

df_filtered <- df_raw |>
  filter(CLASSI_FIN %in% confirmed_codes,
         CRITERIO   %in% c("1", "2"))
cat("Rows after confirmation filter:", nrow(df_filtered), "\n")

# -- Year, age group and state --------------------------------------
df_clean <- df_filtered |>
  mutate(
    year = as.integer(substr(DT_NOTIFIC, 1, 4)),
    age_years = case_when(
      # new layout: 4-digit code (may carry a trailing ".0")
      grepl("^[0-9]{4}(\\.0+)?$", NU_IDADE_N) ~
        parse_age_new(as.character(as.integer(as.numeric(NU_IDADE_N)))),
      # old layout: already converted to age in years
      grepl("^[0-9]{1,3}(\\.0+)?$", NU_IDADE_N) ~ as.numeric(NU_IDADE_N),
      TRUE ~ NA_real_
    ),
    age_group = cut(age_years, breaks = AGE_BREAKS, labels = AGE_LABELS,
                    right = FALSE, include.lowest = TRUE),
    uf = case_when(
      grepl("^[A-Z]{2}$", SG_UF_NOT) ~ SG_UF_NOT,
      TRUE ~ UF_CODES$uf[match(sprintf("%02d", as.integer(SG_UF_NOT)), UF_CODES$co_uf)]
    )
  ) |>
  filter(!is.na(year), year >= YEAR_MIN, year <= YEAR_MAX,
         !is.na(age_group),
         !is.na(uf))

cat("\nConfirmed cases by year:\n")
df_clean |> count(year) |> print(n = 30)

# -- Aggregate: year x age group x state (zero-filled) -------------
case_summary <- df_clean |>
  count(year, age_group, uf, name = "cases") |>
  complete(year = YEAR_MIN:YEAR_MAX, age_group, uf = UF_CODES$uf,
           fill = list(cases = 0L)) |>
  arrange(uf, year, age_group) |>
  mutate(age_group = as.character(age_group)) |>
  select(uf, year, age_group, cases)

stopifnot(n_distinct(case_summary$uf) == 27,
          n_distinct(case_summary$age_group) == length(AGE_LABELS))

dir.create(DIR_PROCESSED, showWarnings = FALSE, recursive = TRUE)
write_csv(case_summary, FILE_CASES)
cat("\nSaved:", FILE_CASES, "\n")
