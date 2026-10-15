# Shared settings for the R pipeline. Every script sources this file and is
# meant to be run from the project root (open denv_cat_sir_2025-2026.Rproj).

if (!file.exists("denv_cat_sir_2025-2026.Rproj")) {
  stop("Run the R scripts from the project root (open denv_cat_sir_2025-2026.Rproj).")
}

YEAR_MIN <- 2000
YEAR_MAX <- 2025

AGE_LABELS <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
                "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
                "65-69", "70-74", "75-79", "80+")
AGE_BREAKS <- c(seq(0, 80, by = 5), Inf)

# IBGE two-digit state code -> state abbreviation
UF_CODES <- data.frame(
  co_uf = c("11", "12", "13", "14", "15", "16", "17",
            "21", "22", "23", "24", "25", "26", "27", "28", "29",
            "31", "32", "33", "35",
            "41", "42", "43",
            "50", "51", "52", "53"),
  uf    = c("RO", "AC", "AM", "RR", "PA", "AP", "TO",
            "MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA",
            "MG", "ES", "RJ", "SP",
            "PR", "SC", "RS",
            "MS", "MT", "GO", "DF"),
  stringsAsFactors = FALSE
)

# Paths
DIR_RAW_SINAN <- file.path("data", "raw", "sinan")
DIR_RAW_POP   <- file.path("data", "raw", "population")
DIR_PROCESSED <- file.path("data", "processed")
DIR_RESULTS   <- "results"
DIR_BIG       <- "results_big"

FILE_CASES      <- file.path(DIR_PROCESSED, "dengue_cases_year_age5_state_2000_2025.csv")
FILE_POP        <- file.path(DIR_PROCESSED, "pop_state_age5_2000_2025.csv")
FILE_STAN_MODEL <- file.path("stan", "dengue_catalytic.stan")
FILE_FIT        <- file.path(DIR_BIG, "catalytic_fit.rds")
FILE_PPC_MAT    <- file.path(DIR_BIG, "catalytic_ppc_draws.mat")
FILE_PARAM_MAT  <- file.path(DIR_BIG, "catalytic_param_draws.mat")
