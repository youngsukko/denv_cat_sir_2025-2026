# ══════════════════════════════════════════════════════════════════
# SINAN Dengue — Year x 5-year Age Group x State case counts
# 2000-2025, confirmed cases only
# 구 양식(2000-2006): NU_IDADE / CON_CLASSI / CON_CRITER
#                     SG_UF_NOT = UF 약자 (AC, AL, ...)
# 신 양식(2007-2025): NU_IDADE_N / CLASSI_FIN / CRITERIO
#                     SG_UF_NOT = 숫자 코드 (11, 12, ...)
# ══════════════════════════════════════════════════════════════════

library(data.table)
library(tidyverse)
library(lubridate)

YEAR_MIN     <- 2000
YEAR_MAX     <- 2025
raw_data_dir <- "data/case/dengue_csv_raw"

age_breaks <- c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,Inf)
age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                "35-39","40-44","45-49","50-54","55-59","60-64",
                "65-69","70-74","75-79","80+")

estado_name_map <- tibble(
  co_estado = c("11","12","13","14","15","16","17",
                "21","22","23","24","25","26","27","28","29",
                "31","32","33","35",
                "41","42","43",
                "50","51","52","53"),
  uf = c("RO","AC","AM","RR","PA","AP","TO",
         "MA","PI","CE","RN","PB","PE","AL","SE","BA",
         "MG","ES","RJ","SP",
         "PR","SC","RS",
         "MS","MT","GO","DF")
)

# ── 나이 파싱 함수 ──────────────────────────────────────────────
# 신 양식: 4자리 숫자 (4xxx = 년, 3xxx = 월, 2xxx = 일, 1xxx = 시)
parse_age_new <- function(x) {
  x     <- suppressWarnings(as.integer(x))
  unit  <- x %/% 1000L
  value <- x %% 1000L
  case_when(
    unit == 4L             ~ as.numeric(value),
    unit %in% c(1L,2L,3L) ~ 0,
    TRUE                   ~ NA_real_
  )
}

# 구 양식: A/M/D/H + 3자리 (A=년, M=월, D=일, H=시)
parse_age_old <- function(x) {
  unit  <- substr(x, 1, 1)
  value <- suppressWarnings(as.integer(substr(x, 2, 4)))
  case_when(
    unit == "A"               ~ as.numeric(value),
    unit %in% c("M","D","H") ~ 0,
    TRUE                      ~ NA_real_
  )
}

# ── 파일별 로드 함수 ────────────────────────────────────────────
read_one_csv <- function(f) {
  message("Reading ", basename(f), " ...")
  all_cols <- names(fread(f, nrows = 0))
  is_old   <- ("NU_IDADE" %in% all_cols) & !("NU_IDADE_N" %in% all_cols)
  
  if (is_old) {
    # 구 양식 (2000-2006)
    need <- intersect(all_cols,
                      c("DT_NOTIFIC","NU_IDADE","SG_UF_NOT",
                        "SG_UF","ID_MN_RESI","CON_CLASSI","CON_CRITER"))
    fread(f, select = need, colClasses = "character") |>
      as_tibble() |>
      rename(CLASSI_FIN = CON_CLASSI,
             CRITERIO   = CON_CRITER) |>
      mutate(NU_IDADE_N = as.character(
        as.integer(parse_age_old(NU_IDADE)))) |>
      select(-NU_IDADE)
  } else {
    # 신 양식 (2007-2025)
    slim <- c("DT_NOTIFIC","NU_IDADE_N","SG_UF_NOT",
              "SG_UF","ID_MN_RESI","CLASSI_FIN","CRITERIO")
    need <- intersect(all_cols, slim)
    fread(f, select = need, colClasses = "character") |>
      as_tibble()
  }
}

# ── 전체 파일 로드 ──────────────────────────────────────────────
csv_files <- list.files(raw_data_dir, pattern = "\\.csv$",
                        full.names = TRUE, ignore.case = TRUE)
if (length(csv_files) == 0)
  stop("CSV 파일 없음: ", raw_data_dir)

cat("파일 수:", length(csv_files), "\n")
df_raw <- map(csv_files, read_one_csv) |>
  bind_rows() |>
  mutate(across(everything(), as.character))
cat("총 행수 (필터 전):", nrow(df_raw), "\n")

# ── Confirmed 필터 ───────────────────────────────────────────────
confirmed_codes <- c("1","2","3","4","10","11","12")

df_filtered <- df_raw |>
  filter(CLASSI_FIN %in% confirmed_codes,
         CRITERIO   %in% c("1","2"))

cat("총 행수 (confirmed 필터 후):", nrow(df_filtered), "\n")

# ── 연도 / 나이 / 주 파싱 ────────────────────────────────────────
df_clean <- df_filtered |>
  mutate(
    year = as.integer(substr(DT_NOTIFIC, 1, 4)),
    
    age_years = case_when(
      # 신 양식: 4자리 숫자 (소수점 포함 가능)
      grepl("^[0-9]{4}(\\.0+)?$", NU_IDADE_N) ~
        parse_age_new(as.character(
          as.integer(as.numeric(NU_IDADE_N)))),
      # 구 양식: 이미 정수 나이(연)로 변환됨
      grepl("^[0-9]{1,3}(\\.0+)?$", NU_IDADE_N) ~
        as.numeric(NU_IDADE_N),
      TRUE ~ NA_real_
    ),
    
    age_group = cut(age_years,
                    breaks         = age_breaks,
                    labels         = age_labels,
                    right          = FALSE,
                    include.lowest = TRUE),
    
    # 구 양식: SG_UF_NOT = UF 약자 (AC, AL, ...)
    # 신 양식: SG_UF_NOT = 숫자 코드 (11, 12, ...)
    uf = case_when(
      grepl("^[A-Z]{2}$", SG_UF_NOT) ~ SG_UF_NOT,
      TRUE ~ estado_name_map$uf[
        match(sprintf("%02d", as.integer(SG_UF_NOT)),
              estado_name_map$co_estado)]
    )
  ) |>
  filter(
    !is.na(year),      year >= YEAR_MIN, year <= YEAR_MAX,
    !is.na(age_group),
    !is.na(uf)
  )

# ── 연도별 케이스 수 확인 ─────────────────────────────────────────
cat("\n연도별 케이스 수:\n")
df_clean |>
  count(year) |>
  print(n = 30)

# ── 집계: year x age_group x uf ──────────────────────────────────
case_summary <- df_clean |>
  count(year, age_group, uf, name = "cases") |>
  complete(
    year      = YEAR_MIN:YEAR_MAX,
    age_group = age_labels,
    uf        = estado_name_map$uf,
    fill      = list(cases = 0L)
  ) |>
  arrange(uf, year, age_group)

# ── 최종 확인 ────────────────────────────────────────────────────
cat("\n최종 집계 확인:\n")
case_summary |>
  group_by(year) |>
  summarise(total = sum(cases), .groups = "drop") |>
  print(n = 30)

cat("\n주 수:", n_distinct(case_summary$uf), "\n")
cat("연령군 수:", n_distinct(case_summary$age_group), "\n")

# ── 저장 ─────────────────────────────────────────────────────────
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(case_summary, "data/processed/dengue_year_age5_state_2000_2025.csv")
saveRDS(case_summary,   "data/processed/dengue_year_age5_state_2000_2025.rds")
cat("\nSaved: data/processed/dengue_year_age5_state_2000_2025.rds\n")