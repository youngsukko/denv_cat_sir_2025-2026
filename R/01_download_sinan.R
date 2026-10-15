# Download the SINAN dengue line lists (one zipped CSV per year) from the
# OpenDataSUS S3 bucket. The DataSUS FTP server is geo-restricted, the S3
# mirror is not. The files are several GB in total and are not tracked by git.
#
# Output: data/raw/sinan/DENGBRyy.csv

source("R/00_common.R")

dir.create(DIR_RAW_SINAN, showWarnings = FALSE, recursive = TRUE)
options(timeout = max(3600, getOption("timeout")))

base_url <- "https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SINAN/Dengue/csv"

for (yr in YEAR_MIN:YEAR_MAX) {
  yy       <- sprintf("%02d", yr %% 100)
  zip_name <- sprintf("DENGBR%s.csv.zip", yy)
  zip_path <- file.path(DIR_RAW_SINAN, zip_name)

  if (file.exists(file.path(DIR_RAW_SINAN, sprintf("DENGBR%s.csv", yy)))) {
    message(yr, ": already downloaded")
    next
  }

  message("Downloading ", yr, " ...")
  ok <- tryCatch(
    download.file(file.path(base_url, zip_name), zip_path, mode = "wb", quiet = TRUE) == 0,
    error = function(e) {
      message("  -> not available (", conditionMessage(e), ")")
      FALSE
    }
  )
  if (ok) unzip(zip_path, exdir = DIR_RAW_SINAN)
}
