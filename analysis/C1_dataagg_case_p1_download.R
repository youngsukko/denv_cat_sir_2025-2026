# Download SINAN dengue per-year CSVs from OpenDataSUS (S3, not the
# geo-restricted DataSUS FTP) into the folder the main script reads from.

out_dir <- "data/case/dengue_csv_raw"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (yr in 2000:2025) {
  yy  <- sprintf("%02d", yr %% 100)
  url <- sprintf("https://s3.sa-east-1.amazonaws.com/ckan.saude.gov.br/SINAN/Dengue/csv/DENGBR%s.csv.zip", yy)
  zip_path <- file.path(out_dir, sprintf("DENGBR%s.csv.zip", yy))
  
  message("Downloading ", yr, " ...")
  ok <- tryCatch(
    download.file(url, zip_path, mode = "wb", quiet = TRUE) == 0,
    error = function(e) {
      message("  -> not available (", conditionMessage(e), ")")
      FALSE
    }
  )
  
  if (ok) unzip(zip_path, exdir = out_dir)
}