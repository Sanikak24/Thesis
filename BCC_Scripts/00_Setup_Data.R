# ============================================================
# Project : sigNATURE Framework
# Script  : 00_Setup_Data.R
# Purpose : Download the public inputs used for BCC atlas mapping.
#
# Outputs:
#   - data/bcc/GSE123813_bcc_scRNA_counts.txt.gz
#   - data/bcc/GSE123813_bcc_tcell_metadata.txt.gz
#   - data/bcc/41591_2019_522_MOESM2_ESM.xlsx
#   - data/reference/CD8_Obj_for_mapping.rds
#   - data/reference/CD4_Obj_for_mapping.rds
# ============================================================

options(timeout = max(3600, getOption("timeout")))

bcc_data_dir <- file.path("data", "bcc")
ref_data_dir <- file.path("data", "reference")

message("Creating data directories...")
dir.create(bcc_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_data_dir, recursive = TRUE, showWarnings = FALSE)

required_dirs <- c(bcc_data_dir, ref_data_dir)
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs) > 0) {
  stop(
    "Could not create required data directory/directories: ",
    paste(missing_dirs, collapse = ", ")
  )
}

download_input <- function(name, url, destination, expected_md5 = NULL) {
  if (file.exists(destination) && file.info(destination)$size > 0) {
    message("Skipping download; nonempty file already exists: ", destination)
  } else {
    partial_file <- paste0(destination, ".part")

    if (file.exists(partial_file)) {
      message("Partial download found; resuming from: ", partial_file)
    }

    message("Downloading ", name, " from ", url)

    curl_path <- Sys.which("curl")
    if (nzchar(curl_path)) {
      download_status <- system2(
        command = curl_path,
        args = c(
          "--location",
          "--fail",
          "--retry", "5",
          "--retry-delay", "10",
          "--retry-all-errors",
          "--continue-at", "-",
          "--output", shQuote(partial_file),
          shQuote(url)
        )
      )

      if (!identical(download_status, 0L)) {
        stop(
          "Failed to download ", name,
          "; curl returned status ", download_status,
          ". Partial download retained at: ", partial_file
        )
      }
    } else {
      message("System curl is unavailable; using download.file().")
      download_status <- tryCatch(
        download.file(
          url = url,
          destfile = partial_file,
          mode = "wb",
          method = "libcurl",
          quiet = FALSE
        ),
        error = function(error) {
          stop(
            "Failed to download ", name, ": ", conditionMessage(error),
            ". Partial download retained at: ", partial_file
          )
        }
      )

      if (!identical(download_status, 0L)) {
        stop(
          "Failed to download ", name,
          "; download.file returned status ", download_status,
          ". Partial download retained at: ", partial_file
        )
      }
    }

    if (!file.exists(partial_file) || file.info(partial_file)$size == 0) {
      stop(
        "Downloaded ", name,
        " is missing or empty. Expected partial file: ", partial_file
      )
    }

    if (!is.null(expected_md5)) {
      downloaded_md5 <- unname(tools::md5sum(partial_file))
      if (!identical(tolower(downloaded_md5), tolower(expected_md5))) {
        stop(
          "MD5 mismatch for downloaded ", name, ". Expected ",
          expected_md5, " but received ", downloaded_md5,
          ". Partial download retained at: ", partial_file
        )
      }
    }

    if (file.exists(destination) && unlink(destination) != 0) {
      stop("Could not remove empty destination file: ", destination)
    }

    if (!file.rename(partial_file, destination)) {
      stop("Validated ", name, " could not be moved to: ", destination)
    }
    message("Saved validated ", name, " to: ", destination)
  }

  if (!is.null(expected_md5)) {
    observed_md5 <- unname(tools::md5sum(destination))
    if (!identical(tolower(observed_md5), tolower(expected_md5))) {
      stop(
        "MD5 mismatch for ", destination, ". Expected ",
        expected_md5, " but received ", observed_md5, "."
      )
    }

    message("MD5 validation passed for: ", destination)
  }
}

# MD Anderson T Cell Map reference objects:
# https://singlecell.mdanderson.org/TCM/
inputs <- list(
  list(
    name = "CD8 reference",
    url = "https://singlecell.mdanderson.org/TCM/download/CD8",
    destination = file.path(ref_data_dir, "CD8_Obj_for_mapping.rds"),
    md5 = "52c6daf010f16020a8b36fb40bb95618"
  ),
  list(
    name = "CD4 reference",
    url = "https://singlecell.mdanderson.org/TCM/download/CD4",
    destination = file.path(ref_data_dir, "CD4_Obj_for_mapping.rds"),
    md5 = "69b4164d2672b39d70c61e574bcdcaa4"
  ),

  # Yost et al. BCC count matrix and T-cell metadata from GEO GSE123813:
  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813
  list(
    name = "GSE123813 BCC count matrix",
    url = paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/",
      "?acc=GSE123813&file=GSE123813_bcc_scRNA_counts.txt.gz&format=file"
    ),
    destination = file.path(
      bcc_data_dir,
      "GSE123813_bcc_scRNA_counts.txt.gz"
    )
  ),
  list(
    name = "GSE123813 BCC T-cell metadata",
    url = paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/",
      "?acc=GSE123813&file=GSE123813_bcc_tcell_metadata.txt.gz&format=file"
    ),
    destination = file.path(
      bcc_data_dir,
      "GSE123813_bcc_tcell_metadata.txt.gz"
    )
  ),

  # Yost et al. Nature Medicine supplementary tables:
  # https://doi.org/10.1038/s41591-019-0522-3
  list(
    name = "Yost clinical supplementary tables",
    url = paste0(
      "https://static-content.springer.com/esm/",
      "art%3A10.1038%2Fs41591-019-0522-3/",
      "MediaObjects/41591_2019_522_MOESM2_ESM.xlsx"
    ),
    destination = file.path(
      bcc_data_dir,
      "41591_2019_522_MOESM2_ESM.xlsx"
    )
  )
)

for (input in inputs) {
  download_input(
    name = input$name,
    url = input$url,
    destination = input$destination,
    expected_md5 = input$md5
  )
}

message("BCC input-data setup complete.")
