#!/usr/bin/env Rscript

# ==============================================================================
# NSCLC data setup
#
# Downloads the Liu et al. NSCLC T-cell data from GEO accession GSE179994.
# Reuses the CD8 and CD4 T-cell reference atlases under data/reference/.
# Downloads the supplementary tables associated with the Liu et al. study:
# https://doi.org/10.1038/s43018-021-00292-8
#
# The response-information worksheet is extracted automatically and saved as:
# data/nsclc/response_info.xlsx
#
# Run from the repository root:
#   Rscript NSCLC_Scripts/00_Setup_Data.R
# ==============================================================================

options(timeout = max(3600, getOption("timeout")))

# ------------------------------------------------------------------------------
# Repository directories
# ------------------------------------------------------------------------------

nsclc_dir <- file.path("data", "nsclc")
reference_dir <- file.path("data", "reference")
results_dir <- file.path("results", "nsclc")

dir.create(nsclc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

is_valid_file <- function(path) {
  file.exists(path) &&
    !is.na(file.info(path)$size) &&
    file.info(path)$size > 0
}

require_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      "The R package '",
      package_name,
      "' is required.\n",
      "Install it with:\n",
      "  install.packages('",
      package_name,
      "')"
    )
  }
}

download_input <- function(url, destination, retries = 3L) {
  if (is_valid_file(destination)) {
    message("Already available: ", destination)
    return(invisible(destination))
  }
  
  part_file <- paste0(destination, ".part")
  
  if (file.exists(part_file)) {
    unlink(part_file)
  }
  
  for (attempt in seq_len(retries)) {
    message(
      "Downloading ",
      basename(destination),
      " (attempt ",
      attempt,
      " of ",
      retries,
      ")"
    )
    
    result <- tryCatch(
      {
        download.file(
          url = url,
          destfile = part_file,
          mode = "wb",
          method = "libcurl",
          quiet = FALSE
        )
        
        TRUE
      },
      error = function(e) {
        warning(
          "Download attempt failed for ",
          basename(destination),
          ": ",
          conditionMessage(e)
        )
        
        FALSE
      }
    )
    
    if (isTRUE(result) && is_valid_file(part_file)) {
      if (file.exists(destination)) {
        unlink(destination)
      }
      
      renamed <- file.rename(part_file, destination)
      
      if (!renamed) {
        stop(
          "Downloaded the file but could not rename:\n  ",
          part_file,
          "\nto:\n  ",
          destination
        )
      }
      
      message("Saved: ", destination)
      return(invisible(destination))
    }
    
    if (file.exists(part_file)) {
      unlink(part_file)
    }
  }
  
  stop(
    "Failed to download ",
    basename(destination),
    " after ",
    retries,
    " attempts."
  )
}

decompress_gzip <- function(gz_file, output_file) {
  if (is_valid_file(output_file)) {
    message("Already decompressed: ", output_file)
    return(invisible(output_file))
  }
  
  if (!is_valid_file(gz_file)) {
    stop("Compressed input is missing or empty: ", gz_file)
  }
  
  message("Decompressing: ", basename(gz_file))
  
  temporary_output <- paste0(output_file, ".part")
  
  if (file.exists(temporary_output)) {
    unlink(temporary_output)
  }
  
  input_connection <- gzfile(gz_file, open = "rb")
  output_connection <- file(temporary_output, open = "wb")
  
  connection_closed <- FALSE
  
  on.exit(
    {
      if (!connection_closed) {
        try(close(input_connection), silent = TRUE)
        try(close(output_connection), silent = TRUE)
      }
    },
    add = TRUE
  )
  
  repeat {
    buffer <- readBin(
      input_connection,
      what = "raw",
      n = 1024L * 1024L
    )
    
    if (length(buffer) == 0L) {
      break
    }
    
    writeBin(buffer, output_connection)
  }
  
  close(input_connection)
  close(output_connection)
  connection_closed <- TRUE
  
  if (!is_valid_file(temporary_output)) {
    stop("Decompression failed: ", temporary_output)
  }
  
  if (file.exists(output_file)) {
    unlink(output_file)
  }
  
  renamed <- file.rename(temporary_output, output_file)
  
  if (!renamed || !is_valid_file(output_file)) {
    stop(
      "Could not finalize decompressed file:\n  ",
      output_file
    )
  }
  
  message("Created: ", output_file)
  
  invisible(output_file)
}

require_reference <- function(path, description) {
  if (!is_valid_file(path)) {
    stop(
      description,
      " was not found at:\n  ",
      path,
      "\n\nRun the existing reference-atlas setup first."
    )
  }
  
  message("Reference available: ", path)
}

normalize_header <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("\\s+", " ", x)
  tolower(x)
}

find_response_table <- function(workbook_path) {
  require_package("readxl")
  
  sample_headers <- c("sample name", "samples")
  sheet_names <- readxl::excel_sheets(workbook_path)
  sheet_names <- c(
    intersect("Supplementary Table 2", sheet_names),
    setdiff(sheet_names, "Supplementary Table 2")
  )
  
  message(
    "Searching supplementary workbook sheets: ",
    paste(sheet_names, collapse = ", ")
  )
  
  for (sheet_name in sheet_names) {
    preview <- tryCatch(
      readxl::read_excel(
        workbook_path,
        sheet = sheet_name,
        col_names = FALSE,
        n_max = 50,
        .name_repair = "minimal"
      ),
      error = function(e) {
        warning(
          "Could not inspect worksheet '",
          sheet_name,
          "': ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (is.null(preview) || nrow(preview) == 0L) {
      next
    }
    
    preview_matrix <- as.matrix(preview)
    
    for (row_index in seq_len(nrow(preview_matrix))) {
      row_values <- normalize_header(preview_matrix[row_index, ])
      
      has_required_headers <-
        any(sample_headers %in% row_values) &&
        "response" %in% row_values
      
      if (!has_required_headers) {
        next
      }
      
      message(
        "Found response table in worksheet '",
        sheet_name,
        "' at Excel row ",
        row_index,
        "."
      )
      
      response_table <- readxl::read_excel(
        workbook_path,
        sheet = sheet_name,
        skip = row_index - 1L,
        .name_repair = "unique"
      )
      
      original_names <- names(response_table)
      normalized_names <- normalize_header(original_names)
      
      sample_column_index <- match(
        TRUE,
        normalized_names %in% sample_headers
      )
      
      response_column_index <- match(
        "response",
        normalized_names
      )
      
      if (
        is.na(sample_column_index) ||
        is.na(response_column_index)
      ) {
        next
      }
      
      names(response_table)[sample_column_index] <- "Sample Name"
      names(response_table)[response_column_index] <- "Response"
      
      response_table <- response_table[
        !is.na(response_table[["Sample Name"]]) &
          trimws(as.character(
            response_table[["Sample Name"]]
          )) != "",
        ,
        drop = FALSE
      ]
      
      if (nrow(response_table) == 0L) {
        stop(
          "The response table was found, but it contains no sample rows."
        )
      }
      
      return(
        list(
          data = response_table,
          sheet = sheet_name,
          header_row = row_index
        )
      )
    }
  }
  
  stop(
    "Could not find a worksheet row containing both required columns:\n",
    "  - Sample Name or Samples\n",
    "  - Response\n\n",
    "Worksheets inspected:\n  ",
    paste(sheet_names, collapse = "\n  ")
  )
}

create_response_file <- function(source_file, destination_file) {
  require_package("readxl")
  require_package("writexl")
  
  response_result <- find_response_table(source_file)
  response_table <- response_result$data

  response_values <- trimws(
    as.character(response_table[["Response"]])
  )
  response_values[response_values == ""] <- NA_character_

  if (length(response_values) > 1L) {
    for (row_index in 2:length(response_values)) {
      if (is.na(response_values[row_index])) {
        response_values[row_index] <-
          response_values[row_index - 1L]
      }
    }
  }

  response_table[["Response"]] <- response_values
  
  required_columns <- c("Sample Name", "Response")
  
  missing_columns <- setdiff(
    required_columns,
    names(response_table)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "The extracted response table is missing required columns:\n  ",
      paste(missing_columns, collapse = "\n  ")
    )
  }
  
  temporary_file <- paste0(destination_file, ".part.xlsx")
  
  if (file.exists(temporary_file)) {
    unlink(temporary_file)
  }
  
  writexl::write_xlsx(
    list(response_info = response_table),
    path = temporary_file
  )
  
  if (!is_valid_file(temporary_file)) {
    stop(
      "Failed to create temporary response file:\n  ",
      temporary_file
    )
  }
  
  if (file.exists(destination_file)) {
    unlink(destination_file)
  }
  
  renamed <- file.rename(
    temporary_file,
    destination_file
  )
  
  if (!renamed || !is_valid_file(destination_file)) {
    stop(
      "Failed to create response metadata file:\n  ",
      destination_file
    )
  }
  
  message(
    "Created response metadata from worksheet '",
    response_result$sheet,
    "', header row ",
    response_result$header_row,
    ": ",
    destination_file
  )
  
  message(
    "Response metadata contains ",
    nrow(response_table),
    " rows and ",
    ncol(response_table),
    " columns."
  )
  
  invisible(destination_file)
}

# ------------------------------------------------------------------------------
# Shared reference atlases
# ------------------------------------------------------------------------------

cd8_reference <- file.path(
  reference_dir,
  "CD8_Obj_for_mapping.rds"
)

cd4_reference <- file.path(
  reference_dir,
  "CD4_Obj_for_mapping.rds"
)

require_reference(
  cd8_reference,
  "The shared CD8 reference atlas"
)

require_reference(
  cd4_reference,
  "The shared CD4 reference atlas"
)

# ------------------------------------------------------------------------------
# Liu et al. NSCLC GEO supplementary files
# ------------------------------------------------------------------------------

geo_base_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE179nnn/GSE179994/suppl/"
)

geo_files <- list(
  list(
    compressed_name =
      "GSE179994_all.Tcell.rawCounts.rds.gz",
    output_name =
      "GSE179994_all.Tcell.rawCounts.rds"
  ),
  list(
    compressed_name =
      "GSE179994_Tcell.metadata.tsv.gz",
    output_name =
      "GSE179994_Tcell.metadata.tsv"
  )
)

for (input_file in geo_files) {
  compressed_path <- file.path(
    nsclc_dir,
    input_file$compressed_name
  )
  
  output_path <- file.path(
    nsclc_dir,
    input_file$output_name
  )
  
  download_input(
    url = paste0(
      geo_base_url,
      input_file$compressed_name
    ),
    destination = compressed_path
  )
  
  decompress_gzip(
    gz_file = compressed_path,
    output_file = output_path
  )
}

# ------------------------------------------------------------------------------
# Clinical-response supplementary workbook
# ------------------------------------------------------------------------------

response_source_name <- "43018_2021_292_MOESM3_ESM.xlsx"

response_source_url <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/",
  "art%3A10.1038%2Fs43018-021-00292-8/",
  "MediaObjects/",
  response_source_name
)

response_source_file <- file.path(
  nsclc_dir,
  response_source_name
)

response_file <- file.path(
  nsclc_dir,
  "response_info.xlsx"
)

if (!is_valid_file(response_source_file)) {
  download_input(
    url = response_source_url,
    destination = response_source_file
  )
} else {
  message(
    "Supplementary workbook available: ",
    response_source_file
  )
}

# Always regenerate response_info.xlsx from the correct worksheet so that
# it contains a clean response table rather than a renamed full workbook.
create_response_file(
  source_file = response_source_file,
  destination_file = response_file
)

# ------------------------------------------------------------------------------
# Final validation
# ------------------------------------------------------------------------------

required_files <- c(
  cd8_reference,
  cd4_reference,
  file.path(
    nsclc_dir,
    "GSE179994_all.Tcell.rawCounts.rds"
  ),
  file.path(
    nsclc_dir,
    "GSE179994_Tcell.metadata.tsv"
  ),
  response_source_file,
  response_file
)

missing_files <- required_files[
  !vapply(
    required_files,
    is_valid_file,
    logical(1)
  )
]

if (length(missing_files) > 0L) {
  stop(
    "Setup did not produce all required NSCLC files:\n",
    paste0(
      "  - ",
      missing_files,
      collapse = "\n"
    )
  )
}

message("")
message("NSCLC data setup completed successfully.")
message("NSCLC inputs: ", nsclc_dir)
message("Shared references: ", reference_dir)
message("Results directory: ", results_dir)
message("")
message("Required files:")
message(
  paste0(
    "  - ",
    required_files,
    collapse = "\n"
  )
)
