################################################################################
### Script 00: Data Preprocessing - Read Counts and TPMs
################################################################################

# Clear workspace
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(janitor)
})

################################################################################
### SECTION 1: PROCESS RAW READ COUNTS
################################################################################

message("\n========================================")
message("Processing Raw Read Count Data")
message("========================================\n")

#' Process mouse read count data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_mouse_reads <- function(input_path = "input", 
                                output_path = "processed_reads") {
  
  message("Processing mouse cell lines...")
  
  # Create output directory
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Load main mouse data
  mouse_file <- file.path(input_path, "mouse_STAR_featureCounts_raw-counts.txt")
  if (!file.exists(mouse_file)) {
    stop(paste("File not found:", mouse_file))
  }
  
  # Read with check.names=FALSE to preserve original column names
  mouse_reads <- read.table(mouse_file, header = TRUE, check.names = FALSE)
  
  # R converts hyphens to dots in column names by default
  # So KRS7_14_aTC1-6_2_Tg6h becomes KRS7_14_aTC1.6_2_Tg6h
  # We need to account for both possibilities
  
  # Handle resequenced sample
  reseq_file <- file.path(input_path, 
                          "KRS7_14_resequenced_mouse_ILMN_1541_Kalwat_mRNAseq54_Nov2022_STAR_featureCounts_raw-counts.txt")
  
  if (file.exists(reseq_file)) {
    message("  Incorporating resequenced KRS7_14 sample...")
    KRS7_14_reads <- read.table(reseq_file, header = TRUE, check.names = FALSE)
    
    # Look for both possible column names (with - or .)
    old_col_hyphen <- "KRS7_14_aTC1-6_2_Tg6h"
    old_col_dot <- "KRS7_14_aTC1.6_2_Tg6h"
    
    # Remove old KRS7_14 column (check both possible names)
    if (old_col_hyphen %in% colnames(mouse_reads)) {
      message("  Removing old KRS7_14 column (with hyphen)...")
      mouse_reads <- mouse_reads[, !colnames(mouse_reads) == old_col_hyphen]
    } else if (old_col_dot %in% colnames(mouse_reads)) {
      message("  Removing old KRS7_14 column (with dot)...")
      mouse_reads <- mouse_reads[, !colnames(mouse_reads) == old_col_dot]
    } else {
      message("  Warning: Old KRS7_14 column not found, checking all KRS7_14 patterns...")
      # Remove any KRS7_14 aTC1 column
      krs7_14_cols <- grep("KRS7_14.*aTC1", colnames(mouse_reads))
      if (length(krs7_14_cols) > 0) {
        message(paste("  Found and removing:", colnames(mouse_reads)[krs7_14_cols]))
        mouse_reads <- mouse_reads[, -krs7_14_cols]
      }
    }
    
    # Add new KRS7_14 data
    KRS7_14_simple <- KRS7_14_reads[, c("Geneid", "KRS7_14")]
    mouse_reads <- full_join(mouse_reads, KRS7_14_simple, by = "Geneid")
    
    # Rename to match expected pattern (using dot for consistency)
    names(mouse_reads)[names(mouse_reads) == "KRS7_14"] <- "KRS7_14_aTC1.6_2_Tg6h"
  }
  
  # Process each mouse cell line
  # Updated reordering indices - these now exclude the Geneid column
  # Old indices assumed column 1 was Geneid, but we're now working with data only
  mouse_lines <- list(
    aTC1= list(
      pattern = "aTC1",
      order = NULL  # determine dynamically
    ),
    GLUTag = list(
      pattern = "GLUTag",
      order = NULL  # determine dynamically
    ),
    MGN3 = list(
      pattern = "MGN3",
      order = NULL  # determine dynamically
    )
  )
  
  for (cell_line in names(mouse_lines)) {
    info <- mouse_lines[[cell_line]]
    
    # Find columns for this cell line (accounting for both - and . in names)
    cols <- grep(info$pattern, colnames(mouse_reads))
    
    if (length(cols) > 0) {
      # Extract data with Geneid column
      cell_data <- mouse_reads[, c(1, cols)]
      
      # Identify columns by treatment and time
      data_cols <- colnames(cell_data)[-1]  # Exclude Geneid
      
      # Identify DMSO, Tg6h, and Tg24h columns
      dmso_idx <- grep("DMSO", data_cols)
      tg6h_idx <- grep("Tg6h", data_cols)
      tg24h_idx <- grep("Tg24h", data_cols)
      
      message(paste("  ", cell_line, ": Found", 
                   length(dmso_idx), "DMSO,",
                   length(tg6h_idx), "Tg6h,",
                   length(tg24h_idx), "Tg24h columns"))
      
      # Reorder: Geneid, DMSO (1-3), Tg6h (1-3), Tg24h (1-3)
      if (length(dmso_idx) == 3 && length(tg6h_idx) == 3 && length(tg24h_idx) == 3) {
        # Sort each group by replicate number
        dmso_cols <- data_cols[dmso_idx]
        tg6h_cols <- data_cols[tg6h_idx]
        tg24h_cols <- data_cols[tg24h_idx]
        
        # Extract replicate numbers and sort
        get_replicate_num <- function(col_names) {
          nums <- as.numeric(gsub(".*_(\\d)_.*", "\\1", col_names))
          col_names[order(nums)]
        }
        
        dmso_cols <- get_replicate_num(dmso_cols)
        tg6h_cols <- get_replicate_num(tg6h_cols)
        tg24h_cols <- get_replicate_num(tg24h_cols)
        
        # Reorder the dataframe
        cell_data <- cell_data[, c("Geneid", dmso_cols, tg6h_cols, tg24h_cols)]
        
        message(paste("  ", cell_line, ": Columns reordered successfully"))
      } else {
        warning(paste("  ", cell_line, ": Unexpected number of columns, keeping original order"))
      }
      
      # Save to file
      output_file <- file.path(output_path, paste0(cell_line, "_genes.txt"))
      write.table(cell_data, file = output_file, sep = '\t', 
                 row.names = FALSE, quote = FALSE)
      message(paste("  ✓", cell_line, "processed -", ncol(cell_data)-1, "samples"))
    }
  }
  
  return(invisible(TRUE))
}

#' Process MIN6 read count data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_MIN6_reads <- function(input_path = "input", 
                              output_path = "processed_reads") {
  
  message("Processing MIN6 cell line...")
  
  min6_file <- file.path(input_path, "Kalwat_TotalRNAseq_processeddata.xlsx")
  
  if (!file.exists(min6_file)) {
    warning(paste("MIN6 file not found:", min6_file))
    return(invisible(FALSE))
  }
  
  # Read MIN6 data
  MIN6_reads <- read_excel(min6_file) %>%
    janitor::clean_names() %>%
    dplyr::rename(Geneid = gene_symbol) %>%
    dplyr::select(contains(c("Geneid", "dmso_24h", "thapsigargin_6h", "thapsigargin_24h")))
  
  # Ensure proper column ordering: Geneid, DMSO_24h (1-3), Tg_6h (1-3), Tg_24h (1-3)
  # The Excel file should already be in the correct order, but let's verify
  col_order <- c("Geneid",
                 "dmso_24h_n_1", "dmso_24h_n_2", "dmso_24h_n_3",
                 "thapsigargin_6h_n_1", "thapsigargin_6h_n_2", "thapsigargin_6h_n_3",
                 "thapsigargin_24h_n_1", "thapsigargin_24h_n_2", "thapsigargin_24h_n_3")
  
  if (all(col_order %in% colnames(MIN6_reads))) {
    MIN6_reads <- MIN6_reads[, col_order]
  }
  
  # Save to file
  output_file <- file.path(output_path, "MIN6_genes.txt")
  write.table(MIN6_reads, file = output_file, sep = '\t', 
             row.names = FALSE, quote = FALSE)
  
  message(paste("  ✓ MIN6 processed -", ncol(MIN6_reads)-1, "samples"))
  
  return(invisible(TRUE))
}

#' Process rat (PCCL3) read count data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_rat_reads <- function(input_path = "input", 
                              output_path = "processed_reads") {
  
  message("Processing PCCL3 (rat) cell line...")
  
  rat_file <- file.path(input_path, "rat_STAR_featureCounts_raw-counts.txt")
  
  if (!file.exists(rat_file)) {
    warning(paste("Rat file not found:", rat_file))
    return(invisible(FALSE))
  }
  
  # Read rat data
  rat_reads <- read.table(rat_file, header = TRUE, check.names = FALSE)
  
  # Extract PCCL3 columns
  PCCL3_cols <- grep("PCCL3", colnames(rat_reads))
  
  if (length(PCCL3_cols) > 0) {
    PCCL3_reads <- rat_reads[, c(1, PCCL3_cols)]
    
    # Identify and reorder columns
    data_cols <- colnames(PCCL3_reads)[-1]
    
    dmso_idx <- grep("DMSO", data_cols)
    tg6h_idx <- grep("Tg6h", data_cols)
    tg24h_idx <- grep("Tg24h", data_cols)
    
    if (length(dmso_idx) == 3 && length(tg6h_idx) == 3 && length(tg24h_idx) == 3) {
      # Sort each group by replicate number
      get_replicate_num <- function(col_names) {
        nums <- as.numeric(gsub(".*_(\\d)_.*", "\\1", col_names))
        col_names[order(nums)]
      }
      
      dmso_cols <- get_replicate_num(data_cols[dmso_idx])
      tg6h_cols <- get_replicate_num(data_cols[tg6h_idx])
      tg24h_cols <- get_replicate_num(data_cols[tg24h_idx])
      
      PCCL3_reads <- PCCL3_reads[, c("Geneid", dmso_cols, tg6h_cols, tg24h_cols)]
    }
    
    # Save to file
    output_file <- file.path(output_path, "PCCL3_genes.txt")
    write.table(PCCL3_reads, file = output_file, sep = '\t', 
               row.names = FALSE, quote = FALSE)
    
    message(paste("  ✓ PCCL3 processed -", ncol(PCCL3_reads)-1, "samples"))
  }
  
  return(invisible(TRUE))
}

#' Process human (QGP1) read count data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_human_reads <- function(input_path = "input", 
                               output_path = "processed_reads") {
  
  message("Processing QGP1 (human) cell line...")
  
  human_file <- file.path(input_path, "human_STAR_featureCounts_raw-counts.txt")
  
  if (!file.exists(human_file)) {
    warning(paste("Human file not found:", human_file))
    return(invisible(FALSE))
  }
  
  # Read human data
  human_reads <- read.table(human_file, header = TRUE, check.names = FALSE)
  
  # QGP1 columns should be at specific positions
  # Identify by sample names
  qgp1_cols <- grep("QGP1", colnames(human_reads))
  
  if (length(qgp1_cols) > 0) {
    QGP1_reads <- human_reads[, c(1, qgp1_cols)]
    
    # Identify and reorder columns
    data_cols <- colnames(QGP1_reads)[-1]
    
    dmso_idx <- grep("DMSO", data_cols)
    tg6h_idx <- grep("Tg6h", data_cols)
    tg24h_idx <- grep("Tg24h", data_cols)
    
    if (length(dmso_idx) == 3 && length(tg6h_idx) == 3 && length(tg24h_idx) == 3) {
      # Sort each group by replicate number
      get_replicate_num <- function(col_names) {
        nums <- as.numeric(gsub(".*_(\\d)_.*", "\\1", col_names))
        col_names[order(nums)]
      }
      
      dmso_cols <- get_replicate_num(data_cols[dmso_idx])
      tg6h_cols <- get_replicate_num(data_cols[tg6h_idx])
      tg24h_cols <- get_replicate_num(data_cols[tg24h_idx])
      
      QGP1_reads <- QGP1_reads[, c("Geneid", dmso_cols, tg6h_cols, tg24h_cols)]
    }
  } else {
    # Fallback to position-based selection if pattern matching fails
    # Based on: columns 7, 10, 13 (DMSO), 8, 11, 14 (Tg6h), 9, 12, 15 (Tg24h)
    if (ncol(human_reads) >= 15) {
      QGP1_reads <- human_reads[, c(1, 7, 10, 13, 8, 11, 14, 9, 12, 15)]
    }
  }
  
  # Save to file
  output_file <- file.path(output_path, "QGP1_genes.txt")
  write.table(QGP1_reads, file = output_file, sep = '\t', 
             row.names = FALSE, quote = FALSE)
  
  message(paste("  ✓ QGP1 processed -", ncol(QGP1_reads)-1, "samples"))
  
  return(invisible(TRUE))
}

################################################################################
### SECTION 2: PROCESS TPM DATA
################################################################################

message("\n========================================")
message("Processing TPM Data")
message("========================================\n")

#' Process MIN6 TPM data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_MIN6_TPMs <- function(input_path = "input/TPMs", 
                             output_path = "processed_TPMs") {
  
  message("Processing MIN6 TPMs...")
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  min6_file <- file.path(input_path, "ILMN_989_Kalwat_IBRI_TotalRNAseq39_Mar2021_fixcolnames.xlsx")
  
  if (!file.exists(min6_file)) {
    # Try alternative path
    min6_file <- file.path("input", "TPMs", "ILMN_989_Kalwat_IBRI_TotalRNAseq39_Mar2021_fixcolnames.xlsx")
  }
  
  if (!file.exists(min6_file)) {
    warning(paste("MIN6 TPM file not found"))
    return(invisible(FALSE))
  }
  
  # Read MIN6 TPMs
  MIN6_TPMs <- read_excel(min6_file, sheet = "Counts") %>%
    janitor::clean_names() %>%
    dplyr::select(contains(c("gene_symbol", "dmso_24h", "thapsigargin_6h", "thapsigargin_24h")))
  
  # Save to file
  output_file <- file.path(output_path, "MIN6_TPMs.txt")
  write.table(MIN6_TPMs, file = output_file, sep = '\t', 
             row.names = FALSE, quote = FALSE)
  
  message("  ✓ MIN6 TPMs processed")
  
  return(invisible(TRUE))
}

#' Process mouse TPM data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_mouse_TPMs <- function(input_path = "input/TPMs", 
                              output_path = "processed_TPMs") {
  
  message("Processing mouse cell line TPMs...")
  
  mouse_file <- file.path(input_path, "mouse_ILMN_1541_Kalwat_mRNAseq54_Nov2022_TPM.txt")
  
  if (!file.exists(mouse_file)) {
    warning(paste("Mouse TPM file not found:", mouse_file))
    return(invisible(FALSE))
  }
  
  # Read mouse TPMs with check.names=FALSE to preserve column names
  mouse_TPMs <- read.table(mouse_file, header = TRUE, check.names = FALSE)
  
  # Handle resequenced sample if exists
  reseq_file <- file.path(input_path, 
                          "KRS7_14_resequenced_mouse_ILMN_1541_Kalwat_mRNAseq54_Nov2022_TPM.txt")
  
  if (file.exists(reseq_file)) {
    KRS7_14_TPMs <- read.table(reseq_file, header = TRUE, check.names = FALSE)
    
    # Remove old KRS7_14 column (accounting for both - and . in names)
    old_cols <- grep("KRS7_14.*aTC1", colnames(mouse_TPMs))
    if (length(old_cols) > 0) {
      mouse_TPMs <- mouse_TPMs[, -old_cols]
    }
    
    mouse_TPMs <- full_join(mouse_TPMs, KRS7_14_TPMs, by = "Geneid")
    names(mouse_TPMs)[names(mouse_TPMs) == "KRS7_14"] <- "KRS7_14_aTC1.6_2_Tg6h"
  }
  
  # Process each mouse cell line
  for (cell_line in c("MGN3", "aTC1", "GLUTag")) {
    cols <- grep(cell_line, colnames(mouse_TPMs))
    
    if (length(cols) > 0) {
      cell_TPMs <- mouse_TPMs[, c(1, cols)]
      
      # Identify and reorder columns
      data_cols <- colnames(cell_TPMs)[-1]
      
      dmso_idx <- grep("DMSO", data_cols)
      tg6h_idx <- grep("Tg6h", data_cols)
      tg24h_idx <- grep("Tg24h", data_cols)
      
      if (length(dmso_idx) == 3 && length(tg6h_idx) == 3 && length(tg24h_idx) == 3) {
        # Sort by replicate number
        get_replicate_num <- function(col_names) {
          nums <- as.numeric(gsub(".*_(\\d)_.*", "\\1", col_names))
          col_names[order(nums)]
        }
        
        dmso_cols <- get_replicate_num(data_cols[dmso_idx])
        tg6h_cols <- get_replicate_num(data_cols[tg6h_idx])
        tg24h_cols <- get_replicate_num(data_cols[tg24h_idx])
        
        cell_TPMs <- cell_TPMs[, c("Geneid", dmso_cols, tg6h_cols, tg24h_cols)]
      }
      
      # Save to file
      output_file <- file.path(output_path, paste0(cell_line, "_TPMs.txt"))
      write.table(cell_TPMs, file = output_file, sep = '\t', 
                 row.names = FALSE, quote = FALSE)
      
      message(paste("  ✓", cell_line, "TPMs processed"))
    }
  }
  
  return(invisible(TRUE))
}

#' Process rat TPM data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_rat_TPMs <- function(input_path = "input/TPMs", 
                            output_path = "processed_TPMs") {
  
  message("Processing PCCL3 (rat) TPMs...")
  
  rat_file <- file.path(input_path, "rat_ILMN_1541_Kalwat_mRNAseq54_Nov2022_TPM.txt")
  
  if (!file.exists(rat_file)) {
    warning(paste("Rat TPM file not found:", rat_file))
    return(invisible(FALSE))
  }
  
  # Read rat TPMs
  rat_TPMs <- read.table(rat_file, header = FALSE, check.names = FALSE)
  colnames(rat_TPMs) <- as.character(rat_TPMs[1, ])
  rat_TPMs <- rat_TPMs[-1, ]
  
  # Extract PCCL3
  PCCL3_cols <- grep("PCCL3", colnames(rat_TPMs))
  
  if (length(PCCL3_cols) > 0) {
    PCCL3_TPMs <- rat_TPMs[, c(1, PCCL3_cols)]
    
    # Identify and reorder columns
    data_cols <- colnames(PCCL3_TPMs)[-1]
    
    dmso_idx <- grep("DMSO", data_cols)
    tg6h_idx <- grep("Tg6h", data_cols)
    tg24h_idx <- grep("Tg24h", data_cols)
    
    if (length(dmso_idx) == 3 && length(tg6h_idx) == 3 && length(tg24h_idx) == 3) {
      # Sort by replicate number
      get_replicate_num <- function(col_names) {
        nums <- as.numeric(gsub(".*_(\\d)_.*", "\\1", col_names))
        col_names[order(nums)]
      }
      
      dmso_cols <- get_replicate_num(data_cols[dmso_idx])
      tg6h_cols <- get_replicate_num(data_cols[tg6h_idx])
      tg24h_cols <- get_replicate_num(data_cols[tg24h_idx])
      
      PCCL3_TPMs <- PCCL3_TPMs[, c("Geneid", dmso_cols, tg6h_cols, tg24h_cols)]
    }
    
    # Save to file
    output_file <- file.path(output_path, "PCCL3_TPMs.txt")
    write.table(PCCL3_TPMs, file = output_file, sep = '\t', 
               row.names = FALSE, quote = FALSE)
    
    message("  ✓ PCCL3 TPMs processed")
  }
  
  return(invisible(TRUE))
}

#' Process human TPM data
#' @param input_path Path to input directory
#' @param output_path Path to output directory
process_human_TPMs <- function(input_path = "input/TPMs", 
                              output_path = "processed_TPMs") {
  
  message("Processing QGP1 (human) TPMs...")
  
  human_file <- file.path(input_path, "human_ILMN_1541_Kalwat_mRNAseq54_Nov2022_TPM.txt")
  
  if (!file.exists(human_file)) {
    warning(paste("Human TPM file not found:", human_file))
    return(invisible(FALSE))
  }
  
  # Read human TPMs
  human_TPMs <- read.table(human_file, header = TRUE, check.names = FALSE)
  
  # Identify QGP1 columns
  qgp1_cols <- grep("QGP1", colnames(human_TPMs))
  
  if (length(qgp1_cols) > 0) {
    QGP1_TPMs <- human_TPMs[, c(1, qgp1_cols)]
    
    # Reorder columns
    data_cols <- colnames(QGP1_TPMs)[-1]
    
    dmso_idx <- grep("DMSO", data_cols)
    tg6h_idx <- grep("Tg6h", data_cols)
    tg24h_idx <- grep("Tg24h", data_cols)
    
    if (length(dmso_idx) == 3 && length(tg6h_idx) == 3 && length(tg24h_idx) == 3) {
      get_replicate_num <- function(col_names) {
        nums <- as.numeric(gsub(".*_(\\d)_.*", "\\1", col_names))
        col_names[order(nums)]
      }
      
      dmso_cols <- get_replicate_num(data_cols[dmso_idx])
      tg6h_cols <- get_replicate_num(data_cols[tg6h_idx])
      tg24h_cols <- get_replicate_num(data_cols[tg24h_idx])
      
      QGP1_TPMs <- QGP1_TPMs[, c("Geneid", dmso_cols, tg6h_cols, tg24h_cols)]
    }
  } else {
    # Fallback to position-based if needed
    if (ncol(human_TPMs) >= 10) {
      QGP1_TPMs <- human_TPMs[, c(1, 2, 5, 8, 3, 6, 9, 4, 7, 10)]
    }
  }
  
  # Save to file
  output_file <- file.path(output_path, "QGP1_TPMs.txt")
  write.table(QGP1_TPMs, file = output_file, sep = '\t', 
             row.names = FALSE, quote = FALSE)
  
  message("  ✓ QGP1 TPMs processed")
  
  return(invisible(TRUE))
}

################################################################################
### SECTION 3: MAIN PREPROCESSING WORKFLOW
################################################################################

#' Run complete preprocessing pipeline
#' @param save_workspace Whether to save R workspace
run_preprocessing <- function(save_workspace = TRUE) {
  
  message("\n========================================")
  message("Starting Data Preprocessing Pipeline")
  message("========================================\n")
  
  # Process read counts
  message("Step 1: Processing read counts...")
  process_mouse_reads()
  process_MIN6_reads()
  process_rat_reads()
  process_human_reads()
  
  # Process TPMs
  message("\nStep 2: Processing TPM data...")
  process_MIN6_TPMs()
  process_mouse_TPMs()
  process_rat_TPMs()
  process_human_TPMs()
  
  # Save workspace if requested
  if (save_workspace) {
    save.image("processed_reads/processed_reads.RData")
    message("\n✓ Workspace saved to processed_reads/processed_reads.RData")
  }
  
  message("\n========================================")
  message("Data Preprocessing Complete!")
  message("========================================")
  
  # Session info
  message("\n=== Session Information ===")
  sessionInfo()
}

################################################################################
### SECTION 4: EXECUTE THE PREPROCESSING
################################################################################

if (interactive()) {
  # Run the preprocessing pipeline
  run_preprocessing(save_workspace = TRUE)
  
  # List output files
  message("\n=== Output Files Created ===")
  
  if (dir.exists("processed_reads")) {
    reads_files <- list.files("processed_reads", pattern = "\\.txt$")
    message("\nProcessed read files:")
    for (f in reads_files) {
      # Check the file to report number of samples
      df <- read.table(file.path("processed_reads", f), header = TRUE, nrows = 1)
      n_samples <- ncol(df) - 1  # Subtract 1 for Geneid column
      message(paste("  -", f, paste0("(", n_samples, " samples)")))
    }
  }
  
  if (dir.exists("processed_TPMs")) {
    tpm_files <- list.files("processed_TPMs", pattern = "\\.txt$")
    message("\nProcessed TPM files:")
    for (f in tpm_files) {
      message(paste("  -", f))
    }
  }
}

writeLines(capture.output(sessionInfo()), "session_info/00_data_preprocessing_sessioninfo.txt")
