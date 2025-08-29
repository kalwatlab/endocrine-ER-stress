################################################################################
### Script 06: Data Preparation for Shiny Deployment
################################################################################

library(tidyverse)
library(jsonlite)

#' Prepare all data for enhanced Shiny deployment
#' @param input_file Path to the merged data file
#' @param output_dir Directory for Shiny app
#' @param include_gsea Include GSEA results if available
#' @param include_upset Include upset plot data if available
#' @param include_unique Include unique gene sets if available
prepare_enhanced_shiny_data <- function(
  input_file = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt",
  output_dir = "RNAseq_Shiny_App",
  include_gsea = TRUE,
  include_upset = TRUE,
  include_unique = TRUE
) {
  
  message("========================================")
  message("Preparing Enhanced Data for Shiny Deployment")
  message("========================================\n")
  
  # Create output directory structure
  dir.create(output_dir, showWarnings = FALSE)
  dir.create(file.path(output_dir, "merged edgeR and TPMs"), 
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "data"), 
             showWarnings = FALSE, recursive = TRUE)
  
  # 1. Process main data (as before)
  message("Loading and optimizing main data...")
  data <- read.delim(input_file, stringsAsFactors = FALSE)
  
  # Remove empty columns
  empty_cols <- sapply(data, function(x) all(is.na(x) | x == ""))
  data <- data[, !empty_cols]
  
  # Reduce precision
  numeric_cols <- sapply(data, is.numeric)
  for (col in names(data)[numeric_cols]) {
    if (!grepl("logFC|FDR|PValue", col)) {
      data[[col]] <- round(data[[col]], 3)
    }
  }
  
  # Save main data
  output_file <- file.path(output_dir, "merged edgeR and TPMs", 
                          "merged_edgeR_TPMs_all_lines.txt")
  write.table(data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("  Main data saved:", round(file.size(output_file)/1024^2, 2), "MB"))
  
  # 2. Copy GSEA results if available
  if (include_gsea) {
    message("\nProcessing GSEA results...")
    
    gsea_dirs <- c("analysis/GSEA_6h", "analysis/GSEA_24h")
    gsea_data_list <- list()
    
    for (gsea_dir in gsea_dirs) {
      if (dir.exists(gsea_dir)) {
        time_point <- ifelse(grepl("6h", gsea_dir), "6h", "24h")
        
        # Look for RDS files with all results
        rds_files <- list.files(gsea_dir, pattern = "all_GSEA_results.rds", 
                               recursive = TRUE, full.names = TRUE)
        
        for (rds_file in rds_files) {
          tryCatch({
            gsea_results <- readRDS(rds_file)
            
            # Skip unique gene results if they exist (shouldn't after script 05 changes)
            if (grepl("unique", rds_file)) {
              message(paste("  Skipping unique gene GSEA file (not applicable):", basename(rds_file)))
              next
            }
            
            # Only process all DEGs
            gene_type <- "all"
            
            # Process and combine results
            for (cell_line in names(gsea_results)) {
              if (!is.null(gsea_results[[cell_line]])) {
                df <- gsea_results[[cell_line]]
                df$cell_line <- cell_line
                df$time_point <- time_point
                df$gene_type <- gene_type
                
                # Store in list
                key <- paste(time_point, gene_type, cell_line, sep = "_")
                gsea_data_list[[key]] <- df
              }
            }
            
            message(paste("  Loaded GSEA results:", time_point, gene_type))
          }, error = function(e) {
            message(paste("  Warning: Could not load", rds_file))
          })
        }
      }
    }
    
    # Combine and save GSEA results
    if (length(gsea_data_list) > 0) {
      gsea_combined <- do.call(rbind, gsea_data_list)
      
      # Convert list columns to character for saving
      list_cols <- sapply(gsea_combined, is.list)
      for (col in names(gsea_combined)[list_cols]) {
        gsea_combined[[col]] <- vapply(gsea_combined[[col]], 
                                       function(x) {
                                         if (is.null(x)) return("")
                                         paste(x, collapse = ", ")
                                       }, 
                                       character(1))
      }
      
      # Remove rows with NA pathway
      gsea_combined <- gsea_combined %>%
        filter(!is.na(pathway) & pathway != "") %>%
        filter(!is.na(NES)) %>%
        filter(!is.na(padj))
      
      # Save as RDS for Shiny (preserves structure better)
      saveRDS(gsea_combined, file.path(output_dir, "data", "gsea_results.rds"))
      
      # Also save as text for reference
      write.table(gsea_combined, 
                 file.path(output_dir, "data", "gsea_results.txt"),
                 sep = "\t", row.names = FALSE, quote = FALSE)
      
      message(paste("  GSEA results saved:", nrow(gsea_combined), "pathways"))
    }
  }
 
  # 3. Process upset plot data
  if (include_upset) {
    message("\nProcessing upset plot data...")
    
    # Create upset data from the main data
    cell_lines <- c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1")
    upset_data <- list()
    
    for (cell_line in cell_lines) {
      for (time_point in c("6h", "24h")) {
        logfc_col <- paste0("logFC_TgvDMSO_", time_point, "_", cell_line)
        fdr_col <- paste0("FDR_TgvDMSO_", time_point, "_", cell_line)
        
        if (all(c(logfc_col, fdr_col) %in% colnames(data))) {
          # Get significantly changed genes
          sig_genes <- data$genes[
            !is.na(data[[logfc_col]]) & 
            !is.na(data[[fdr_col]]) &
            abs(data[[logfc_col]]) > log2(1.5) &  # 1.5-fold change
            data[[fdr_col]] < 0.05
          ]
          
          key <- paste(cell_line, time_point, sep = "_")
          upset_data[[key]] <- sig_genes
        }
      }
    }
    
    # Save upset data
    saveRDS(upset_data, file.path(output_dir, "data", "upset_data.rds"))
    
    # Create summary
    upset_summary <- data.frame(
      set = names(upset_data),
      n_genes = sapply(upset_data, length)
    )
    write.table(upset_summary, 
               file.path(output_dir, "data", "upset_summary.txt"),
               sep = "\t", row.names = FALSE, quote = FALSE)
    
    message(paste("  Upset data saved:", length(upset_data), "gene sets"))
  }
  
  if (include_unique) {
    message("\nProcessing unique gene sets...")
    
    unique_gene_sets <- list()
    
    # Process both 6h and 24h directories
    for (time_point in c("6h", "24h")) {
      unique_dir <- file.path("analysis/custom_gene_sets", time_point)
      
      if (dir.exists(unique_dir)) {
        message(paste("  Processing", time_point, "unique gene sets..."))
        
        # Option 1: Read the RDS file with all gene sets (recommended)
        rds_file <- file.path(unique_dir, paste0("all_unique_gene_sets_", time_point, ".rds"))
        if (file.exists(rds_file)) {
          time_gene_sets <- readRDS(rds_file)
          
          # Add time point to the set names if not already included
          for (set_name in names(time_gene_sets)) {
            # Check if time point is already in the name
            if (!grepl(time_point, set_name)) {
              new_name <- paste0(set_name, "_", time_point)
            } else {
              new_name <- set_name
            }
            unique_gene_sets[[new_name]] <- time_gene_sets[[set_name]]
          }
          
          message(paste("    Loaded", length(time_gene_sets), "gene sets from", time_point))
        } else {
          # Option 2: Fall back to reading individual text files
          message(paste("    RDS file not found, reading individual files for", time_point))
          
          unique_files <- list.files(unique_dir, pattern = "\\.txt$", full.names = TRUE)
          
          # Exclude the summary file
          unique_files <- unique_files[!grepl("summary", unique_files)]
          
          for (file in unique_files) {
            set_name <- gsub("\\.txt$", "", basename(file))
            genes <- readLines(file)
            genes <- genes[genes != "" & !is.na(genes)]
            
            if (length(genes) > 0) {
              unique_gene_sets[[set_name]] <- genes
            }
          }
          
          message(paste("    Loaded", length(unique_files), "gene set files from", time_point))
        }
      } else {
        message(paste("  Warning: Directory not found:", unique_dir))
      }
    }
    
    if (length(unique_gene_sets) > 0) {
      # Save combined unique gene sets
      saveRDS(unique_gene_sets, 
              file.path(output_dir, "data", "unique_gene_sets.rds"))
      
      # Create summary with better organization
      unique_summary <- data.frame(
        gene_set = names(unique_gene_sets),
        n_genes = sapply(unique_gene_sets, length),
        stringsAsFactors = FALSE
      )
      
      # Extract metadata from gene set names
      unique_summary$cell_line <- gsub("_unique.*", "", unique_summary$gene_set)
      unique_summary$time_point <- ifelse(grepl("_6h", unique_summary$gene_set), "6h", 
                                          ifelse(grepl("_24h", unique_summary$gene_set), "24h", "unknown"))
      unique_summary$direction <- ifelse(grepl("_UP_", unique_summary$gene_set), "UP",
                                         ifelse(grepl("_DOWN_", unique_summary$gene_set), "DOWN",
                                                ifelse(grepl("_combined_", unique_summary$gene_set), "combined",
                                                       ifelse(grepl("_ALL_", unique_summary$gene_set), "ALL", "mixed"))))
      
      # Sort by cell line, time point, and direction
      unique_summary <- unique_summary[order(unique_summary$cell_line, 
                                             unique_summary$time_point, 
                                             unique_summary$direction), ]
      
      write.table(unique_summary, 
                  file.path(output_dir, "data", "unique_gene_sets_summary.txt"),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      message(paste("  Unique gene sets saved:", length(unique_gene_sets), "sets total"))
      message(paste("    6h sets:", sum(grepl("_6h", names(unique_gene_sets)))))
      message(paste("    24h sets:", sum(grepl("_24h", names(unique_gene_sets)))))
    } else {
      message("  Warning: No unique gene sets found")
    }
  }
  
  # 5. Copy volcano plots if available
  volcano_files <- c(
    "analysis/volcano_plots/Volcano_Grid_All_Lines.png",
    "analysis/volcano_plots/Volcano_Grid_All_Lines.pdf"
  )
  
  for (vfile in volcano_files) {
    if (file.exists(vfile)) {
      file.copy(vfile, file.path(output_dir, "data", basename(vfile)), 
               overwrite = TRUE)
      message(paste("  Copied:", basename(vfile)))
    }
  }
  
  # 6. Create app configuration file
  config <- list(
    has_gsea = include_gsea && length(gsea_data_list) > 0,
    gsea_note = "GSEA performed on all DEGs (full ranked gene lists)",
    has_upset = include_upset && length(upset_data) > 0,
    has_unique = include_unique && exists("unique_gene_sets"),
    cell_lines = c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1"),
    time_points = c("6h", "24h"),
    data_version = Sys.Date()
  )
  
  saveRDS(config, file.path(output_dir, "data", "app_config.rds"))
  
  message("\n========================================")
  message("Enhanced Shiny data preparation complete!")
  message(paste("Location:", output_dir))
  message("\nData includes:")
  if (config$has_gsea) message("  ✓ GSEA results")
  if (config$has_upset) message("  ✓ Upset plot data")
  if (config$has_unique) message("  ✓ Unique gene sets")
  message("  ✓ Main expression data")
  message("  ✓ Configuration file")
  message("========================================\n")
  
  return(invisible(config))
}

################################################################################
### RUN THE ENHANCED PREPARATION
################################################################################

if (interactive()) {
  # Run the enhanced data preparation
  config <- prepare_enhanced_shiny_data(
    input_file = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt",
    output_dir = "RNAseq_Shiny_App",
    include_gsea = TRUE,
    include_upset = TRUE,
    include_unique = TRUE
  )
  
  # Check what was included
  print(config)
}

writeLines(capture.output(sessionInfo()), "session_info/06A_prepare_data_for_shiny_sessioninfo.txt")
