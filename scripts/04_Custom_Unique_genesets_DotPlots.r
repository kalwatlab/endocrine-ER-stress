################################################################################
### Script 04: Create Custom Gene Sets and dot plots from Unique Genes per Cell Line
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# run_dotplot_analysis() function is in script 03. If scripts are run in order, the function will be loaded and will work.

################################################################################
### Section 1: FUNCTION TO CREATE TOP REGULATED GENE SETS PER TIME POINT
################################################################################

#' Create gene sets from top unique genes per cell line and time point
#' @param input_path Path to unique gene files from script 02
#' @param output_path Path to save gene set files
#' @param n_top Number of top up/down genes to select (default 16)
#' @param time_point Which time point to analyze ("6h" or "24h")
#' @return List of gene sets ready for dot plot analysis
create_top_unique_gene_sets_timepoint <- function(input_path = "edgeR_output",
                                                  output_path = "analysis/custom_gene_sets",
                                                  n_top = 16,
                                                  time_point = "6h") {
  
  message("\n========================================")
  message(paste("Creating Custom Gene Sets from Unique Genes -", time_point))
  message("========================================\n")
  
  # Validate time point
  if (!time_point %in% c("6h", "24h")) {
    stop("time_point must be either '6h' or '24h'")
  }
  
  # Create output directory with time point subdirectory
  output_dir <- file.path(output_path, time_point)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  cell_lines <- c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1")
  
  # Initialize gene sets
  gene_sets_up <- list()
  gene_sets_down <- list()
  gene_sets_combined <- list()
  gene_sets_all <- list()
  
  for (cell_line in cell_lines) {
    # Read time-point-specific unique genes file
    unique_file <- file.path(input_path, paste0("Unique_genes_", cell_line, "_", time_point, ".txt"))
    
    if (!file.exists(unique_file)) {
      message(paste("  Warning: File not found for", cell_line, "at", time_point))
      next
    }
    
    # Load data
    unique_data <- read.delim(unique_file, stringsAsFactors = FALSE)
    
    if (nrow(unique_data) == 0) {
      message(paste("  No unique genes for", cell_line, "at", time_point))
      next
    }
    
    message(paste("\nProcessing", cell_line, "-", nrow(unique_data), "unique genes at", time_point))
    
    # Use specific time point columns
    logfc_col <- paste0("logFC_TgvDMSO_", time_point, "_", cell_line)
    fdr_col <- paste0("FDR_TgvDMSO_", time_point, "_", cell_line)
    
    # Check if columns exist
    if (!logfc_col %in% colnames(unique_data)) {
      message(paste("  LogFC column not found:", logfc_col))
      next
    }
    
    # Convert to numeric and handle NAs
    unique_data[[logfc_col]] <- as.numeric(unique_data[[logfc_col]])
    unique_data <- unique_data[!is.na(unique_data[[logfc_col]]), ]
    
    # Add FDR filter if column exists
    if (fdr_col %in% colnames(unique_data)) {
      unique_data[[fdr_col]] <- as.numeric(unique_data[[fdr_col]])
      # Filter for significant genes only (should already be significant if in unique edgeR output)
      unique_data <- unique_data[!is.na(unique_data[[fdr_col]]) & 
                                   unique_data[[fdr_col]] < 0.05, ]
    }
    
    if (nrow(unique_data) == 0) {
      message(paste("  No significant unique genes for", cell_line, "at", time_point))
      next
    }
    
    # Sort by absolute log fold change
    unique_data <- unique_data[order(abs(unique_data[[logfc_col]]), decreasing = TRUE), ]
    
    # Get top upregulated genes
    up_genes <- unique_data[unique_data[[logfc_col]] > 0, ]
    top_up <- character(0)  # Initialize as empty
    if (nrow(up_genes) > 0) {
      top_up <- head(up_genes$genes, n_top)
      gene_set_name <- paste0(cell_line, "_unique_UP_", time_point)
      gene_sets_up[[gene_set_name]] <- top_up
      
      message(paste("  Top", length(top_up), "upregulated genes selected at", time_point))
      
      # Save to file
      writeLines(top_up, file.path(output_dir, 
                                   paste0(cell_line, "_unique_top", n_top, "_UP_", time_point, ".txt")))
    }
    
    # Get top downregulated genes
    down_genes <- unique_data[unique_data[[logfc_col]] < 0, ]
    top_down <- character(0)  # Initialize as empty
    if (nrow(down_genes) > 0) {
      top_down <- head(down_genes$genes, n_top)
      gene_set_name <- paste0(cell_line, "_unique_DOWN_", time_point)
      gene_sets_down[[gene_set_name]] <- top_down
      
      message(paste("  Top", length(top_down), "downregulated genes selected at", time_point))
      
      # Save to file
      writeLines(top_down, file.path(output_dir, 
                                     paste0(cell_line, "_unique_top", n_top, "_DOWN_", time_point, ".txt")))
    }
    
    # Combined top regulated (both up AND down) - CORRECTED
    # Combine the top up and top down genes instead of just taking top by absolute value
    top_combined <- c(top_up, top_down)
    
    if (length(top_combined) > 0) {
      gene_set_name <- paste0(cell_line, "_unique_top", length(top_combined), "_", time_point)
      gene_sets_combined[[gene_set_name]] <- top_combined
      
      message(paste("  Combined:", length(top_up), "up +", length(top_down), "down =", 
                    length(top_combined), "total genes"))
      
      # Save combined to file
      writeLines(top_combined, file.path(output_dir, 
                                         paste0(cell_line, "_unique_top", length(top_combined), "_combined_", time_point, ".txt")))
    }
    
    # Also create a set with ALL unique genes (for reference)
    gene_set_name <- paste0(cell_line, "_unique_ALL_", time_point)
    gene_sets_all[[gene_set_name]] <- unique_data$genes
  }
  
  # Combine all gene sets
  all_gene_sets <- c(gene_sets_up, gene_sets_down, gene_sets_combined)
  
  # Save a master file with all gene sets for this time point
  saveRDS(all_gene_sets, file.path(output_dir, paste0("all_unique_gene_sets_", time_point, ".rds")))
  
  # Create summary
  summary_df <- data.frame(
    gene_set = names(all_gene_sets),
    n_genes = sapply(all_gene_sets, length),
    time_point = time_point,
    stringsAsFactors = FALSE
  )
  
  write.table(summary_df, 
              file.path(output_dir, paste0("gene_sets_summary_", time_point, ".txt")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("\n=== Summary for", time_point, "===")
  print(summary_df)
  
  message(paste("\nGene sets saved to:", output_dir))
  
  return(all_gene_sets)
}

################################################################################
### Section 2: WRAPPER FUNCTION TO PROCESS TIME POINTS
################################################################################

#' Process unique gene sets for specified time points
#' @param input_path Path to unique gene files
#' @param output_path Path to save gene set files
#' @param n_top Number of top genes to select
#' @param time_points Vector of time points to analyze (default c("6h", "24h"))
#' @return List of all gene sets organized by time point
process_all_timepoints <- function(input_path = "edgeR_output",
                                   output_path = "analysis/custom_gene_sets",
                                   n_top = 16,
                                   time_points = c("6h", "24h")) {
  
  all_results <- list()
  
  for (tp in time_points) {
    message("\n========================================")
    message(paste("Processing", tp, "time point"))
    message("========================================")
    
    gene_sets <- create_top_unique_gene_sets_timepoint(
      input_path = input_path,
      output_path = output_path,
      n_top = n_top,
      time_point = tp
    )
    
    all_results[[tp]] <- gene_sets
  }
  
  return(all_results)
}

################################################################################
### Section 3: RUN THE ANALYSIS AND CREATE DOT PLOTS
################################################################################

if (interactive()) {
  
  # Step 1: Create the custom gene sets from unique genes for all time points
  message("\n=== STEP 1: Creating custom gene sets for all time points ===")
  
  all_gene_sets <- process_all_timepoints(
    input_path = "edgeR_output",
    output_path = "analysis/custom_gene_sets",
    n_top = 16  # Will select top 16 up and top 16 down
  )
  
  # Step 2: Load the dot plot script if not already loaded
  if (!exists("run_dotplot_analysis")) {
    message("\n=== STEP 2: Loading dot plot analysis script ===")
    source("03_Granular_Analysis_DotPlots.r")
  }
  
  # Step 3: Run dot plot analysis for each time point
  # Note: We only want log2FC plots, not regulation or FDR plots
  message("\n=== STEP 3: Creating dot plots for unique genes (log2FC only) ===")
  
  # Check if run_dotplot_analysis supports plot_types parameter
  # If not, this would need to be added to script 03
  
  # Process 6h time point
  if ("6h" %in% names(all_gene_sets)) {
    message("\n--- Creating dot plots for 6h unique genes ---")
    dotplot_results_6h <- run_dotplot_analysis(
      gene_sets = all_gene_sets[["6h"]],
      output_path = "analysis/unique_gene_dotplots/6h",
      plot_types = "log2FC"  # Only create log2FC plots
    )
  }
  
  # Process 24h time point
  if ("24h" %in% names(all_gene_sets)) {
    message("\n--- Creating dot plots for 24h unique genes ---")
    dotplot_results_24h <- run_dotplot_analysis(
      gene_sets = all_gene_sets[["24h"]],
      output_path = "analysis/unique_gene_dotplots/24h",
      plot_types = "log2FC"  # Only create log2FC plots
    )
  }
  
  message("\n========================================")
  message("Analysis Complete!")
  message("Dot plots saved to:")
  message("  - analysis/unique_gene_dotplots/6h")
  message("  - analysis/unique_gene_dotplots/24h")
  message("\nOnly log2FC plots have been generated as requested.")
  message("========================================")
}

writeLines(capture.output(sessionInfo()), "session_info/04_custom_unique_genesets_dotplots_sessioninfo.txt")