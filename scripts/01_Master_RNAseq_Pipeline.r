################################################################################
### Script 01: Master RNAseq Analysis Pipeline for Thapsigargin Treatment Study
################################################################################

# Clear workspace
rm(list = ls())

# Load all required libraries
suppressPackageStartupMessages({
  library(edgeR)        # v4.0+ for differential expression
  library(limma)        # for linear modeling
  library(tidyverse)    # for data manipulation
  library(ggplot2)      # for plotting
  library(ggrepel)      # for label repelling in plots
  library(orthogene)    # for ortholog conversion
  library(clusterProfiler) # for GO analysis
  library(AnnotationDbi)
  library(org.Mm.eg.db) # Mouse annotation
  library(org.Rn.eg.db) # Rat annotation
  library(org.Hs.eg.db) # Human annotation
  library(readxl)       # for Excel files
  library(janitor)      # for cleaning names
  library(patchwork)    # for combining plots
})

# Set global parameters
CELL_LINES <- list(
  mouse = c("MIN6", "aTC1", "MGN3", "GLUTag"),
  rat = c("PCCL3"),
  human = c("QGP1")
)

# Define experimental design (common for all cell lines)
EXPERIMENTAL_DESIGN <- data.frame(
  replicate = rep(paste0("N=", 1:3), 3),
  treatment = rep(c("DMSO", "Tg", "Tg"), each = 3),
  time = c(rep("24h", 3), rep("6h", 3), rep("24h", 3))
)

# Fold change and p-value cutoffs
FC_CUTOFF_HIGH <- log2(2)
FC_CUTOFF_LOW <- log2(0.5)
PVAL_CUTOFF <- 0.05

################################################################################
### SECTION 1: DATA IMPORT AND PREPROCESSING FUNCTIONS
################################################################################

#' Import and organize read count data for a specific cell line
#' @param cell_line Character string specifying the cell line
#' @param data_path Path to the processed reads file
#' @return Matrix of read counts with proper column names
import_reads <- function(cell_line, data_path = "processed_reads") {
  
  file_path <- file.path(data_path, paste0(cell_line, "_genes.txt"))
  
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Read data
  data <- read.table(file_path, sep = '\t', header = TRUE, row.names = 1)
  
  # Standardize column names
  colnames(data) <- c("DMSO.24h.1", "DMSO.24h.2", "DMSO.24h.3",
                       "Tg.6h.1", "Tg.6h.2", "Tg.6h.3",
                       "Tg.24h.1", "Tg.24h.2", "Tg.24h.3")
  
  return(as.matrix(data))
}

#' Import TPM data for a specific cell line
#' @param cell_line Character string specifying the cell line
#' @param data_path Path to the processed TPMs file
#' @return Data frame of TPM values
import_tpms <- function(cell_line, data_path = "processed_TPMs") {
  
  file_path <- file.path(data_path, paste0(cell_line, "_TPMs.txt"))
  
  if (!file.exists(file_path)) {
    warning(paste("TPM file not found:", file_path))
    return(NULL)
  }
  
  return(read.table(file_path, sep = '\t', header = TRUE))
}

################################################################################
### SECTION 2: EDGER ANALYSIS FUNCTION
################################################################################

#' Perform edgeR analysis for a single cell line
#' @param cell_line Character string specifying the cell line
#' @param reads_path Path to reads directory
#' @param tpm_path Path to TPM directory
#' @param output_path Path to output directory
#' @return List containing edgeR results and merged data
run_edgeR_analysis <- function(cell_line, 
                               reads_path = "processed_reads",
                               tpm_path = "processed_TPMs",
                               output_path = "edgeR_output") {
  
  message(paste("\n=== Running edgeR analysis for", cell_line, "==="))
  
  # Import read counts
  read_counts <- import_reads(cell_line, reads_path)
  
  # Create DGEList object
  group_factor <- factor(
    c(rep("DMSO.24h", 3), rep("Tg.6h", 3), rep("Tg.24h", 3)),
    levels = c("DMSO.24h", "Tg.6h", "Tg.24h")
  )
  
  dge_list <- DGEList(
    counts = read_counts,
    genes = rownames(read_counts),
    group = group_factor
  )
  
  # Filter lowly expressed genes
  keep <- filterByExpr(dge_list)
  dge_list <- dge_list[keep, , keep.lib.sizes = FALSE]
  message(paste("  Kept", sum(keep), "genes after filtering"))
  
  # Normalize
  dge_list <- calcNormFactors(dge_list)
  
  # Design matrix
  design <- model.matrix(~0 + group, data = dge_list$samples)
  colnames(design) <- levels(dge_list$samples$group)
  
  # Estimate dispersion
  dge_list <- estimateDisp(dge_list, design)
  
  # Fit model
  fit <- glmQLFit(dge_list, design)
  
  # Define contrasts
  contrasts <- makeContrasts(
    TgvDMSO.6h = Tg.6h - DMSO.24h,
    TgvDMSO.24h = Tg.24h - DMSO.24h,
    Tg.24hv6h = Tg.24h - Tg.6h,
    levels = design
  )
  
  # Perform tests
  results <- list()
  for (contrast_name in colnames(contrasts)) {
    qlf_test <- glmQLFTest(fit, contrast = contrasts[, contrast_name])
    results[[contrast_name]] <- as.data.frame(topTags(qlf_test, n = Inf))
    
    # Add cell line suffix to column names
    cols_to_rename <- 2:6
    colnames(results[[contrast_name]])[cols_to_rename] <- 
      paste(colnames(results[[contrast_name]])[cols_to_rename], 
            gsub("\\.", "_", contrast_name), cell_line, sep = "_")
  }
  
  # Merge results
  merged_results <- purrr::reduce(results, full_join, by = "genes")
  
  # Import and merge TPMs if available
  tpm_data <- import_tpms(cell_line, tpm_path)
  if (!is.null(tpm_data)) {
    # Handle different column names for gene IDs
    gene_col <- ifelse("Geneid" %in% colnames(tpm_data), "Geneid", "gene_symbol")
    merged_with_tpm <- full_join(merged_results, tpm_data, 
                                  by = c("genes" = gene_col))
  } else {
    merged_with_tpm <- merged_results
  }
  
  # Save outputs
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Full results
  write.table(merged_results, 
              file = file.path(output_path, paste0(cell_line, "_edgeR_output.txt")),
              sep = '\t', row.names = FALSE, quote = FALSE)
  
  # Brief results (logFC and FDR only)
  brief_cols <- c(1, grep("logFC|FDR", colnames(merged_results)))
  brief_results <- merged_results[, brief_cols]
  write.table(brief_results,
              file = file.path(output_path, paste0(cell_line, "_edgeR_brief_output.txt")),
              sep = '\t', row.names = FALSE, quote = FALSE)
  
  # Results with TPM
  if (!is.null(tpm_data)) {
    write.table(merged_with_tpm,
                file = file.path(output_path, paste0(cell_line, "_edgeR_TPM_output.txt")),
                sep = '\t', row.names = FALSE, quote = FALSE)
  }
  
  return(list(
    results = merged_results,
    brief = brief_results,
    with_tpm = merged_with_tpm,
    dge_list = dge_list
  ))
}

################################################################################
### SECTION 3: VOLCANO PLOT FUNCTIONS
################################################################################

#' Create volcano plot for a specific comparison (NO LEGEND VERSION)
#' @param data Data frame with logFC and FDR columns
#' @param logfc_col Name of log fold change column
#' @param fdr_col Name of FDR column
#' @param title Plot title
#' @param cell_line Cell line name for labeling
#' @return ggplot object
create_volcano_plot <- function(data, logfc_col, fdr_col, title, cell_line) {
  
  # Add significance categories
  plot_data <- data %>%
    mutate(
      Change = case_when(
        .data[[logfc_col]] > FC_CUTOFF_HIGH & .data[[fdr_col]] < PVAL_CUTOFF ~ "Increased",
        .data[[logfc_col]] < FC_CUTOFF_LOW & .data[[fdr_col]] < PVAL_CUTOFF ~ "Decreased",
        TRUE ~ "Unchanged"
      ),
      neg_log10_fdr = -log10(.data[[fdr_col]])
    )
  
  # Select top genes for labeling
  top_genes <- plot_data %>%
    filter(Change == "Increased") %>%
    mutate(distance = neg_log10_fdr + abs(.data[[logfc_col]])) %>%
    top_n(10, distance)  # Reduced from 20 to 10 for cleaner plots
  
  bottom_genes <- plot_data %>%
    filter(Change == "Decreased") %>%
    mutate(distance = neg_log10_fdr + abs(.data[[logfc_col]])) %>%
    top_n(10, distance)  # Reduced from 20 to 10
  
  label_genes <- rbind(top_genes, bottom_genes)
  
  # Define colors
  colors <- c("Decreased" = "steelblue1", "Increased" = "tomato1", "Unchanged" = "grey70")
  
  # Create plot WITHOUT legend
  p <- ggplot(plot_data, aes(x = .data[[logfc_col]], y = neg_log10_fdr)) +
    geom_point(aes(color = Change), alpha = 0.6, size = 1.5) +
    geom_point(data = label_genes, 
               aes(fill = Change), 
               shape = 21, color = "black", size = 2) +
    geom_text_repel(data = label_genes,
                    aes(label = genes),
                    size = 2.5,  # Smaller text for grid
                    box.padding = 0.3,
                    max.overlaps = 20) +
    scale_color_manual(values = colors, guide = "none") +  # Remove legend
    scale_fill_manual(values = colors, guide = "none") +   # Remove legend
    geom_vline(xintercept = c(FC_CUTOFF_LOW, FC_CUTOFF_HIGH), 
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = -log10(PVAL_CUTOFF), 
               linetype = "dashed", color = "red", alpha = 0.5) +
    labs(
      title = paste(cell_line, "-", title),
      x = expression('Log'[2]*' Fold Change'),
      y = expression('-Log'[10]*' FDR')
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    ) +
    xlim(-8, 10)
  
  return(p)
}

#' Generate all volcano plots and create grid layout
#' @param edger_results Results from run_edgeR_analysis
#' @param cell_line Cell line name
#' @param output_path Path to save plots
#' @return List of ggplot objects
generate_volcano_plots <- function(edger_results, cell_line, 
                                   output_path = "analysis/volcano_plots") {
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  data <- edger_results$results
  plots <- list()
  
  # Find logFC and FDR columns
  comparisons <- list(
    "6h" = list(logfc = paste0("logFC_TgvDMSO_6h_", cell_line),
                fdr = paste0("FDR_TgvDMSO_6h_", cell_line),
                title = "6h Tg vs DMSO"),
    "24h" = list(logfc = paste0("logFC_TgvDMSO_24h_", cell_line),
                 fdr = paste0("FDR_TgvDMSO_24h_", cell_line),
                 title = "24h Tg vs DMSO")
  )
  
  for (comp_name in names(comparisons)) {
    comp <- comparisons[[comp_name]]
    
    if (comp$logfc %in% colnames(data) && comp$fdr %in% colnames(data)) {
      p <- create_volcano_plot(data, comp$logfc, comp$fdr, comp$title, cell_line)
      plots[[comp_name]] <- p
      
      # Save individual plot
      ggsave(
        filename = file.path(output_path, paste0("Volcano_", cell_line, "_", comp_name, ".svg")),
        plot = p,
        width = 6,
        height = 5,
        dpi = 300
      )
    }
  }
  
  return(plots)
}

#' Create combined volcano plot grid for all cell lines
#' @param all_volcano_plots List of volcano plots for all cell lines
#' @param output_path Path to save the combined plot
#' @return Combined plot
create_volcano_grid <- function(all_volcano_plots, 
                                output_path = "analysis/volcano_plots") {
  
  library(patchwork)
  
  # Extract plots in the specified order
  # Row 1: MIN6, aTC1
  # Row 2: QGP1, GLUTag
  # Row 3: MGN3, PCCL3
  
  plot_list <- list()
  
  # Collect plots in order
  cell_order <- c("MIN6", "aTC1", "QGP1", "GLUTag", "MGN3", "PCCL3")
  
  for (cell_line in cell_order) {
    if (cell_line %in% names(all_volcano_plots)) {
      # Add 6h plot
      if ("6h" %in% names(all_volcano_plots[[cell_line]])) {
        plot_list[[paste0(cell_line, "_6h")]] <- all_volcano_plots[[cell_line]][["6h"]]
      } else {
        # Create empty plot if missing
        plot_list[[paste0(cell_line, "_6h")]] <- ggplot() + 
          theme_void() + 
          labs(title = paste(cell_line, "- 6h (No data)"))
      }
      
      # Add 24h plot
      if ("24h" %in% names(all_volcano_plots[[cell_line]])) {
        plot_list[[paste0(cell_line, "_24h")]] <- all_volcano_plots[[cell_line]][["24h"]]
      } else {
        plot_list[[paste0(cell_line, "_24h")]] <- ggplot() + 
          theme_void() + 
          labs(title = paste(cell_line, "- 24h (No data)"))
      }
    }
  }
  
  # Create the grid layout
  # Each row has 2 cell lines × 2 time points = 4 plots
  combined_plot <- (
    (plot_list[["MIN6_6h"]] | plot_list[["MIN6_24h"]] | 
       plot_list[["aTC1_6h"]] | plot_list[["aTC1_24h"]]) /
      (plot_list[["QGP1_6h"]] | plot_list[["QGP1_24h"]] | 
         plot_list[["GLUTag_6h"]] | plot_list[["GLUTag_24h"]]) /
      (plot_list[["MGN3_6h"]] | plot_list[["MGN3_24h"]] | 
         plot_list[["PCCL3_6h"]] | plot_list[["PCCL3_24h"]])
  ) + 
    plot_annotation(
      title = "Volcano Plots: Thapsigargin vs DMSO Treatment",
      subtitle = "All cell lines at 6h and 24h time points",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  # Save the combined plot
  ggsave(
    filename = file.path(output_path, "Volcano_Grid_All_Lines.pdf"),
    plot = combined_plot,
    width = 16,
    height = 12,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_path, "Volcano_Grid_All_Lines.png"),
    plot = combined_plot,
    width = 16,
    height = 12,
    dpi = 300
  )
  
  message("✓ Volcano plot grid saved")
  
  return(combined_plot)
}

################################################################################
### SECTION 4: GO ANALYSIS FUNCTIONS
################################################################################

#' Perform GO enrichment analysis
#' @param genes Vector of gene symbols
#' @param species Species for annotation ("mouse", "rat", or "human")
#' @param direction Direction of change ("UP" or "DOWN")
#' @param ont GO ontology ("BP", "MF", or "CC")
#' @return enrichResult object
perform_go_analysis <- function(genes, species = "mouse", direction = "ALL", 
                                ont = "BP") {
  
  # Select appropriate OrgDb
  org_db <- switch(species,
                   mouse = org.Mm.eg.db,
                   rat = org.Rn.eg.db,
                   human = org.Hs.eg.db,
                   org.Mm.eg.db)  # Default to mouse
  
  # Run enrichGO
  go_result <- enrichGO(
    gene = genes,
    OrgDb = org_db,
    keyType = "SYMBOL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  return(go_result)
}

################################################################################
### SECTION 5: MAIN ANALYSIS WORKFLOW
################################################################################

#' Run complete analysis pipeline for all cell lines
run_complete_pipeline <- function() {
  
  message("\n========================================")
  message("Starting Complete RNAseq Analysis Pipeline")
  message("========================================\n")
  
  # Store results for all cell lines
  all_results <- list()
  all_volcano_plots <- list()
  
  # Process each cell line
  all_cell_lines <- unlist(CELL_LINES)
  
  for (cell_line in all_cell_lines) {
    
    tryCatch({
      # Run edgeR analysis
      edger_results <- run_edgeR_analysis(cell_line)
      all_results[[cell_line]] <- edger_results
      
      # Generate volcano plots
      volcano_plots <- generate_volcano_plots(edger_results, cell_line)
      all_volcano_plots[[cell_line]] <- volcano_plots
      
      message(paste("✓ Completed analysis for", cell_line))
      
    }, error = function(e) {
      message(paste("✗ Error processing", cell_line, ":", e$message))
    })
  }
  
  # Create the combined grid after all individual plots are generated
  if (length(all_volcano_plots) > 0) {
    message("\nCreating volcano plot grid...")
    volcano_grid <- create_volcano_grid(all_volcano_plots)
  }
  
  return(all_results)

    # Save workspace
  save.image("edgeR_output/Master_RNAseq_Pipeline.RData")
  message("\n✓ Workspace saved to edgeR_output/Master_RNAseq_Pipeline.RData")
  
}

################################################################################
### SECTION 6: RUN THE PIPELINE
################################################################################

# Execute the pipeline
if (interactive()) {
  # Run the complete pipeline
  results <- run_complete_pipeline()
  
  # Print summary
  message("\n========================================")
  message("Pipeline Execution Complete")
  message(paste("Successfully processed:", 
                paste(names(results), collapse = ", ")))
  message("========================================\n")
  
  # Session information
  message("\n=== Session Information ===")
  sessionInfo()
  writeLines(capture.output(sessionInfo()), "session_info/01_master_RNAseq_pipeline_sessioninfo.txt")
  
}
