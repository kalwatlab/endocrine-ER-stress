################################################################################
### Script 03: Granular Analysis - Gene Expression Dot Plots Across Cell Lines
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
  library(pheatmap)
  library(readxl)
})

################################################################################
### SECTION 1: DATA LOADING AND PREPARATION
################################################################################

#' Load merged data from comparative analysis
#' @param data_path Path to merged data file
#' @return Data frame with merged results
load_merged_data <- function(data_path = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt") {
  
  if (!file.exists(data_path)) {
    stop(paste("Merged data file not found:", data_path, 
               "\nPlease run scripts 01 and 02 first."))
  }
  
  merged_data <- read.delim(data_path, stringsAsFactors = FALSE)
  message(paste("Loaded data:", nrow(merged_data), "genes x", ncol(merged_data), "columns"))
  
  return(merged_data)
}

################################################################################
### SECTION 2: GENE SET DEFINITIONS
################################################################################

#' Define default gene sets for analysis
#' @return List of gene sets
define_gene_sets <- function() {
  
  gene_sets <- list(
    # ER stress markers
    ER_stress = c("Ddit3", "Hspa5", "Atf4", "Atf6", "Xbp1", "Ern1", 
                  "Eif2ak3", "Atf3", "Gadd45a", "Herpud1"),
    
    # UPR pathway
    UPR_pathway = c("Hspa5", "Calr", "Pdia4", "Pdia6", "Dnajb9", "Dnajc3",
                    "Hyou1", "Sdf2l1", "Manf", "Creld2"),
    
    # Apoptosis markers
    Apoptosis = c("Casp3", "Casp7", "Casp8", "Casp9", "Bax", "Bcl2",
                  "Bcl2l1", "Mcl1", "Parp1", "Cycs"),
    
    # Beta cell identity
    Beta_cell = c("Ins1", "Ins2", "Pdx1", "Nkx6-1", "Nkx2-2", "Pax6",
                  "Neurod1", "Mafa", "Mafb", "Foxa2"),
    
    # Alpha cell identity  
    Alpha_cell = c("Gcg", "Arx", "Irx1", "Irx2", "Mafb", "Pou3f4"),
    
    # Cell cycle
    Cell_cycle = c("Ccnd1", "Ccnd2", "Ccne1", "Ccna2", "Ccnb1", "Ccnb2",
                   "Cdk1", "Cdk2", "Cdk4", "Cdk6"),
    
    # Oxidative stress
    Oxidative_stress = c("Sod1", "Sod2", "Cat", "Gpx1", "Gpx4", "Gsr",
                         "Nfe2l2", "Keap1", "Hmox1", "Nqo1"),
    
    # Inflammation
    Inflammation = c("Il1b", "Il6", "Tnf", "Cxcl1", "Cxcl2", "Cxcl10",
                     "Ccl2", "Ccl5", "Nfkb1", "Rela")
  )
  
  return(gene_sets)
}

#' Load custom gene set from file
#' @param file_path Path to file with gene list (one gene per line or Excel)
#' @param set_name Name for the gene set
#' @return Named list with gene set
load_custom_gene_set <- function(file_path, set_name = NULL) {
  
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Determine file type and read accordingly
  if (grepl("\\.xlsx?$", file_path)) {
    # Excel file
    genes <- read_excel(file_path)[[1]]  # Take first column
  } else {
    # Text file
    genes <- readLines(file_path)
  }
  
  # Clean gene names
  genes <- trimws(genes)
  genes <- genes[genes != "" & !is.na(genes)]
  
  # Create named list
  if (is.null(set_name)) {
    set_name <- tools::file_path_sans_ext(basename(file_path))
  }
  
  gene_set <- list()
  gene_set[[set_name]] <- genes
  
  message(paste("Loaded", length(genes), "genes for set:", set_name))
  
  return(gene_set)
}

################################################################################
### SECTION 3: DATA EXTRACTION FUNCTIONS
################################################################################

#' Extract expression and statistics for gene sets
#' @param merged_data Merged data frame
#' @param gene_sets List of gene sets
#' @param time_point "6h" or "24h"
#' @param use_log10 Whether to log10 transform TPMs
#' @return Data frame ready for plotting
prepare_dotplot_data <- function(merged_data, gene_sets, 
                                 time_point = "6h", use_log10 = TRUE) {
  
  cell_lines <- c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1")
  plot_data <- list()
  
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    
    for (gene in genes) {
      # Find gene in data (case-insensitive)
      gene_rows <- which(tolower(merged_data$genes) == tolower(gene))
      
      if (length(gene_rows) == 0) {
        message(paste("  Gene not found:", gene))
        next
      }
      
      gene_data <- merged_data[gene_rows[1], ]
      
      for (cell_line in cell_lines) {
        # Extract TPM values (average across replicates)
        tpm_cols <- grep(paste0(cell_line, "_(DMSO|Tg)_"), colnames(gene_data), value = TRUE)
        
        # Separate by condition
        dmso_cols <- grep("DMSO", tpm_cols, value = TRUE)
        tg_cols <- grep(paste0("Tg_", time_point), tpm_cols, value = TRUE)
        
        # Calculate averages
        dmso_tpm <- 0
        tg_tpm <- 0
        
        if (length(dmso_cols) > 0) {
          dmso_vals <- as.numeric(gene_data[, dmso_cols])
          dmso_vals <- dmso_vals[!is.na(dmso_vals)]
          if (length(dmso_vals) > 0) dmso_tpm <- mean(dmso_vals)
        }
        
        if (length(tg_cols) > 0) {
          tg_vals <- as.numeric(gene_data[, tg_cols])
          tg_vals <- tg_vals[!is.na(tg_vals)]
          if (length(tg_vals) > 0) tg_tpm <- mean(tg_vals)
        }
        
        # Average TPM across conditions
        avg_tpm <- mean(c(dmso_tpm, tg_tpm), na.rm = TRUE)
        
        # Get p-value
        pval_col <- paste0("PValue_TgvDMSO_", time_point, "_", cell_line)
        fdr_col <- paste0("FDR_TgvDMSO_", time_point, "_", cell_line)
        logfc_col <- paste0("logFC_TgvDMSO_", time_point, "_", cell_line)
        
        pval <- NA
        fdr <- NA
        logfc <- NA
        
        if (pval_col %in% colnames(gene_data)) {
          pval <- as.numeric(gene_data[[pval_col]])
        }
        if (fdr_col %in% colnames(gene_data)) {
          fdr <- as.numeric(gene_data[[fdr_col]])
        }
        if (logfc_col %in% colnames(gene_data)) {
          logfc <- as.numeric(gene_data[[logfc_col]])
        }
        
        # Store data
        plot_data[[length(plot_data) + 1]] <- data.frame(
          gene = gene_data$genes,
          gene_set = set_name,
          cell_line = cell_line,
          avg_tpm = avg_tpm,
          dmso_tpm = dmso_tpm,
          tg_tpm = tg_tpm,
          pvalue = pval,
          fdr = fdr,
          logfc = logfc,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine all data
  if (length(plot_data) > 0) {
    plot_df <- do.call(rbind, plot_data)
    
    # Transform TPMs if requested
    if (use_log10) {
      plot_df$expression <- log10(plot_df$avg_tpm + 1)  # Add 1 to avoid log(0)
      plot_df$expression_label <- "log10(TPM + 1)"
    } else {
      plot_df$expression <- plot_df$avg_tpm
      plot_df$expression_label <- "TPM"
    }
    
    # Calculate -log10(p-value) for coloring
    plot_df$neg_log10_pval <- -log10(pmax(plot_df$pvalue, 1e-300))  # Avoid Inf
    plot_df$neg_log10_fdr <- -log10(pmax(plot_df$fdr, 1e-300))
    
    # Add significance categories
    plot_df$significance <- case_when(
      plot_df$fdr < 0.001 ~ "***",
      plot_df$fdr < 0.01 ~ "**",
      plot_df$fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
    
    # Add regulation direction
    plot_df$regulation <- case_when(
      plot_df$logfc > 0 & plot_df$fdr < 0.05 ~ "Up",
      plot_df$logfc < 0 & plot_df$fdr < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
    
    return(plot_df)
  } else {
    stop("No data found for the specified gene sets")
  }
}

################################################################################
### SECTION 4: PLOTTING FUNCTIONS
################################################################################

#' Create dot plot for gene expression
#' Added black borders to dots to ensure visibility on white background
#' @param plot_data Data frame from prepare_dotplot_data
#' @param color_by "pvalue", "fdr", "logfc", or "regulation"
#' @param facet_by "gene_set" or NULL
#' @param title Plot title
#' @return ggplot object
create_dotplot <- function(plot_data, color_by = "fdr", 
                           facet_by = "gene_set", title = NULL) {
  
  # Set color scale based on choice
  if (color_by == "pvalue") {
    plot_data$color_value <- plot_data$neg_log10_pval
    color_label <- "-log10(p-value)"
    
    # Use shape 21 (filled circle with border) for better visibility
    p <- ggplot(plot_data, aes(x = cell_line, y = gene)) +
      geom_point(aes(size = expression, fill = color_value), 
                 shape = 21, stroke = 0.5, color = "black") +
      scale_fill_viridis(name = color_label, option = "plasma", 
                         direction = -1, limits = c(0, NA)) +
      scale_size_continuous(name = unique(plot_data$expression_label),
                            range = c(0.5, 8))
    
  } else if (color_by == "fdr") {
    plot_data$color_value <- plot_data$neg_log10_fdr
    color_label <- "-log10(FDR)"
    
    p <- ggplot(plot_data, aes(x = cell_line, y = gene)) +
      geom_point(aes(size = expression, fill = color_value), 
                 shape = 21, stroke = 0.5, color = "black") +
      scale_fill_viridis(name = color_label, option = "plasma", 
                         direction = -1, limits = c(0, NA)) +
      scale_size_continuous(name = unique(plot_data$expression_label),
                            range = c(0.5, 8))
    
  } else if (color_by == "logfc") {
    plot_data$color_value <- plot_data$logfc
    color_label <- "log2(FC)"
    max_fc <- max(abs(plot_data$logfc), na.rm = TRUE)
    
    p <- ggplot(plot_data, aes(x = cell_line, y = gene)) +
      geom_point(aes(size = expression, fill = color_value), 
                 shape = 21, stroke = 0.5, color = "black") +
      scale_fill_gradient2(name = color_label, 
                           low = "blue", mid = "white", high = "red",
                           midpoint = 0, limits = c(-max_fc, max_fc)) +
      scale_size_continuous(name = unique(plot_data$expression_label),
                            range = c(0.5, 8))
    
  } else if (color_by == "regulation") {
    plot_data$color_value <- plot_data$regulation
    color_label <- "Regulation"
    
    p <- ggplot(plot_data, aes(x = cell_line, y = gene)) +
      geom_point(aes(size = expression, fill = color_value), 
                 shape = 21, stroke = 0.5, color = "black") +
      scale_fill_manual(name = color_label,
                        values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
      scale_size_continuous(name = unique(plot_data$expression_label),
                            range = c(0.5, 8))
  }
  
  # Add common theme elements
  p <- p +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "black")
    ) +
    labs(
      x = "Cell Line",
      y = "Gene",
      title = title
    )
  
  # Add faceting if requested
  if (!is.null(facet_by) && facet_by %in% colnames(plot_data)) {
    p <- p + facet_wrap(~ get(facet_by), scales = "free_y", ncol = 1)
  }
  
  # Add significance markers
  if (any(plot_data$significance != "", na.rm = TRUE)) {
    sig_data <- plot_data[plot_data$significance != "" & !is.na(plot_data$significance), ]
    p <- p + geom_text(data = sig_data, aes(label = significance), 
                       size = 3, vjust = 0.5, hjust = -0.2)
  }
  
  return(p)
}

#' Create heatmap-style plot
#' @param plot_data Data frame from prepare_dotplot_data
#' @param value_type "expression", "logfc", or "tg_tpm"
#' @param title Plot title
#' @return ggplot object
create_heatmap_plot <- function(plot_data, value_type = "logfc", title = NULL) {
  
  # Prepare matrix for heatmap
  if (value_type == "expression") {
    value_col <- "expression"
    value_label <- unique(plot_data$expression_label)
  } else if (value_type == "logfc") {
    value_col <- "logfc"
    value_label <- "log2(Fold Change)"
  } else if (value_type == "tg_tpm") {
    value_col <- "tg_tpm"
    value_label <- "Tg TPM"
  }
  
  # Reshape data to wide format
  heatmap_data <- plot_data %>%
    dplyr::select(gene, cell_line, all_of(value_col)) %>%
    pivot_wider(names_from = cell_line, values_from = all_of(value_col))
  
  # Convert to matrix
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  rownames(heatmap_matrix) <- heatmap_data$gene
  
  # Handle missing values
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Create heatmap
  if (value_type == "logfc") {
    # Use diverging color palette for fold changes
    colors <- colorRampPalette(c("blue", "white", "red"))(100)
  } else {
    # Use sequential palette for expression
    colors <- colorRampPalette(c("white", "orange", "red"))(100)
  }
  
  pheatmap(heatmap_matrix,
           color = colors,
           main = title,
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 10,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           border_color = "grey60",
           cellwidth = 30,
           angle_col = 45)
}

################################################################################
### SECTION 5: MULTI-PANEL PLOTS WITH PLOT TYPE CONTROL
################################################################################

#' Create comprehensive multi-panel plot with control over plot types
#' @param merged_data Merged data frame
#' @param gene_sets List of gene sets
#' @param output_path Path to save plots
#' @param plot_types Character vector specifying which plot types to create.
#'                   Options: "log2FC", "regulation", "FDR", or "all" (default)
#' @return List of plots
create_comprehensive_plots <- function(merged_data, gene_sets, 
                                       output_path = "analysis/dotplots",
                                       plot_types = "all") {
  
  # Validate and process plot_types parameter
  valid_types <- c("log2FC", "regulation", "FDR", "all")
  
  if ("all" %in% plot_types) {
    plot_types <- c("log2FC", "regulation", "FDR")
  } else {
    # Convert to lowercase for case-insensitive matching
    plot_types_lower <- tolower(plot_types)
    valid_types_lower <- tolower(valid_types)
    
    # Map the input types to standard names
    plot_types_mapped <- character(0)
    if ("log2fc" %in% plot_types_lower || "logfc" %in% plot_types_lower) {
      plot_types_mapped <- c(plot_types_mapped, "log2FC")
    }
    if ("regulation" %in% plot_types_lower) {
      plot_types_mapped <- c(plot_types_mapped, "regulation")
    }
    if ("fdr" %in% plot_types_lower) {
      plot_types_mapped <- c(plot_types_mapped, "FDR")
    }
    
    plot_types <- plot_types_mapped
    
    if (length(plot_types) == 0) {
      message("No valid plot types specified. Defaulting to 'log2FC'")
      plot_types <- "log2FC"
    }
  }
  
  message(paste("Plot types to generate:", paste(plot_types, collapse = ", ")))
  
  # Create output directory
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Create subdirectories for each plot type
  if ("log2FC" %in% plot_types) {
    dir.create(file.path(output_path, "log2FC_plots"), showWarnings = FALSE)
  }
  if ("regulation" %in% plot_types) {
    dir.create(file.path(output_path, "regulation_plots"), showWarnings = FALSE)
  }
  if ("FDR" %in% plot_types) {
    dir.create(file.path(output_path, "FDR_plots"), showWarnings = FALSE)
  }
  
  plots <- list()
  
  # Generate plots for both time points
  for (time_point in c("6h", "24h")) {
    
    message(paste("\nGenerating plots for", time_point, "time point..."))
    
    # Prepare data
    plot_data <- prepare_dotplot_data(merged_data, gene_sets, 
                                      time_point = time_point, use_log10 = TRUE)
    
    # Generate only requested plot types
    
    # 1. Dot plot colored by FDR
    if ("FDR" %in% plot_types) {
      p1 <- create_dotplot(plot_data, color_by = "fdr", facet_by = "gene_set",
                           title = paste("Gene Expression -", time_point, "Tg vs DMSO (FDR colored)"))
      plots[[paste0("dotplot_fdr_", time_point)]] <- p1
      
      # Save plot
      ggsave(file.path(output_path, "FDR_plots", paste0("dotplot_fdr_", time_point, ".pdf")),
             plot = p1, width = 12, height = length(unique(plot_data$gene)) * 0.3 + 3,
             dpi = 300, limitsize = FALSE)
      
      # Create separate plots for each gene set
      for (set_name in unique(plot_data$gene_set)) {
        set_data <- plot_data[plot_data$gene_set == set_name, ]
        
        p_set <- create_dotplot(set_data, color_by = "fdr", facet_by = NULL,
                                title = paste(set_name, "-", time_point, "Tg vs DMSO (FDR)"))
        
        plots[[paste0(set_name, "_fdr_", time_point)]] <- p_set
        
        clean_name <- gsub("[^A-Za-z0-9_-]", "_", set_name)
        ggsave(file.path(output_path, "FDR_plots", paste0(clean_name, "_fdr_", time_point, ".pdf")),
               plot = p_set, width = 10, height = length(unique(set_data$gene)) * 0.4 + 2,
               dpi = 300, limitsize = FALSE)
      }
    }
    
    # 2. Dot plot colored by log fold change
    if ("log2FC" %in% plot_types) {
      p2 <- create_dotplot(plot_data, color_by = "logfc", facet_by = "gene_set",
                           title = paste("Gene Expression -", time_point, "Tg vs DMSO (log2FC colored)"))
      plots[[paste0("dotplot_logfc_", time_point)]] <- p2
      
      ggsave(file.path(output_path, "log2FC_plots", paste0("dotplot_logfc_", time_point, ".pdf")),
             plot = p2, width = 12, height = length(unique(plot_data$gene)) * 0.3 + 3,
             dpi = 300, limitsize = FALSE)
      
      # Create separate plots for each gene set
      for (set_name in unique(plot_data$gene_set)) {
        set_data <- plot_data[plot_data$gene_set == set_name, ]
        
        p_set <- create_dotplot(set_data, color_by = "logfc", facet_by = NULL,
                                title = paste(set_name, "-", time_point, "Tg vs DMSO (log2FC)"))
        
        plots[[paste0(set_name, "_logfc_", time_point)]] <- p_set
        
        clean_name <- gsub("[^A-Za-z0-9_-]", "_", set_name)
        ggsave(file.path(output_path, "log2FC_plots", paste0(clean_name, "_logfc_", time_point, ".pdf")),
               plot = p_set, width = 10, height = length(unique(set_data$gene)) * 0.4 + 2,
               dpi = 300, limitsize = FALSE)
      }
    }
    
    # 3. Dot plot colored by regulation
    if ("regulation" %in% plot_types) {
      p3 <- create_dotplot(plot_data, color_by = "regulation", facet_by = "gene_set",
                           title = paste("Gene Regulation -", time_point, "Tg vs DMSO"))
      plots[[paste0("dotplot_regulation_", time_point)]] <- p3
      
      ggsave(file.path(output_path, "regulation_plots", paste0("dotplot_regulation_", time_point, ".pdf")),
             plot = p3, width = 12, height = length(unique(plot_data$gene)) * 0.3 + 3,
             dpi = 300, limitsize = FALSE)
      
      # Create separate plots for each gene set
      for (set_name in unique(plot_data$gene_set)) {
        set_data <- plot_data[plot_data$gene_set == set_name, ]
        
        p_set <- create_dotplot(set_data, color_by = "regulation", facet_by = NULL,
                                title = paste(set_name, "-", time_point, "Tg vs DMSO (Regulation)"))
        
        plots[[paste0(set_name, "_regulation_", time_point)]] <- p_set
        
        clean_name <- gsub("[^A-Za-z0-9_-]", "_", set_name)
        ggsave(file.path(output_path, "regulation_plots", paste0(clean_name, "_regulation_", time_point, ".pdf")),
               plot = p_set, width = 10, height = length(unique(set_data$gene)) * 0.4 + 2,
               dpi = 300, limitsize = FALSE)
      }
    }
  }
  
  message(paste("\nPlots saved to:", output_path))
  message(paste("  Generated", length(plots), "plots"))
  
  return(plots)
}

################################################################################
### SECTION 6: SUMMARY STATISTICS
################################################################################

#' Generate summary statistics for gene sets
#' @param merged_data Merged data frame
#' @param gene_sets List of gene sets
#' @param output_path Path to save summary
#' @return Data frame with summary statistics
generate_summary_stats <- function(merged_data, gene_sets,
                                   output_path = "analysis/dotplots") {
  
  cell_lines <- c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1")
  summary_list <- list()
  
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    
    for (time_point in c("6h", "24h")) {
      for (cell_line in cell_lines) {
        
        # Count significant genes
        sig_genes <- 0
        up_genes <- 0
        down_genes <- 0
        
        for (gene in genes) {
          gene_rows <- which(tolower(merged_data$genes) == tolower(gene))
          
          if (length(gene_rows) > 0) {
            gene_data <- merged_data[gene_rows[1], ]
            
            fdr_col <- paste0("FDR_TgvDMSO_", time_point, "_", cell_line)
            logfc_col <- paste0("logFC_TgvDMSO_", time_point, "_", cell_line)
            
            if (fdr_col %in% colnames(gene_data) && logfc_col %in% colnames(gene_data)) {
              fdr <- as.numeric(gene_data[[fdr_col]])
              logfc <- as.numeric(gene_data[[logfc_col]])
              
              if (!is.na(fdr) && !is.na(logfc) && fdr < 0.05) {
                sig_genes <- sig_genes + 1
                if (logfc > 0) up_genes <- up_genes + 1
                if (logfc < 0) down_genes <- down_genes + 1
              }
            }
          }
        }
        
        summary_list[[length(summary_list) + 1]] <- data.frame(
          gene_set = set_name,
          time_point = time_point,
          cell_line = cell_line,
          total_genes = length(genes),
          sig_genes = sig_genes,
          up_genes = up_genes,
          down_genes = down_genes,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine summary
  summary_df <- do.call(rbind, summary_list)
  
  # Save summary
  write.table(summary_df,
              file = file.path(output_path, "gene_set_summary_statistics.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("\nSummary statistics saved")
  
  return(summary_df)
}

################################################################################
### SECTION 7: MAIN WORKFLOW WITH PLOT TYPE CONTROL
################################################################################

#' Run complete dot plot analysis with control over plot types
#' @param gene_sets List of gene sets or "default" to use predefined sets
#' @param custom_gene_files List of paths to custom gene set files
#' @param output_path Path to save all outputs
#' @param plot_types Character vector specifying which plot types to create. 
#'                   Options: "log2FC", "regulation", "FDR", or "all" (default)
#'                   Can specify multiple: c("log2FC", "FDR")
#' @return List of results
run_dotplot_analysis <- function(gene_sets = "default",
                                 custom_gene_files = NULL,
                                 output_path = "analysis/dotplots",
                                 plot_types = "all") {
  
  message("\n========================================")
  message("Starting Dot Plot Analysis")
  message("========================================\n")
  
  # Validate plot_types parameter
  valid_types <- c("log2FC", "regulation", "FDR", "all")
  if ("all" %in% plot_types) {
    plot_types <- c("log2FC", "regulation", "FDR")
  } else {
    # Convert to lowercase for case-insensitive matching
    plot_types_lower <- tolower(plot_types)
    
    # Map the input types to standard names
    plot_types_mapped <- character(0)
    if ("log2fc" %in% plot_types_lower || "logfc" %in% plot_types_lower) {
      plot_types_mapped <- c(plot_types_mapped, "log2FC")
    }
    if ("regulation" %in% plot_types_lower) {
      plot_types_mapped <- c(plot_types_mapped, "regulation")
    }
    if ("fdr" %in% plot_types_lower) {
      plot_types_mapped <- c(plot_types_mapped, "FDR")
    }
    
    plot_types <- plot_types_mapped
    
    if (length(plot_types) == 0) {
      message("No valid plot types specified. Defaulting to 'log2FC'")
      plot_types <- "log2FC"
    }
  }
  
  message(paste("Plot types to generate:", paste(plot_types, collapse = ", ")))
  
  # Load merged data
  message("Loading merged data...")
  merged_data <- load_merged_data()
  
  # Prepare gene sets
  if (is.character(gene_sets) && gene_sets == "default") {
    message("\nUsing default gene sets...")
    gene_sets <- define_gene_sets()
  }
  
  # Add custom gene sets if provided
  if (!is.null(custom_gene_files)) {
    message("\nLoading custom gene sets...")
    for (file in custom_gene_files) {
      custom_set <- load_custom_gene_set(file)
      gene_sets <- c(gene_sets, custom_set)
    }
  }
  
  message(paste("\nTotal gene sets:", length(gene_sets)))
  for (set_name in names(gene_sets)) {
    message(paste("  ", set_name, ":", length(gene_sets[[set_name]]), "genes"))
  }
  
  # Generate plots with specified types only
  message("\nGenerating plots...")
  plots <- create_comprehensive_plots(merged_data, gene_sets, output_path, plot_types)
  
  # Generate summary statistics
  message("\nGenerating summary statistics...")
  summary_stats <- generate_summary_stats(merged_data, gene_sets, output_path)
  
  # Print summary
  message("\n=== Summary ===")
  summary_by_set <- summary_stats %>%
    group_by(gene_set) %>%
    summarise(
      total_sig_6h = sum(sig_genes[time_point == "6h"]),
      total_sig_24h = sum(sig_genes[time_point == "24h"]),
      .groups = 'drop'
    )
  print(summary_by_set)
  
  message("\n========================================")
  message("Dot Plot Analysis Complete!")
  message(paste("Results saved to:", output_path))
  message(paste("Plot types generated:", paste(plot_types, collapse = ", ")))
  message("========================================\n")
  
  return(list(
    plots = plots,
    summary = summary_stats,
    gene_sets = gene_sets,
    plot_types = plot_types
  ))
}

################################################################################
### SECTION 8: EXECUTE ANALYSIS
################################################################################

if (interactive()) {
  # Default run with all plot types
  # results <- run_dotplot_analysis(gene_sets = "default")
  
  # Example: Run with only log2FC plots
  results <- run_dotplot_analysis(gene_sets = "default", plot_types = "log2FC")
  
  # Example: Run with multiple specific plot types
  # results <- run_dotplot_analysis(gene_sets = "default", plot_types = c("log2FC", "FDR"))
  
  message("To run analysis, uncomment one of the examples above or call run_dotplot_analysis() directly")
}

writeLines(capture.output(sessionInfo()), "session_info/03_granular_analysis_dotplots_sessioninfo.txt")