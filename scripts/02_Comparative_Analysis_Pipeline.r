################################################################################
### Script 02: Comparative Analysis Pipeline - Cross-Species Gene Comparison
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(orthogene)
  library(UpSetR)
  library(VennDiagram)
  library(ggVennDiagram)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Hs.eg.db)
})

# Figure parameters for publication
FIGURE_CONFIG <- list(
  upset_text_scale = 1.3,      # Increase for larger upset plot text
  volcano_label_size = 3,       # Gene label size in volcano plots  
  dotplot_gene_size = 4,        # Gene name size in dot plots
  axis_text_size = 12,          # Axis text size
  title_size = 14,               # Plot title size
  base_text_scale = c(1.5, 1.3, 1, 1, 1.3, 1), # font sizes in upset plots
  base_text_scale_combined = c(1.3, 1.1, 0.9, 0.9, 1.1, 0.9) # font sizes in combined upset plot
)

# Source the main pipeline if needed
#source("scripts/01_Master_RNAseq_Pipeline.R")

################################################################################
### SECTION 1: ORTHOLOG CONVERSION FUNCTIONS
################################################################################

#' Convert gene names to mouse orthologs - COMPREHENSIVE VERSION
#' Tries multiple methods and keeps all genes, even unmapped ones
#' @param gene_df Data frame with gene names and other data
#' @param input_species Source species ("human" or "rat")
#' @param gene_column Name of column containing gene symbols
#' @param biomart_file Optional local file with BiomaRt mappings (if server issues)
#' @param verbose Print detailed progress messages
#' @return Data frame with mouse orthologs and all original data columns

convert_to_mouse_orthologs <- function(gene_df, 
                                       input_species, 
                                       gene_column = "genes",
                                       biomart_file = NULL,
                                       verbose = TRUE) {
  
  if (verbose) {
    message(paste("\n=== Converting", input_species, "genes to mouse orthologs (COMPREHENSIVE) ==="))
    message(paste("  Starting with", nrow(gene_df), "genes"))
  }
  
  # Store original data and gene order
  original_df <- gene_df
  original_genes <- gene_df[[gene_column]]
  
  # Initialize tracking dataframe
  mapping_results <- data.frame(
    original_gene = original_genes,
    mouse_gene = NA_character_,
    mapping_method = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Remove any duplicate genes in input
  mapping_results <- mapping_results[!duplicated(mapping_results$original_gene), ]
  
  # Track unmapped genes
  unmapped_genes <- original_genes
  
  # ============================================================================
  # STEP 1: Direct Symbol Matching (case-insensitive)
  # ============================================================================
  if (verbose) message("\n  Step 1: Direct symbol matching...")
  
  # For rat-to-mouse, many genes have identical symbols except for case
  if (input_species == "rat") {
    # Try exact match first
    exact_matches <- mapping_results$original_gene
    mapping_results$mouse_gene[!is.na(exact_matches)] <- exact_matches[!is.na(exact_matches)]
    mapping_results$mapping_method[!is.na(mapping_results$mouse_gene)] <- "direct_exact"
    
    # Try case conversion (Rat genes often just need first letter capitalized)
    still_unmapped <- is.na(mapping_results$mouse_gene)
    if (any(still_unmapped)) {
      case_converted <- stringr::str_to_title(mapping_results$original_gene[still_unmapped])
      mapping_results$mouse_gene[still_unmapped] <- case_converted
      mapping_results$mapping_method[still_unmapped] <- "direct_case"
    }
    
    mapped_count <- sum(!is.na(mapping_results$mouse_gene))
    if (verbose) message(paste("    Mapped", mapped_count, "genes by direct matching"))
    unmapped_genes <- mapping_results$original_gene[is.na(mapping_results$mouse_gene)]
  }
  
  # For human-to-mouse, try case variations
  if (input_species == "human") {
    # HUMAN genes are often ALL CAPS, mouse are Title Case
    case_converted <- stringr::str_to_title(mapping_results$original_gene)
    mapping_results$mouse_gene <- case_converted
    mapping_results$mapping_method <- "direct_case"
    
    # Keep track of what we've tentatively mapped
    # We'll verify these with ortholog databases
    tentative_mappings <- mapping_results
  }
  
  # ============================================================================
  # STEP 2: Manual Mappings for Known Problem Genes
  # ============================================================================
  if (verbose) message("\n  Step 2: Manual mappings for known genes...")
  
  # Define manual mappings for genes known to have issues
  manual_mappings <- data.frame(
    human = c("NEAT1", "MALAT1", "XIST", "TSIX", "GAS5", "SNHG1", "TUG1", 
              "H19", "MEG3", "UCA1", "HOTAIR", "KCNQ1OT1", "AIRN",
              "INS", "GCG", "SST", "PPY", "GHRL", "IAPP"),
    mouse = c("Neat1", "Malat1", "Xist", "Tsix", "Gas5", "Snhg1", "Tug1",
              "H19", "Meg3", "Uca1", "Hotair", "Kcnq1ot1", "Airn",
              "Ins1", "Gcg", "Sst", "Ppy", "Ghrl", "Iapp"),
    rat = c("Neat1", "Malat1", "Xist", "Tsix", "Gas5", "Snhg1", "Tug1",
            "H19", "Meg3", "Uca1", "Hotair", "Kcnq1ot1", "Airn",
            "Ins1", "Gcg", "Sst", "Ppy", "Ghrl", "Iapp"),
    stringsAsFactors = FALSE
  )
  
  # Note: INS maps to Ins1 in mouse/rat (they also have Ins2)
  
  # Apply manual mappings
  source_col <- if (input_species == "human") "human" else "rat"
  
  for (i in 1:nrow(manual_mappings)) {
    source_gene <- manual_mappings[[source_col]][i]
    target_gene <- manual_mappings[["mouse"]][i]
    
    # Find this gene in our unmapped list
    gene_idx <- which(mapping_results$original_gene == source_gene & is.na(mapping_results$mouse_gene))
    
    if (length(gene_idx) > 0) {
      mapping_results$mouse_gene[gene_idx] <- target_gene
      mapping_results$mapping_method[gene_idx] <- "manual"
      if (verbose) message(paste("    Manually mapped:", source_gene, "->", target_gene))
    }
  }
  
  unmapped_genes <- mapping_results$original_gene[is.na(mapping_results$mouse_gene)]
  
  # ============================================================================
  # STEP 3: Orthogene with Multiple Methods
  # ============================================================================
  if (length(unmapped_genes) > 0) {
    if (verbose) message(paste("\n  Step 3: Orthogene mapping for", length(unmapped_genes), "remaining genes..."))
    
    # Try different orthogene methods
    orthogene_methods <- c("gprofiler", "homologene", "babelgene")
    
    for (method in orthogene_methods) {
      if (length(unmapped_genes) == 0) break
      
      if (verbose) message(paste("    Trying", method, "..."))
      
      tryCatch({
        # Create a dataframe with just unmapped genes
        unmapped_df <- data.frame(genes = unmapped_genes, stringsAsFactors = FALSE)
        
        # Try orthogene with this method
        converted <- orthogene::convert_orthologs(
          gene_df = unmapped_df,
          gene_input = "genes",
          gene_output = "columns",  # Get both input and output
          input_species = input_species,
          output_species = "mouse",
          non121_strategy = "kp_popular",  # Keep popular in many-to-many
          drop_nonorths = FALSE,  # Don't drop non-orthologous genes
          method = method
        )
        
        # Process results
        if (!is.null(converted) && nrow(converted) > 0) {
          if ("ortholog_gene" %in% colnames(converted)) {
            # Update our mapping results
            for (i in 1:nrow(converted)) {
              orig_gene <- converted$input_gene[i]
              if (is.na(orig_gene)) orig_gene <- converted$genes[i]
              mouse_gene <- converted$ortholog_gene[i]
              
              if (!is.na(orig_gene) && !is.na(mouse_gene)) {
                idx <- which(mapping_results$original_gene == orig_gene & is.na(mapping_results$mouse_gene))
                if (length(idx) > 0) {
                  mapping_results$mouse_gene[idx[1]] <- mouse_gene
                  mapping_results$mapping_method[idx[1]] <- paste0("orthogene_", method)
                }
              }
            }
            
            newly_mapped <- sum(!is.na(mapping_results$mouse_gene)) - 
              sum(!is.na(mapping_results$mouse_gene) | mapping_results$mapping_method != paste0("orthogene_", method))
            if (verbose) message(paste("      Mapped", nrow(converted), "genes with", method))
          }
        }
        
      }, error = function(e) {
        if (verbose) message(paste("      Warning:", method, "failed:", e$message))
      })
      
      # Update unmapped list
      unmapped_genes <- mapping_results$original_gene[is.na(mapping_results$mouse_gene)]
    }
  }
  
  # ============================================================================
  # STEP 4: BiomaRt or Local File
  # ============================================================================
  unmapped_genes <- mapping_results$original_gene[is.na(mapping_results$mouse_gene)]
  
  if (length(unmapped_genes) > 0 && !is.null(biomart_file)) {
    if (verbose) message(paste("\n  Step 4: Using local BiomaRt file for", length(unmapped_genes), "remaining genes..."))
    
    if (file.exists(biomart_file)) {
      # Read the local mapping file
      # Expected format: columns for human_symbol, mouse_symbol, rat_symbol
      biomart_mappings <- read.csv(biomart_file, stringsAsFactors = FALSE)
      
      source_col <- if (input_species == "human") "human_symbol" else "rat_symbol"
      
      if (source_col %in% colnames(biomart_mappings) && "mouse_symbol" %in% colnames(biomart_mappings)) {
        for (gene in unmapped_genes) {
          matches <- biomart_mappings[biomart_mappings[[source_col]] == gene, "mouse_symbol"]
          matches <- unique(matches[!is.na(matches) & matches != ""])
          
          if (length(matches) > 0) {
            # Take first match or most common if multiple
            idx <- which(mapping_results$original_gene == gene)
            mapping_results$mouse_gene[idx] <- matches[1]
            mapping_results$mapping_method[idx] <- "biomart_local"
          }
        }
        
        newly_mapped <- sum(mapping_results$mapping_method == "biomart_local", na.rm = TRUE)
        if (verbose) message(paste("    Mapped", newly_mapped, "genes from local BiomaRt file"))
      }
    } else {
      if (verbose) message(paste("    Warning: BiomaRt file not found:", biomart_file))
    }
  } else if (length(unmapped_genes) > 0) {
    # Try live BiomaRt connection
    if (verbose) message(paste("\n  Step 4: Trying BiomaRt for", length(unmapped_genes), "remaining genes..."))
    
    tryCatch({
      library(biomaRt)
      
      # Set up marts
      if (input_species == "human") {
        source_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        source_attr <- "hgnc_symbol"
      } else {  # rat
        source_mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
        source_attr <- "rgd_symbol"
      }
      mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      
      # Get orthologs
      orthologs <- getLDS(
        attributes = c(source_attr),
        filters = source_attr,
        values = unmapped_genes,
        mart = source_mart,
        attributesL = c("mgi_symbol"),
        martL = mouse_mart
      )
      
      if (nrow(orthologs) > 0) {
        colnames(orthologs) <- c("source_gene", "mouse_gene")
        
        for (i in 1:nrow(orthologs)) {
          idx <- which(mapping_results$original_gene == orthologs$source_gene[i] & 
                         is.na(mapping_results$mouse_gene))
          if (length(idx) > 0) {
            mapping_results$mouse_gene[idx[1]] <- orthologs$mouse_gene[i]
            mapping_results$mapping_method[idx[1]] <- "biomart_live"
          }
        }
        
        if (verbose) message(paste("    Mapped", nrow(orthologs), "genes via BiomaRt"))
      }
      
    }, error = function(e) {
      if (verbose) message(paste("    Warning: BiomaRt connection failed:", e$message))
      if (verbose) message("    Consider providing a local BiomaRt file with biomart_file parameter")
    })
  }
  
  # ============================================================================
  # STEP 5: Keep Unmapped Genes with Original Names
  # ============================================================================
  unmapped_genes <- mapping_results$original_gene[is.na(mapping_results$mouse_gene)]
  
  if (length(unmapped_genes) > 0) {
    if (verbose) message(paste("\n  Step 5: Keeping", length(unmapped_genes), "unmapped genes with original names"))
    
    # For unmapped genes, keep original name
    unmapped_idx <- is.na(mapping_results$mouse_gene)
    mapping_results$mouse_gene[unmapped_idx] <- mapping_results$original_gene[unmapped_idx]
    mapping_results$mapping_method[unmapped_idx] <- "no_ortholog"
    
    # Print some examples if verbose
    if (verbose && length(unmapped_genes) > 0) {
      examples <- head(unmapped_genes, 5)
      message(paste("    Examples:", paste(examples, collapse = ", ")))
    }
  }
  
  # ============================================================================
  # STEP 6: Merge Results Back with Original Data
  # ============================================================================
  if (verbose) message("\n  Step 6: Merging results with original data...")
  
  # Join mapping results with original dataframe
  result_df <- original_df %>%
    dplyr::left_join(mapping_results, by = c("genes" = "original_gene")) %>%
    dplyr::mutate(
      original_gene = genes,
      genes = mouse_gene
    ) %>%
    dplyr::select(-mouse_gene)  # Remove temporary column
  
  # Add mapping statistics column
  result_df$ortholog_confidence <- case_when(
    mapping_results$mapping_method == "direct_exact" ~ "high",
    mapping_results$mapping_method == "manual" ~ "high",
    mapping_results$mapping_method == "direct_case" ~ "high",
    mapping_results$mapping_method %in% c("orthogene_gprofiler", "orthogene_homologene", "orthogene_babelgene") ~ "high",
    mapping_results$mapping_method %in% c("biomart_live", "biomart_local") ~ "medium",
    mapping_results$mapping_method == "no_ortholog" ~ "none",
    TRUE ~ "unknown"
  )
  
  # ============================================================================
  # Summary Report
  # ============================================================================
  if (verbose) {
    message("\n=== Ortholog Conversion Summary ===")
    message(paste("  Input genes:", nrow(original_df)))
    message(paste("  Output genes:", nrow(result_df)))
    
    # Count by method
    method_counts <- table(mapping_results$mapping_method)
    message("\n  Mapping methods used:")
    for (method in names(method_counts)) {
      message(paste("    ", method, ":", method_counts[method], "genes"))
    }
    
    # Check for specific genes of interest
    genes_of_interest <- c("NEAT1", "Neat1", "MALAT1", "Malat1", "XIST", "Xist")
    for (gene in genes_of_interest) {
      if (gene %in% original_df[[gene_column]]) {
        idx <- which(mapping_results$original_gene == gene)
        if (length(idx) > 0) {
          mapped_to <- mapping_results$mouse_gene[idx[1]]
          method <- mapping_results$mapping_method[idx[1]]
          message(paste("  ", gene, "->", mapped_to, "(", method, ")"))
        }
      }
    }
  }
  
  return(result_df)
}

################################################################################
### SECTION 2: DATA MERGING FUNCTIONS
################################################################################

#' Load and process edgeR results for all cell lines
#' @param cell_lines Vector of cell line names
#' @param input_path Path to edgeR output directory
#' @param mouse_symbol_file Optional file for updating mouse symbols
#' @return List of processed data frames
load_all_edger_results <- function(cell_lines = c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1"),
                                   input_path = "edgeR_output",
                                   mouse_symbol_file = NULL) {
  
  results_list <- list()
  
  for (cell_line in cell_lines) {
    file_path <- file.path(input_path, paste0(cell_line, "_edgeR_TPM_output.txt"))
    
    if (file.exists(file_path)) {
      data <- read.delim(file_path)
      
      # Store original dimensions for debugging
      orig_dims <- dim(data)
      
      # Determine species
      species <- case_when(
        cell_line == "QGP1" ~ "human",
        cell_line == "PCCL3" ~ "rat",
        TRUE ~ "mouse"
      )
      
      # Convert non-mouse to mouse orthologs
      if (species != "mouse") {
        set.seed(123)
        data <- convert_to_mouse_orthologs(data, species)
      }
      
      results_list[[cell_line]] <- data
      message(paste("✓ Loaded", cell_line, 
                   "- Original:", orig_dims[1], "x", orig_dims[2],
                   "Final:", nrow(data), "x", ncol(data)))
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  
  return(results_list)
}

#' Merge all cell line results into a simplified data frame
#' SIMPLIFIED: Only keeps genes, logFC, FDR, and TPM columns
#' @param results_list List of data frames from load_all_edger_results
#' @return Merged data frame with reduced columns
merge_all_results <- function(results_list) {
  
  # Process each data frame to select relevant columns
  processed_list <- list()
  
  for (cell_line in names(results_list)) {
    df <- results_list[[cell_line]]
    
    # Debug: Print column names for troubleshooting
    message(paste("  Processing", cell_line, "with", ncol(df), "columns"))
    
    # Start with genes column
    cols_to_keep <- "genes"
    
    # SIMPLIFIED: Only keep logFC and FDR columns (exclude logCPM, F, PValue)
    logfc_6h_col <- paste0("logFC_TgvDMSO_6h_", cell_line)
    logfc_24h_col <- paste0("logFC_TgvDMSO_24h_", cell_line)
    logfc_24v6h_col <- paste0("logFC_Tg_24hv6h_", cell_line)
    
    fdr_6h_col <- paste0("FDR_TgvDMSO_6h_", cell_line)
    fdr_24h_col <- paste0("FDR_TgvDMSO_24h_", cell_line)
    fdr_24v6h_col <- paste0("FDR_Tg_24hv6h_", cell_line)
    
    # Also catch any variations in column naming
    logfc_cols <- grep(paste0("logFC.*", cell_line), colnames(df), value = TRUE)
    fdr_cols <- grep(paste0("FDR.*", cell_line), colnames(df), value = TRUE)
    
    # Explicitly exclude CPM, F, and PValue columns
    # Combine only logFC and FDR columns
    simplified_cols <- unique(c(logfc_6h_col, logfc_24h_col, logfc_24v6h_col,
                                fdr_6h_col, fdr_24h_col, fdr_24v6h_col,
                                logfc_cols, fdr_cols))
    
    # Keep only columns that exist
    simplified_cols <- simplified_cols[simplified_cols %in% colnames(df)]
    cols_to_keep <- c(cols_to_keep, simplified_cols)
    
    # Handle TPM columns
    if (cell_line == "MIN6") {
      # MIN6 has different column naming from the original Excel file
      tpm_cols <- grep("dmso_24h|thapsigargin", colnames(df), value = TRUE)
      if (length(tpm_cols) > 0) {
        # Rename MIN6 TPM columns to match pattern
        for (i in seq_along(tpm_cols)) {
          old_name <- tpm_cols[i]
          new_name <- old_name
          
          # Replace patterns to standardize naming
          new_name <- gsub("dmso_24h", paste0(cell_line, "_DMSO_24h"), new_name)
          new_name <- gsub("thapsigargin_6h", paste0(cell_line, "_Tg_6h"), new_name)
          new_name <- gsub("thapsigargin_24h", paste0(cell_line, "_Tg_24h"), new_name)
          
          colnames(df)[colnames(df) == old_name] <- new_name
          tpm_cols[i] <- new_name
        }
        cols_to_keep <- c(cols_to_keep, tpm_cols)
      }
    } else {
      # For other cell lines, look for TPM columns
      krs_cols <- grep("^KRS", colnames(df), value = TRUE)
      if (length(krs_cols) > 0) {
        # Rename KRS columns to standard format
        new_tpm_cols <- character(length(krs_cols))
        for (i in seq_along(krs_cols)) {
          col <- krs_cols[i]
          # Parse the column name to extract treatment, time, and replicate
          if (grepl("DMSO", col)) {
            rep_num <- gsub(".*_(\\d)_DMSO.*", "\\1", col)
            new_name <- paste0(cell_line, "_DMSO_24h_", rep_num)
          } else if (grepl("Tg6h", col)) {
            rep_num <- gsub(".*_(\\d)_Tg6h.*", "\\1", col)
            new_name <- paste0(cell_line, "_Tg_6h_", rep_num)
          } else if (grepl("Tg24h", col)) {
            rep_num <- gsub(".*_(\\d)_Tg24h.*", "\\1", col)
            new_name <- paste0(cell_line, "_Tg_24h_", rep_num)
          } else {
            new_name <- col
          }
          colnames(df)[colnames(df) == col] <- new_name
          new_tpm_cols[i] <- new_name
        }
        cols_to_keep <- c(cols_to_keep, new_tpm_cols)
      }
      
      # Also look for already formatted TPM columns
      standard_tpm_cols <- grep(paste0("^", cell_line, "_(DMSO|Tg)_"), colnames(df), value = TRUE)
      cols_to_keep <- c(cols_to_keep, standard_tpm_cols)
    }
    
    # Keep unique columns that exist in the dataframe
    cols_to_keep <- unique(cols_to_keep[cols_to_keep %in% colnames(df)])
    
    if (length(cols_to_keep) > 0) {
      processed_list[[cell_line]] <- df[, cols_to_keep, drop = FALSE]
      message(paste("    Kept", length(cols_to_keep), "columns for", cell_line))
      message(paste("      Columns:", paste(head(cols_to_keep, 10), collapse = ", "),
                    ifelse(length(cols_to_keep) > 10, "...", "")))
    } else {
      warning(paste("No columns selected for", cell_line))
    }
  }
  
  # Merge all data frames
  if (length(processed_list) > 0) {
    merged_df <- Reduce(function(x, y) {
      merge(x, y, by = "genes", all = TRUE)
    }, processed_list)
    
    # Replace NA with empty string for character columns only
    for (col in colnames(merged_df)) {
      if (is.character(merged_df[[col]])) {
        merged_df[[col]][is.na(merged_df[[col]])] <- ""
      }
    }
    
    message(paste("  Final merged data:", nrow(merged_df), "genes x", ncol(merged_df), "columns"))
    message("  SIMPLIFIED: Excluded logCPM, F, and PValue columns")
  } else {
    stop("No data frames to merge")
  }
  
  return(merged_df)
}

################################################################################
### SECTION 3: DIFFERENTIAL GENE IDENTIFICATION
################################################################################

#' Identify significantly changed genes for each cell line - this identifies changed genes a 6h and 24h, then combines those two lists
#' @param merged_df Merged data frame from merge_all_results
#' @param fc_cutoff_high Upper fold change cutoff (log2)
#' @param fc_cutoff_low Lower fold change cutoff (log2)
#' @param fdr_cutoff FDR cutoff
#' @param debug Print debugging information
#' @return List of gene sets
identify_changed_genes <- function(merged_df, 
                                  fc_cutoff_high = log2(2),
                                  fc_cutoff_low = log2(0.5),
                                  fdr_cutoff = 0.05,
                                  debug = FALSE) {
  
  cell_lines <- c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1")
  changed_genes <- list()
  
  for (cell_line in cell_lines) {
    # Use the correct column names from edgeR output
    
    # Get genes changed at 6h
    logfc_6h_col <- paste0("logFC_TgvDMSO_6h_", cell_line)
    fdr_6h_col <- paste0("FDR_TgvDMSO_6h_", cell_line)
    
    if (debug) {
      message(paste("\n  Looking for 6h columns for", cell_line, ":"))
      message(paste("    LogFC column:", logfc_6h_col, 
                   "- Found:", logfc_6h_col %in% colnames(merged_df)))
      message(paste("    FDR column:", fdr_6h_col, 
                   "- Found:", fdr_6h_col %in% colnames(merged_df)))
    }
    
    genes_6h <- character(0)
    if (logfc_6h_col %in% colnames(merged_df) && fdr_6h_col %in% colnames(merged_df)) {
      # Convert columns to numeric if they aren't already
      logfc_6h_values <- as.numeric(merged_df[[logfc_6h_col]])
      fdr_6h_values <- as.numeric(merged_df[[fdr_6h_col]])
      
      # Filter for significantly changed genes
      valid_6h <- !is.na(logfc_6h_values) & !is.na(fdr_6h_values) & 
                  merged_df[[fdr_6h_col]] != ""
      
      genes_6h <- merged_df$genes[
        valid_6h &
        (logfc_6h_values > fc_cutoff_high | logfc_6h_values < fc_cutoff_low) &
        fdr_6h_values < fdr_cutoff
      ]
      
      if (debug && length(genes_6h) > 0) {
        message(paste("    Found", length(genes_6h), "significant genes at 6h"))
      }
    }
    
    # Get genes changed at 24h
    logfc_24h_col <- paste0("logFC_TgvDMSO_24h_", cell_line)
    fdr_24h_col <- paste0("FDR_TgvDMSO_24h_", cell_line)
    
    if (debug) {
      message(paste("  Looking for 24h columns for", cell_line, ":"))
      message(paste("    LogFC column:", logfc_24h_col, 
                   "- Found:", logfc_24h_col %in% colnames(merged_df)))
      message(paste("    FDR column:", fdr_24h_col, 
                   "- Found:", fdr_24h_col %in% colnames(merged_df)))
    }
    
    genes_24h <- character(0)
    if (logfc_24h_col %in% colnames(merged_df) && fdr_24h_col %in% colnames(merged_df)) {
      # Convert columns to numeric if they aren't already
      logfc_24h_values <- as.numeric(merged_df[[logfc_24h_col]])
      fdr_24h_values <- as.numeric(merged_df[[fdr_24h_col]])
      
      # Filter for significantly changed genes
      valid_24h <- !is.na(logfc_24h_values) & !is.na(fdr_24h_values) & 
                   merged_df[[fdr_24h_col]] != ""
      
      genes_24h <- merged_df$genes[
        valid_24h &
        (logfc_24h_values > fc_cutoff_high | logfc_24h_values < fc_cutoff_low) &
        fdr_24h_values < fdr_cutoff
      ]
      
      if (debug && length(genes_24h) > 0) {
        message(paste("    Found", length(genes_24h), "significant genes at 24h"))
      }
    }
    
    # Combine 6h and 24h
    changed_genes[[cell_line]] <- unique(c(genes_6h, genes_24h))
    
    message(paste(cell_line, ":", length(changed_genes[[cell_line]]), "changed genes"))
  }
  
  return(changed_genes)
}

#' Identify significantly changed genes for each cell line BY TIME POINT - this keeps the lists separate for better UpSet plot and unique cell line DEG list extraction
#' @param merged_df Merged data frame from merge_all_results
#' @param time_point Specific time point ("6h" or "24h")
#' @param fc_cutoff Fold change cutoff (linear scale, default 1.5)
#' @param fdr_cutoff FDR cutoff
#' @return List of gene sets for the specified time point
identify_changed_genes_timepoint <- function(merged_df, 
                                             time_point = "6h",
                                             fc_cutoff = 1.5, #a linear number, the function takes the log2 of this, if set to 2, log2(2) = 1, so genes >1 or <-1. 
                                             fdr_cutoff = 0.05) {
  
  cell_lines <- c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1")
  changed_genes <- list()
  
  for (cell_line in cell_lines) {
    # Get column names for this time point
    logfc_col <- paste0("logFC_TgvDMSO_", time_point, "_", cell_line)
    fdr_col <- paste0("FDR_TgvDMSO_", time_point, "_", cell_line)
    
    if (logfc_col %in% colnames(merged_df) && fdr_col %in% colnames(merged_df)) {
      # Convert columns to numeric
      logfc_values <- as.numeric(merged_df[[logfc_col]])
      fdr_values <- as.numeric(merged_df[[fdr_col]])
      
      # Filter for significantly changed genes
      valid_rows <- !is.na(logfc_values) & !is.na(fdr_values) & 
        merged_df[[fdr_col]] != ""
      
      sig_genes <- merged_df$genes[
        valid_rows &
          abs(logfc_values) >= log2(fc_cutoff) &
          fdr_values < fdr_cutoff
      ]
      
      changed_genes[[cell_line]] <- sig_genes
      message(paste("  ", cell_line, "@", time_point, ":", 
                    length(sig_genes), "changed genes"))
    }
  }
  
  return(changed_genes)
}

################################################################################
### Section 4: VENN DIAGRAM FUNCTIONS
################################################################################

#' Create Venn diagram for gene comparisons
#' @param gene_sets List of gene sets
#' @param output_path Path to save the plot
#' @return ggplot object
create_venn_diagram <- function(gene_sets, output_path = "analysis/venn_diagrams") {
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Create Venn diagram
  venn_plot <- ggVennDiagram(
    gene_sets,
    label_color = "black",
    set_size = 5,
    label_alpha = 0,
    set_color = c("blue", "red", "green", "purple", "orange", "brown")
  ) +
    scale_fill_gradient(low = "white", high = "red") +
    labs(
      title = "Differentially Expressed Genes Across Cell Lines",
      subtitle = "Tg treatment (6h and 24h combined)"
    ) +
    theme_minimal()
  
  # Save plot
  ggsave(
    file.path(output_path, "Venn_diagram_all_cell_lines.svg"),
    plot = venn_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  return(venn_plot)
}

################################################################################
### Section 5: UPSET PLOT FUNCTIONS
################################################################################

#' Generate upset plots with consistent gene counting
#' @param merged_file Path to merged data file
#' @param fc_cutoff Fold change cutoff (linear scale)
#' @param fdr_cutoff FDR cutoff
#' @param output_path Path to save plots
generate_consistent_upset_plots <- function(
    merged_file = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt",
    fc_cutoff = 1.5,
    fdr_cutoff = 0.05,
    output_path = "analysis/upset_plots",
    text_scale_factor = 1.3
) {
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Load data
  message("\nLoading merged data...")
  data <- read.delim(merged_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Convert numeric columns
  numeric_cols <- grep("logFC|FDR", colnames(data))
  for (col in numeric_cols) {
    data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
  }
  
  cell_lines <- c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1")
  
  # Create gene sets for each time point separately
  gene_sets_6h <- list()
  gene_sets_24h <- list()
  gene_sets_both <- list()
  
  for (cell_line in cell_lines) {
    # 6h gene set
    logfc_col_6h <- paste0("logFC_TgvDMSO_6h_", cell_line)
    fdr_col_6h <- paste0("FDR_TgvDMSO_6h_", cell_line)
    
    if (all(c(logfc_col_6h, fdr_col_6h) %in% colnames(data))) {
      sig_genes_6h <- data$genes[
        !is.na(data[[logfc_col_6h]]) & 
          !is.na(data[[fdr_col_6h]]) &
          abs(data[[logfc_col_6h]]) > log2(fc_cutoff) &
          data[[fdr_col_6h]] < fdr_cutoff
      ]
      gene_sets_6h[[paste0(cell_line, "_6h")]] <- sig_genes_6h
      gene_sets_both[[paste0(cell_line, "_6h")]] <- sig_genes_6h
    }
    
    # 24h gene set
    logfc_col_24h <- paste0("logFC_TgvDMSO_24h_", cell_line)
    fdr_col_24h <- paste0("FDR_TgvDMSO_24h_", cell_line)
    
    if (all(c(logfc_col_24h, fdr_col_24h) %in% colnames(data))) {
      sig_genes_24h <- data$genes[
        !is.na(data[[logfc_col_24h]]) & 
          !is.na(data[[fdr_col_24h]]) &
          abs(data[[logfc_col_24h]]) > log2(fc_cutoff) &
          data[[fdr_col_24h]] < fdr_cutoff
      ]
      gene_sets_24h[[paste0(cell_line, "_24h")]] <- sig_genes_24h
      gene_sets_both[[paste0(cell_line, "_24h")]] <- sig_genes_24h
    }
  }
  
  # Print summary
  message("\n=== Gene Set Sizes ===")
  message("\n6h time point:")
  for (set_name in names(gene_sets_6h)) {
    message(sprintf("  %s: %d genes", set_name, length(gene_sets_6h[[set_name]])))
  }
  
  message("\n24h time point:")
  for (set_name in names(gene_sets_24h)) {
    message(sprintf("  %s: %d genes", set_name, length(gene_sets_24h[[set_name]])))
  }
  
  # Create upset plots using UpSetR
  if (length(gene_sets_6h) > 0) {
    # Convert to matrix format for UpSetR
    all_genes_6h <- unique(unlist(gene_sets_6h))
    upset_matrix_6h <- matrix(0, nrow = length(all_genes_6h), 
                              ncol = length(gene_sets_6h))
    rownames(upset_matrix_6h) <- all_genes_6h
    colnames(upset_matrix_6h) <- names(gene_sets_6h)
    
    for (i in seq_along(gene_sets_6h)) {
      upset_matrix_6h[all_genes_6h %in% gene_sets_6h[[i]], i] <- 1
    }
    
    upset_df_6h <- as.data.frame(upset_matrix_6h)
    
    # Create plot
    pdf(file.path(output_path, "upset_plot_6h.pdf"), width = 12, height = 8)
    print(upset(upset_df_6h, 
          sets = names(gene_sets_6h),
          order.by = "freq",
          decreasing = TRUE,
          number.angles = 30,
          point.size = 3,
          line.size = 1,
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size",
          text.scale = FIGURE_CONFIG$base_text_scale * text_scale_factor,
          mb.ratio = c(0.65, 0.35),
          nintersects = 43)
          )
    dev.off()
    
    message("✓ 6h upset plot saved")
  }
  
  if (length(gene_sets_24h) > 0) {
    # Convert to matrix format for UpSetR
    all_genes_24h <- unique(unlist(gene_sets_24h))
    upset_matrix_24h <- matrix(0, nrow = length(all_genes_24h), 
                               ncol = length(gene_sets_24h))
    rownames(upset_matrix_24h) <- all_genes_24h
    colnames(upset_matrix_24h) <- names(gene_sets_24h)
    
    for (i in seq_along(gene_sets_24h)) {
      upset_matrix_24h[all_genes_24h %in% gene_sets_24h[[i]], i] <- 1
    }
    
    upset_df_24h <- as.data.frame(upset_matrix_24h)
    
    # Create plot
    pdf(file.path(output_path, "upset_plot_24h.pdf"), width = 12, height = 8)
    print(upset(upset_df_24h, 
          sets = names(gene_sets_24h),
          order.by = "freq",
          decreasing = TRUE,
          number.angles = 30,
          point.size = 3,
          line.size = 1,
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size",
          text.scale = FIGURE_CONFIG$base_text_scale * text_scale_factor,
          mb.ratio = c(0.65, 0.35),
          nintersects = 43)
        )
    dev.off()
    
    message("✓ 24h upset plot saved")
  }
  
  if (length(gene_sets_both) > 0) {
    # Convert to matrix format for UpSetR
    all_genes_both <- unique(unlist(gene_sets_both))
    upset_matrix_both <- matrix(0, nrow = length(all_genes_both), 
                                ncol = length(gene_sets_both))
    rownames(upset_matrix_both) <- all_genes_both
    colnames(upset_matrix_both) <- names(gene_sets_both)
    
    for (i in seq_along(gene_sets_both)) {
      upset_matrix_both[all_genes_both %in% gene_sets_both[[i]], i] <- 1
    }
    
    upset_df_both <- as.data.frame(upset_matrix_both)
    
    # Create plot
    pdf(file.path(output_path, "upset_plot_both.pdf"), width = 16, height = 10)
    print(upset(upset_df_both, 
          sets = names(gene_sets_both),
          order.by = "freq",
          decreasing = TRUE,
          number.angles = 30,
          point.size = 2.5,
          line.size = 0.8,
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size",
          text.scale = FIGURE_CONFIG$base_text_scale_combined * text_scale_factor,
          mb.ratio = c(0.65, 0.35),
          nintersects = 60)
          )
    dev.off()
    
    message("✓ Combined upset plot saved")
  }
  
  # Save gene sets for reference
  saveRDS(list(
    gene_sets_6h = gene_sets_6h,
    gene_sets_24h = gene_sets_24h,
    gene_sets_both = gene_sets_both
  ), file.path(output_path, "upset_gene_sets.rds"))
  
  # Create summary statistics
  summary_stats <- data.frame(
    Set = character(),
    Size = integer(),
    Time = character(),
    stringsAsFactors = FALSE
  )
  
  for (set_name in names(gene_sets_both)) {
    parts <- strsplit(set_name, "_")[[1]]
    summary_stats <- rbind(summary_stats, data.frame(
      Set = set_name,
      Size = length(gene_sets_both[[set_name]]),
      Time = parts[2],
      CellLine = parts[1],
      stringsAsFactors = FALSE
    ))
  }
  
  write.csv(summary_stats, file.path(output_path, "upset_summary_stats.csv"), 
            row.names = FALSE)
  
  message("\n=== Upset Plot Generation Complete ===")
  
  return(list(
    gene_sets_6h = gene_sets_6h,
    gene_sets_24h = gene_sets_24h,
    gene_sets_both = gene_sets_both,
    summary = summary_stats
  ))
}

# Print detailed intersection analysis
analyze_intersections <- function(gene_sets) {
  all_genes <- unique(unlist(gene_sets))
  intersection_matrix <- sapply(gene_sets, function(set) all_genes %in% set)
  rownames(intersection_matrix) <- all_genes
  
  # Calculate intersection sizes
  intersection_ids <- apply(intersection_matrix, 1, function(row) paste(as.integer(row), collapse = ""))
  intersection_counts <- table(intersection_ids)
  
  # Create detailed report
  intersection_report <- data.frame()
  
  for (id in names(intersection_counts)) {
    binary <- as.integer(strsplit(id, "")[[1]])
    sets_in_intersection <- names(gene_sets)[binary == 1]
    
    intersection_report <- rbind(intersection_report, data.frame(
      ID = id,
      Sets = paste(sets_in_intersection, collapse = ", "),
      N_Sets = length(sets_in_intersection),
      Size = as.numeric(intersection_counts[id]),
      stringsAsFactors = FALSE
    ))
  }
  
  intersection_report <- intersection_report %>%
    arrange(desc(Size))
  
  return(intersection_report)

}


################################################################################
### SECTION 6: UNIQUE GENE ANALYSIS
################################################################################
  
#' Extract unique genes for each cell line PER TIME POINT
#' @param merged_df Merged data frame with all results
#' @param fc_cutoff Fold change cutoff (linear scale)
#' @param fdr_cutoff FDR cutoff
#' @param output_path Path to save unique gene files
#' @return List containing unique genes for both time points
extract_unique_genes_per_timepoint <- function(merged_df,
                                               fc_cutoff = 1.5,
                                               fdr_cutoff = 0.05,
                                               output_path = "edgeR_output") {
  
  message("\n=== Extracting Unique Genes Per Time Point ===")
  
  # Create lists to store results
  unique_genes_6h <- list()
  unique_genes_24h <- list()
  
  # Process 6h time point
  message("\nProcessing 6h time point...")
  gene_sets_6h <- identify_changed_genes_timepoint(merged_df, "6h", fc_cutoff, fdr_cutoff)
  
  for (cell_line in names(gene_sets_6h)) {
    # Get genes unique to this cell line at 6h
    other_lines <- setdiff(names(gene_sets_6h), cell_line)
    other_genes <- unique(unlist(gene_sets_6h[other_lines]))
    unique_genes <- setdiff(gene_sets_6h[[cell_line]], other_genes)
    
    if (length(unique_genes) > 0) {
      # Extract data for unique genes
      unique_data <- merged_df %>%
        dplyr::filter(genes %in% unique_genes) %>%
        dplyr::select(genes, contains(cell_line))
      
      # Add average TPM if columns exist
      tpm_cols <- grep(paste0(cell_line, "_(DMSO|Tg)_"), colnames(unique_data), value = TRUE)
      if (length(tpm_cols) > 0) {
        unique_data <- unique_data %>%
          dplyr::mutate(avg_expression = rowMeans(dplyr::select(., all_of(tpm_cols)), na.rm = TRUE))
      }
      
      unique_genes_6h[[cell_line]] <- unique_data
      
      # Save to file with time point in filename
      write.table(
        unique_data,
        file = file.path(output_path, paste0("Unique_genes_", cell_line, "_6h.txt")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      
      message(paste("  ", cell_line, ": ", length(unique_genes), " unique genes at 6h"))
    } else {
      message(paste("  ", cell_line, ": 0 unique genes at 6h"))
    }
  }
  
  # Process 24h time point
  message("\nProcessing 24h time point...")
  gene_sets_24h <- identify_changed_genes_timepoint(merged_df, "24h", fc_cutoff, fdr_cutoff)
  
  for (cell_line in names(gene_sets_24h)) {
    # Get genes unique to this cell line at 24h
    other_lines <- setdiff(names(gene_sets_24h), cell_line)
    other_genes <- unique(unlist(gene_sets_24h[other_lines]))
    unique_genes <- setdiff(gene_sets_24h[[cell_line]], other_genes)
    
    if (length(unique_genes) > 0) {
      # Extract data for unique genes
      unique_data <- merged_df %>%
        dplyr::filter(genes %in% unique_genes) %>%
        dplyr::select(genes, contains(cell_line))
      
      # Add average TPM if columns exist
      tpm_cols <- grep(paste0(cell_line, "_(DMSO|Tg)_"), colnames(unique_data), value = TRUE)
      if (length(tpm_cols) > 0) {
        unique_data <- unique_data %>%
          dplyr::mutate(avg_expression = rowMeans(dplyr::select(., all_of(tpm_cols)), na.rm = TRUE))
      }
      
      unique_genes_24h[[cell_line]] <- unique_data
      
      # Save to file with time point in filename
      write.table(
        unique_data,
        file = file.path(output_path, paste0("Unique_genes_", cell_line, "_24h.txt")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      
      message(paste("  ", cell_line, ": ", length(unique_genes), " unique genes at 24h"))
    } else {
      message(paste("  ", cell_line, ": 0 unique genes at 24h"))
    }
  }
  
  return(list(
    unique_6h = unique_genes_6h,
    unique_24h = unique_genes_24h,
    gene_sets_6h = gene_sets_6h,
    gene_sets_24h = gene_sets_24h
  ))
}

################################################################################
### SECTION 7: GO ANALYSIS FOR UNIQUE GENES
################################################################################

#' Perform GO analysis on unique genes
#' @param unique_genes_list List from extract_unique_genes
#' @param output_path Path to save GO results
#' @return List of GO results
analyze_unique_genes_go <- function(unique_genes_list, 
                                   output_path = "analysis/GO_unique") {
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  go_results <- list()
  
  for (cell_line in names(unique_genes_list)) {
    # Check if unique_genes_list[[cell_line]] exists and has genes column
    if (is.null(unique_genes_list[[cell_line]]) || 
        !"genes" %in% colnames(unique_genes_list[[cell_line]])) {
      message(paste("  Skipping", cell_line, "- no unique genes data"))
      next
    }
    
    genes <- unique_genes_list[[cell_line]]$genes
    
    # Filter out empty strings and NAs
    genes <- genes[!is.na(genes) & genes != "" & genes != "NA"]
    
    message(paste("  Analyzing", length(genes), "unique genes for", cell_line))
    
    if (length(genes) > 0) {
      # Determine species for annotation
      org_db <- org.Mm.eg.db  # Default to mouse
      
      # Run GO analysis for each ontology
      for (ont in c("BP", "MF", "CC")) {
        tryCatch({
          go_result <- enrichGO(
            gene = genes,
            OrgDb = org_db,
            keyType = "SYMBOL",
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05
          )
          
          if (!is.null(go_result) && nrow(go_result@result) > 0) {
            message(paste("    Found", nrow(go_result@result), "enriched", ont, "terms"))
            
            # Save results
            write.table(
              go_result@result,
              file = file.path(output_path, 
                             paste0("GO_", ont, "_", cell_line, "_unique.txt")),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE
            )
            
            # Create and save plot
            tryCatch({
              p <- barplot(go_result, showCategory = min(15, nrow(go_result@result)), 
                          title = paste("GO", ont, "-", cell_line, "Unique Genes"))
              
              ggsave(
                file.path(output_path, 
                         paste0("GO_", ont, "_", cell_line, "_unique.svg")),
                plot = p,
                width = 10,
                height = 8,
                dpi = 300
              )
            }, error = function(e) {
              message(paste("      Could not create plot:", e$message))
            })
            
            go_results[[paste(cell_line, ont, sep = "_")]] <- go_result
          } else {
            message(paste("    No significant", ont, "terms found"))
          }
        }, error = function(e) {
          message(paste("    Error in", ont, "analysis:", e$message))
        })
      }
    } else {
      message(paste("  No valid genes for GO analysis in", cell_line))
    }
  }
  
  if (length(go_results) == 0) {
    message("  No GO enrichment results obtained for any cell line")
  }
  
  return(go_results)
}

################################################################################
###  Section 8: FUNCTIONS for Pairwise comparisons of log2FC for all cell lines
################################################################################

#' Create scatter plot comparing two cell lines with gene labels
#' @param data Merged data with all results
#' @param cell_x First cell line (x-axis)
#' @param cell_y Second cell line (y-axis) 
#' @param time_point Time point ("6h" or "24h")
#' @param fdr_cutoff FDR cutoff for significance
#' @param n_labels Number of top genes to label (default: 15)
#' @return ggplot object
create_comparison_plot <- function(data, cell_x, cell_y, time_point, fdr_cutoff = 0.05, n_labels = 15) {
  
  # Get column names
  logfc_x <- paste0("logFC_TgvDMSO_", time_point, "_", cell_x)
  logfc_y <- paste0("logFC_TgvDMSO_", time_point, "_", cell_y)
  fdr_x <- paste0("FDR_TgvDMSO_", time_point, "_", cell_x)
  fdr_y <- paste0("FDR_TgvDMSO_", time_point, "_", cell_y)
  
  # Check if columns exist
  if (!all(c(logfc_x, logfc_y, fdr_x, fdr_y) %in% colnames(data))) {
    warning(paste("Missing columns for comparison:", cell_x, "vs", cell_y, "at", time_point))
    return(NULL)
  }
  
  # Prepare comparison data
  comp_data <- data %>%
    dplyr::select(genes, all_of(c(logfc_x, logfc_y, fdr_x, fdr_y))) %>%
    dplyr::rename(x = !!logfc_x, y = !!logfc_y, fdr_x = !!fdr_x, fdr_y = !!fdr_y) %>%
    dplyr::filter(!is.na(x), !is.na(y), !is.na(fdr_x), !is.na(fdr_y))
  
  # Add significance categories
  comp_data <- comp_data %>%
    mutate(
      significant = case_when(
        fdr_x < fdr_cutoff & fdr_y < fdr_cutoff ~ "Both Significant",
        fdr_x < fdr_cutoff ~ paste(cell_x, "Only"),
        fdr_y < fdr_cutoff ~ paste(cell_y, "Only"),
        TRUE ~ "Not Significant"
      )
    )
  
  # Select genes for labeling using smart selection strategy
  label_genes <- comp_data %>%
    filter(significant != "Not Significant") %>%  # Only significant genes
    mutate(
      # Enhanced distance metric: emphasize fold change magnitude + significance
      distance = abs(x) + abs(y) +  # Individual fold change magnitudes
        sqrt(x^2 + y^2) +   # Combined fold change magnitude (Euclidean)
        (-log10(pmax(fdr_x, 1e-300, na.rm = TRUE))) + 
        (-log10(pmax(fdr_y, 1e-300, na.rm = TRUE)))
    ) %>%
    top_n(n_labels, distance)  # Select top N genes for labeling
  
  # Add gene_label column
  comp_data <- comp_data %>%
    mutate(gene_label = ifelse(genes %in% label_genes$genes, genes, ""))
  
  # Calculate correlation
  cor_pearson <- cor(comp_data$x, comp_data$y, use = "complete.obs", method = "pearson")
  cor_spearman <- cor(comp_data$x, comp_data$y, use = "complete.obs", method = "spearman")
  
  # Define colors
  colors <- c(
    "Both Significant" = "red",
    "Not Significant" = "grey70"
  )
  colors[paste(cell_x, "Only")] <- "blue"
  colors[paste(cell_y, "Only")] <- "green"
  
  # Create plot
  p <- ggplot(comp_data, aes(x = x, y = y)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1) +
    # Highlight labeled genes with black border
    geom_point(data = comp_data[comp_data$gene_label != "", ], 
               aes(fill = significant), 
               shape = 21, color = "black", size = 2) +
    scale_color_manual(values = colors, name = "Significance") +
    scale_fill_manual(values = colors, guide = "none") +  # No separate legend for fill
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
      x = paste(cell_x, "Log2 FC"),
      y = paste(cell_y, "Log2 FC"),
      title = paste(cell_x, "vs", cell_y, "-", time_point),
      subtitle = sprintf("Pearson r = %.3f, Spearman ρ = %.3f", cor_pearson, cor_spearman)
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.position = "none"  # Remove legend for grid layout
    ) +
    coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8))
  
  # Add gene labels for selected genes only
  labeled_data <- comp_data[comp_data$gene_label != "", ]
  if (nrow(labeled_data) > 0) {
    p <- p + geom_text_repel(
      data = labeled_data,
      aes(label = gene_label),
      size = 2.5,
      box.padding = 0.3,
      max.overlaps = 20,
      color = "black",
      segment.color = "grey50",
      segment.size = 0.3
    )
  }
  
  return(p)
}

#' Create comparison plot with larger labels for individual viewing
#' @param data Merged data with all results
#' @param cell_x First cell line (x-axis)
#' @param cell_y Second cell line (y-axis) 
#' @param time_point Time point ("6h" or "24h")
#' @param fdr_cutoff FDR cutoff for significance
#' @param n_labels Number of top genes to label (default: 20)
#' @return ggplot object
create_comparison_plot_large <- function(data, cell_x, cell_y, time_point, fdr_cutoff = 0.05, n_labels = 20) {
  
  # Get column names
  logfc_x <- paste0("logFC_TgvDMSO_", time_point, "_", cell_x)
  logfc_y <- paste0("logFC_TgvDMSO_", time_point, "_", cell_y)
  fdr_x <- paste0("FDR_TgvDMSO_", time_point, "_", cell_x)
  fdr_y <- paste0("FDR_TgvDMSO_", time_point, "_", cell_y)
  
  # Check if columns exist
  if (!all(c(logfc_x, logfc_y, fdr_x, fdr_y) %in% colnames(data))) {
    warning(paste("Missing columns for comparison:", cell_x, "vs", cell_y, "at", time_point))
    return(NULL)
  }
  
  # Prepare comparison data
  comp_data <- data %>%
    dplyr::select(genes, all_of(c(logfc_x, logfc_y, fdr_x, fdr_y))) %>%
    dplyr::rename(x = !!logfc_x, y = !!logfc_y, fdr_x = !!fdr_x, fdr_y = !!fdr_y) %>%
    dplyr::filter(!is.na(x), !is.na(y), !is.na(fdr_x), !is.na(fdr_y))
  
  # Add significance categories
  comp_data <- comp_data %>%
    mutate(
      significant = case_when(
        fdr_x < fdr_cutoff & fdr_y < fdr_cutoff ~ "Both Significant",
        fdr_x < fdr_cutoff ~ paste(cell_x, "Only"),
        fdr_y < fdr_cutoff ~ paste(cell_y, "Only"),
        TRUE ~ "Not Significant"
      )
    )
  
  # Select genes for labeling
  label_genes <- comp_data %>%
    filter(significant != "Not Significant") %>%
    mutate(
      distance = abs(x) + abs(y) + sqrt(x^2 + y^2) +
        (-log10(pmax(fdr_x, 1e-300, na.rm = TRUE))) + 
        (-log10(pmax(fdr_y, 1e-300, na.rm = TRUE)))
    ) %>%
    top_n(n_labels, distance)
  
  comp_data <- comp_data %>%
    mutate(gene_label = ifelse(genes %in% label_genes$genes, genes, ""))
  
  # Calculate correlation
  cor_pearson <- cor(comp_data$x, comp_data$y, use = "complete.obs", method = "pearson")
  cor_spearman <- cor(comp_data$x, comp_data$y, use = "complete.obs", method = "spearman")
  
  # Define colors
  colors <- c(
    "Both Significant" = "red",
    "Not Significant" = "grey70"
  )
  colors[paste(cell_x, "Only")] <- "blue"
  colors[paste(cell_y, "Only")] <- "green"
  
  # Create larger plot for individual viewing
  p <- ggplot(comp_data, aes(x = x, y = y)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    geom_point(data = comp_data[comp_data$gene_label != "", ], 
               aes(fill = significant), 
               shape = 21, color = "black", size = 3) +
    scale_color_manual(values = colors, name = "Significance") +
    scale_fill_manual(values = colors, guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    labs(
      x = paste(cell_x, "Log2 Fold Change"),
      y = paste(cell_y, "Log2 Fold Change"),
      title = paste("Cell Line Comparison:", cell_x, "vs", cell_y, "-", time_point),
      subtitle = sprintf("Pearson r = %.3f, Spearman ρ = %.3f", cor_pearson, cor_spearman)
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ) +
    coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8))
  
  # Add gene labels
  labeled_data <- comp_data[comp_data$gene_label != "", ]
  if (nrow(labeled_data) > 0) {
    p <- p + geom_text_repel(
      data = labeled_data,
      aes(label = gene_label),
      size = 3.5,
      box.padding = 0.5,
      max.overlaps = 25,
      color = "black",
      segment.color = "grey50",
      segment.size = 0.4
    )
  }
  
  return(p)
}

#' Generate all pairwise comparison plots
#' @param data Merged data with all results
#' @param cell_lines Vector of cell line names
#' @param output_path Path to save plots
#' @return List of plots
generate_comparison_plots <- function(data, 
                                      cell_lines = c("MIN6", "aTC1", "MGN3", "GLUTag", "PCCL3", "QGP1"),
                                      output_path = "analysis/comparison_plots") {
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Generate all pairwise combinations (excluding self-comparisons)
  pairs <- combn(cell_lines, 2, simplify = FALSE)
  
  plots_6h <- list()
  plots_24h <- list()
  
  message("\n=== Generating Cell Line Comparison Plots ===")
  
  for (pair in pairs) {
    cell_x <- pair[1]
    cell_y <- pair[2]
    
    # 6h comparison
    p_6h <- create_comparison_plot(data, cell_x, cell_y, "6h")
    if (!is.null(p_6h)) {
      plots_6h[[paste(cell_x, cell_y, sep = "_vs_")]] <- p_6h
      
      # Save individual plot
      ggsave(
        filename = file.path(output_path, paste0("Comparison_", cell_x, "_vs_", cell_y, "_6h.pdf")),
        plot = p_6h,
        width = 5,
        height = 5,
        dpi = 300
      )
    }
    
    # 24h comparison
    p_24h <- create_comparison_plot(data, cell_x, cell_y, "24h")
    if (!is.null(p_24h)) {
      plots_24h[[paste(cell_x, cell_y, sep = "_vs_")]] <- p_24h
      
      # Save individual plot
      ggsave(
        filename = file.path(output_path, paste0("Comparison_", cell_x, "_vs_", cell_y, "_24h.pdf")),
        plot = p_24h,
        width = 5,
        height = 5,
        dpi = 300
      )
    }
    
    message(paste("  ✓ Generated comparison:", cell_x, "vs", cell_y))
  }
  
  return(list(plots_6h = plots_6h, plots_24h = plots_24h))
}

#' Create combined comparison plot grid
#' @param comparison_plots List of comparison plots from generate_comparison_plots
#' @param output_path Path to save the combined plot
#' @return Combined plot
create_comparison_grid <- function(comparison_plots, 
                                   output_path = "analysis/comparison_plots") {
  
  library(patchwork)
  
  # Extract plots
  plots_6h <- comparison_plots$plots_6h
  plots_24h <- comparison_plots$plots_24h
  
  # Define the order for a nice layout (5x3 grid for each time point)
  # Row 1: MIN6 vs others
  # Row 2: aTC1 vs others (excluding MIN6)
  # Row 3: MGN3 vs others (excluding MIN6, aTC1)
  # Row 4: GLUTag vs others (excluding MIN6, aTC1, MGN3)
  # Row 5: PCCL3 vs QGP1
  
  plot_order_6h <- c(
    "MIN6_vs_aTC1", "MIN6_vs_MGN3", "MIN6_vs_GLUTag",
    "MIN6_vs_PCCL3", "MIN6_vs_QGP1", "aTC1_vs_MGN3",
    "aTC1_vs_GLUTag", "aTC1_vs_PCCL3", "aTC1_vs_QGP1",
    "MGN3_vs_GLUTag", "MGN3_vs_PCCL3", "MGN3_vs_QGP1",
    "GLUTag_vs_PCCL3", "GLUTag_vs_QGP1", "PCCL3_vs_QGP1"
  )
  
  # Create 6h grid
  if (length(plots_6h) > 0) {
    # Get plots in order
    ordered_plots_6h <- list()
    for (pair_name in plot_order_6h) {
      if (pair_name %in% names(plots_6h)) {
        ordered_plots_6h[[pair_name]] <- plots_6h[[pair_name]]
      }
    }
    
    # Combine plots using patchwork
    combined_6h <- wrap_plots(ordered_plots_6h, ncol = 3) +
      plot_annotation(
        title = "Cell Line Comparisons - 6h Thapsigargin Treatment",
        subtitle = "Log2 Fold Change Correlations",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
      )
    
    # Save 6h grid
    ggsave(
      filename = file.path(output_path, "Comparison_Grid_6h.pdf"),
      plot = combined_6h,
      width = 15,
      height = 20,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(output_path, "Comparison_Grid_6h.png"),
      plot = combined_6h,
      width = 15,
      height = 20,
      dpi = 300
    )
    
    message("✓ 6h comparison grid saved")
  }
  
  # Create 24h grid
  if (length(plots_24h) > 0) {
    ordered_plots_24h <- list()
    for (pair_name in plot_order_6h) {  # Same order as 6h
      if (pair_name %in% names(plots_24h)) {
        ordered_plots_24h[[pair_name]] <- plots_24h[[pair_name]]
      }
    }
    
    combined_24h <- wrap_plots(ordered_plots_24h, ncol = 3) +
      plot_annotation(
        title = "Cell Line Comparisons - 24h Thapsigargin Treatment",
        subtitle = "Log2 Fold Change Correlations",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
      )
    
    # Save 24h grid
    ggsave(
      filename = file.path(output_path, "Comparison_Grid_24h.pdf"),
      plot = combined_24h,
      width = 15,
      height = 20,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(output_path, "Comparison_Grid_24h.png"),
      plot = combined_24h,
      width = 15,
      height = 20,
      dpi = 300
    )
    
    message("✓ 24h comparison grid saved")
  }
  
  # Create combined grid with both time points
  if (length(plots_6h) > 0 && length(plots_24h) > 0) {
    # Select key comparisons for paper figure
    key_pairs <- c("MIN6_vs_aTC1", "MIN6_vs_QGP1", "aTC1_vs_QGP1",
                   "MIN6_vs_GLUTag", "aTC1_vs_GLUTag", "GLUTag_vs_QGP1")
    
    key_plots <- list()
    for (pair in key_pairs) {
      if (pair %in% names(plots_6h)) {
        key_plots[[paste0(pair, "_6h")]] <- plots_6h[[pair]]
      }
      if (pair %in% names(plots_24h)) {
        key_plots[[paste0(pair, "_24h")]] <- plots_24h[[pair]]
      }
    }
    
    combined_key <- wrap_plots(key_plots, ncol = 3) +
      plot_annotation(
        title = "Key Cell Line Comparisons",
        subtitle = "6h and 24h Thapsigargin Treatment",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
      )
    
    ggsave(
      filename = file.path(output_path, "Comparison_Grid_Key.pdf"),
      plot = combined_key,
      width = 15,
      height = 10,
      dpi = 300
    )
    
    message("✓ Key comparison grid saved")
  }
  
  return(list(combined_6h = combined_6h, combined_24h = combined_24h))
}

#' Run comparison analysis
#' @param merged_file Path to merged data file
run_comparison_analysis <- function(merged_file = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt") {
  
  # Load data
  message("\nLoading merged data...")
  data <- read.delim(merged_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Convert numeric columns
  numeric_cols <- grep("logFC|FDR|logCPM|TPM", colnames(data))
  for (col in numeric_cols) {
    data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
  }
  
  # Generate comparison plots
  comparison_plots <- generate_comparison_plots(data)
  
  # Create combined grids
  create_comparison_grid(comparison_plots)
  
  message("\n=== Comparison Analysis Complete ===")
  
  return(comparison_plots)
}

#' Extract concordant and discordant genes between cell line pairs
#' @param data Merged data frame with all results
#' @param cell_line_pairs List of cell line pairs to analyze (default: only MIN6 vs aTC1)
#' @param fdr_cutoff FDR cutoff for significance (default: 0.05)
#' @param logfc_diff_cutoff Minimum absolute difference in log2FC between cell lines for DISCORDANT genes only (default: 1.5)
#' @param concordant_logfc_cutoff Minimum |log2FC| for CONCORDANT genes in both cell lines (default: 1)
#' @param output_path Path to save Excel file
#' @return List of data frames with concordant/discordant genes
extract_concordant_discordant_genes <- function(data,
                                                cell_line_pairs = list(c("MIN6", "aTC1")),
                                                fdr_cutoff = 0.05,
                                                logfc_diff_cutoff = 1.5,
                                                concordant_logfc_cutoff = 1,
                                                output_path = "analysis/comparison_plots") {
  
  library(openxlsx)
  
  message("\n=== Extracting Concordant/Discordant Genes ===")
  message(paste("  FDR cutoff:", fdr_cutoff))
  message(paste("  Log2FC difference cutoff (discordant):", logfc_diff_cutoff))
  message(paste("  Minimum |Log2FC| (concordant):", concordant_logfc_cutoff))
  
  # Create output directory if it doesn't exist
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize list to store all results
  all_results <- list()
  
  # Process each cell line pair
  for (pair in cell_line_pairs) {
    cell1 <- pair[1]
    cell2 <- pair[2]
    
    message(paste("\nProcessing:", cell1, "vs", cell2))
    
    # Process each time point
    for (time_point in c("6h", "24h")) {
      
      # Get column names
      logfc_col1 <- paste0("logFC_TgvDMSO_", time_point, "_", cell1)
      logfc_col2 <- paste0("logFC_TgvDMSO_", time_point, "_", cell2)
      fdr_col1 <- paste0("FDR_TgvDMSO_", time_point, "_", cell1)
      fdr_col2 <- paste0("FDR_TgvDMSO_", time_point, "_", cell2)
      
      # Check if columns exist
      required_cols <- c(logfc_col1, logfc_col2, fdr_col1, fdr_col2)
      if (!all(required_cols %in% colnames(data))) {
        warning(paste("Missing columns for", cell1, "vs", cell2, "at", time_point))
        next
      }
      
      # Filter for significant genes in both cell lines (base filtering)
      sig_data <- data %>%
        dplyr::select(genes, all_of(required_cols)) %>%
        dplyr::rename(
          logFC1 = !!logfc_col1,
          logFC2 = !!logfc_col2,
          FDR1 = !!fdr_col1,
          FDR2 = !!fdr_col2
        ) %>%
        dplyr::filter(
          !is.na(logFC1), !is.na(logFC2),
          !is.na(FDR1), !is.na(FDR2),
          FDR1 < fdr_cutoff,
          FDR2 < fdr_cutoff
        )
      
      # Categorize genes into quadrants
      quadrants <- list()
      
      # Quadrant 1: cell1_down_cell2_up (DISCORDANT - apply difference filter)
      quadrant_name <- paste0(cell1, "_down_", cell2, "_up_", time_point)
      quadrants[[quadrant_name]] <- sig_data %>%
        dplyr::filter(
          logFC1 < 0, 
          logFC2 > 0,
          abs(logFC1 - logFC2) > logfc_diff_cutoff  # Apply difference filter for discordant
        ) %>%
        dplyr::mutate(
          category = quadrant_name,
          logFC_diff = logFC2 - logFC1
        ) %>%
        dplyr::arrange(desc(abs(logFC_diff))) %>%
        dplyr::select(
          genes, 
          !!paste0(cell1, "_logFC") := logFC1,
          !!paste0(cell1, "_FDR") := FDR1,
          !!paste0(cell2, "_logFC") := logFC2,
          !!paste0(cell2, "_FDR") := FDR2,
          logFC_diff
        )
      
      # Quadrant 2: cell1_up_cell2_up (CONCORDANT up - apply minimum logFC filter)
      quadrant_name <- paste0(cell1, "_up_", cell2, "_up_", time_point)
      quadrants[[quadrant_name]] <- sig_data %>%
        dplyr::filter(
          logFC1 > 0, 
          logFC2 > 0,
          abs(logFC1) > concordant_logfc_cutoff,  # Both must exceed minimum threshold
          abs(logFC2) > concordant_logfc_cutoff
        ) %>%
        dplyr::mutate(
          category = quadrant_name,
          logFC_combined = abs(logFC1 + logFC2)
        ) %>%
        dplyr::arrange(desc(abs(logFC_combined))) %>%  # Sort by combined magnitude BEFORE renaming
        dplyr::select(
          genes,
          !!paste0(cell1, "_logFC") := logFC1,
          !!paste0(cell1, "_FDR") := FDR1,
          !!paste0(cell2, "_logFC") := logFC2,
          !!paste0(cell2, "_FDR") := FDR2,
        )
      
      # Quadrant 3: cell1_up_cell2_down (DISCORDANT - apply difference filter)
      quadrant_name <- paste0(cell1, "_up_", cell2, "_down_", time_point)
      quadrants[[quadrant_name]] <- sig_data %>%
        dplyr::filter(
          logFC1 > 0, 
          logFC2 < 0,
          abs(logFC1 - logFC2) > logfc_diff_cutoff  # Apply difference filter for discordant
        ) %>%
        dplyr::mutate(
          category = quadrant_name,
          logFC_diff = logFC1 - logFC2
        ) %>%
        dplyr::arrange(desc(abs(logFC_diff))) %>%
        dplyr::select(
          genes,
          !!paste0(cell1, "_logFC") := logFC1,
          !!paste0(cell1, "_FDR") := FDR1,
          !!paste0(cell2, "_logFC") := logFC2,
          !!paste0(cell2, "_FDR") := FDR2,
          logFC_diff
        )
      
      # Quadrant 4: cell1_down_cell2_down (CONCORDANT down - apply minimum logFC filter)
      quadrant_name <- paste0(cell1, "_down_", cell2, "_down_", time_point)
      quadrants[[quadrant_name]] <- sig_data %>%
        dplyr::filter(
          logFC1 < 0, 
          logFC2 < 0,
          abs(logFC1) > concordant_logfc_cutoff,  # Both must exceed minimum threshold
          abs(logFC2) > concordant_logfc_cutoff
        ) %>%
        dplyr::mutate(
          category = quadrant_name,
          logFC_combined = abs(logFC1 + logFC2)
        ) %>%
        dplyr::arrange(desc(abs(logFC_combined))) %>%  # Sort by combined magnitude
        dplyr::select(
          genes,
          !!paste0(cell1, "_logFC") := logFC1,
          !!paste0(cell1, "_FDR") := FDR1,
          !!paste0(cell2, "_logFC") := logFC2,
          !!paste0(cell2, "_FDR") := FDR2,
        ) 
      
      # Add to results
      for (q_name in names(quadrants)) {
        all_results[[q_name]] <- quadrants[[q_name]]
        message(paste("    ", q_name, ":", nrow(quadrants[[q_name]]), "genes"))
      }
    }
  }
  
  # Export to Excel
  if (length(all_results) > 0) {
    # Create workbook
    wb <- createWorkbook()
    
    # Add each quadrant as a sheet
    for (sheet_name in names(all_results)) {
      # Excel has a 31 character limit for sheet names
      short_name <- sheet_name
      if (nchar(sheet_name) > 31) {
        # Shorten the name while keeping it informative
        short_name <- gsub("MIN6", "M6", sheet_name)
        short_name <- gsub("aTC1", "aT", short_name)
        if (nchar(short_name) > 31) {
          short_name <- substr(short_name, 1, 31)
        }
      }
      
      addWorksheet(wb, short_name)
      writeData(wb, short_name, all_results[[sheet_name]])
      
      # Add formatting
      headerStyle <- createStyle(
        fontSize = 12,
        fontColour = "#FFFFFF",
        halign = "center",
        fgFill = "#4472C4",
        border = "TopBottomLeftRight",
        textDecoration = "bold"
      )
      addStyle(wb, short_name, headerStyle, rows = 1, 
               cols = 1:ncol(all_results[[sheet_name]]), gridExpand = TRUE)
      
      # Auto-width columns
      setColWidths(wb, short_name, cols = 1:ncol(all_results[[sheet_name]]), widths = "auto")
    }
    
    # Add summary sheet
    summary_data <- data.frame(
      Category = names(all_results),
      N_Genes = sapply(all_results, nrow),
      stringsAsFactors = FALSE
    )
    
    addWorksheet(wb, "Summary")
    writeData(wb, "Summary", summary_data)
    addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
    setColWidths(wb, "Summary", cols = 1:2, widths = "auto")
    
    # Save workbook
    excel_file <- file.path(output_path, paste0(cell_line_pairs[[1]][1], "_vs_", 
                                                cell_line_pairs[[1]][2], 
                                                "_concordant_discordant_genes.xlsx"))
    saveWorkbook(wb, excel_file, overwrite = TRUE)
    message(paste("\n✓ Results saved to:", excel_file))
  }
  
  return(all_results)
}

#' Run concordant/discordant analysis with visualization
#' @param merged_file Path to merged data file
#' @param cell_line_pairs List of cell line pairs to analyze
#' @param fdr_cutoff FDR cutoff
#' @param logfc_diff_cutoff Minimum log2FC difference between cell lines for discordant genes
#' @param concordant_logfc_cutoff Minimum |log2FC| for concordant genes
#' @param create_plots Whether to create visualization plots
run_concordant_discordant_analysis <- function(
    merged_file = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt",
    cell_line_pairs = list(c("MIN6", "aTC1")),  # Add more pairs as needed: list(c("MIN6", "aTC1"), c("MIN6", "QGP1"))
    fdr_cutoff = 0.05,
    logfc_diff_cutoff = 1.5,
    concordant_logfc_cutoff = 1,
    create_plots = TRUE) {
  
  # Load data
  message("\nLoading merged data...")
  data <- read.delim(merged_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Convert numeric columns
  numeric_cols <- grep("logFC|FDR", colnames(data))
  for (col in numeric_cols) {
    data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
  }
  
  # Extract concordant/discordant genes
  results <- extract_concordant_discordant_genes(
    data = data,
    cell_line_pairs = cell_line_pairs,
    fdr_cutoff = fdr_cutoff,
    logfc_diff_cutoff = logfc_diff_cutoff,
    concordant_logfc_cutoff = concordant_logfc_cutoff
  )
  
  # Optional: Create summary visualization
  if (create_plots && length(results) > 0) {
    library(ggplot2)
    library(patchwork)
    
    # Create bar plot of gene counts
    summary_df <- data.frame(
      Category = names(results),
      Count = sapply(results, nrow),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        Time = ifelse(grepl("_6h$", Category), "6h", "24h"),
        Type = case_when(
          grepl("_up_.*_up_", Category) ~ "Concordant Up",
          grepl("_down_.*_down_", Category) ~ "Concordant Down",
          grepl("_up_.*_down_", Category) ~ "Discordant (Up/Down)",
          grepl("_down_.*_up_", Category) ~ "Discordant (Down/Up)",
          TRUE ~ "Other"
        )
      )
    
    p_summary <- ggplot(summary_df, aes(x = Type, y = Count, fill = Time)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("6h" = "#E69F00", "24h" = "#56B4E9")) +
      labs(
        title = paste("Concordant/Discordant Genes:", cell_line_pairs[[1]][1], "vs", cell_line_pairs[[1]][2]),
        subtitle = paste("FDR <", fdr_cutoff, 
                         "| Discordant: |ΔLog2FC| >", logfc_diff_cutoff,
                         "| Concordant: |Log2FC| >", concordant_logfc_cutoff),
        x = "Category",
        y = "Number of Genes"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    # Save plot
    output_path <- "analysis/comparison_plots"
    ggsave(
      filename = file.path(output_path, paste0(cell_line_pairs[[1]][1], "_vs_", 
                                               cell_line_pairs[[1]][2], 
                                               "_concordant_discordant_summary.pdf")),
      plot = p_summary,
      width = 8,
      height = 6,
      dpi = 300
    )
    
    message("✓ Summary plot saved")
  }
  
  return(results)
}

################################################################################
### SECTION 9: MAIN COMPARATIVE ANALYSIS WORKFLOW
################################################################################
#' Run complete comparative analysis
#' @param cell_lines Vector of cell line names to analyze
#' @param debug Whether to print debugging information
#' @param mouse_symbol_file Optional file for updating mouse gene symbols
#' @return List of all results

run_comparative_analysis <- function(cell_lines = c("MIN6", "aTC1", "MGN3", 
                                                    "GLUTag", "PCCL3", "QGP1"),
                                     debug = FALSE,
                                     mouse_symbol_file = NULL,
                                     run_GO_analysis = FALSE) {
  
  message("\n========================================")
  message("Starting Comparative Analysis Pipeline")
  message("========================================\n")
  
  # Load all edgeR results
  message("Step 1: Loading edgeR results...")
  all_results <- load_all_edger_results(cell_lines, 
                                        mouse_symbol_file = mouse_symbol_file)
  
  if (debug) {
    message("\nDebug: Loaded data dimensions:")
    for (cl in names(all_results)) {
      message(paste("  ", cl, ":", nrow(all_results[[cl]]), "rows x", 
                    ncol(all_results[[cl]]), "columns"))
    }
  }
  
  # Merge all results
  message("\nStep 2: Merging results...")
  merged_df <- merge_all_results(all_results)
  
  # Create output directory if it doesn't exist
  dir.create("merged edgeR and TPMs", showWarnings = FALSE, recursive = TRUE)
  
  # Save merged results
  write.table(
    merged_df,
    file = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  message(paste("  Merged data:", nrow(merged_df), "genes"))
  
  if (debug) {
    message("\nDebug: Column names in merged data:")
    message(paste("  ", paste(head(colnames(merged_df), 20), collapse = ", ")))
    if (ncol(merged_df) > 20) {
      message(paste("  ... and", ncol(merged_df) - 20, "more columns"))
    }
  }
  
  # Identify changed genes (for upset plots - still uses combined)
  message("\nStep 3: Identifying changed genes...")
  changed_genes <- identify_changed_genes(merged_df, debug = debug)
  
  # Initialize plot variables
  upset_plot <- NULL
  venn_plot <- NULL
  
  # upset plots - updated
  upset_results <- generate_consistent_upset_plots(
    merged_file = "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt",
    fc_cutoff = 1.5,
    fdr_cutoff = 0.05,
    output_path = "analysis/upset_plots"
  )
  
  # Analyze each time point
  message("\n=== Intersection Analysis ===")
  
  if (length(upset_results$gene_sets_6h) > 0) {
    report_6h <- analyze_intersections(upset_results$gene_sets_6h)
    write.csv(report_6h, "analysis/upset_plots/intersection_report_6h.csv", row.names = FALSE)
    message("\n6h - Top 10 intersections:")
    print(head(report_6h, 10))
  }
  
  if (length(upset_results$gene_sets_24h) > 0) {
    report_24h <- analyze_intersections(upset_results$gene_sets_24h)
    write.csv(report_24h, "analysis/upset_plots/intersection_report_24h.csv", row.names = FALSE)
    message("\n24h - Top 10 intersections:")
    print(head(report_24h, 10))
  }
  
  if (length(upset_results$gene_sets_both) > 0) {
    report_both <- analyze_intersections(upset_results$gene_sets_both)
    write.csv(report_both, "analysis/upset_plots/intersection_report_both.csv", row.names = FALSE)
    message("\nCombined - Top 10 intersections:")
    print(head(report_both, 10))
  }
  
  # Create Venn diagram
  message("\nStep 5: Creating Venn diagram...")
  tryCatch({
    venn_plot <- create_venn_diagram(changed_genes)
  }, error = function(e) {
    message(paste("  Warning: Could not create Venn diagram:", e$message))
    venn_plot <- NULL
  })
  
  # Extract unique genes PER TIME POINT
  message("\nStep 6: Extracting unique genes per time point...")
  unique_genes_results <- extract_unique_genes_per_timepoint(
    merged_df,
    fc_cutoff = 1.5,  # Use 1.5 fold change
    fdr_cutoff = 0.05,
    output_path = "edgeR_output"
  )
  
if(run_GO_analysis == TRUE){
  # GO analysis on unique genes - need to handle the new structure
  message("\nStep 7: Performing GO analysis on unique genes...")
  go_results_all <- list()
  
  # Analyze 6h unique genes
  if (length(unique_genes_results$unique_6h) > 0) {
    message("\n  Analyzing 6h unique genes...")
    tryCatch({
      go_results_6h <- analyze_unique_genes_go(
        unique_genes_results$unique_6h,
        output_path = "analysis/GO_unique/6h"
      )
      go_results_all[["6h"]] <- go_results_6h
    }, error = function(e) {
      message(paste("    Warning: GO analysis failed for 6h:", e$message))
    })
  }
  
  # Analyze 24h unique genes
  if (length(unique_genes_results$unique_24h) > 0) {
    message("\n  Analyzing 24h unique genes...")
    tryCatch({
      go_results_24h <- analyze_unique_genes_go(
        unique_genes_results$unique_24h,
        output_path = "analysis/GO_unique/24h"
      )
      go_results_all[["24h"]] <- go_results_24h
    }, error = function(e) {
      message(paste("    Warning: GO analysis failed for 24h:", e$message))
    })
  }
}
  # Compile results with new structure
  results <- list(
    merged_data = merged_df,
    changed_genes = changed_genes,
    unique_genes_6h = unique_genes_results$unique_6h,      # Separate 6h
    unique_genes_24h = unique_genes_results$unique_24h,    # Separate 24h
    gene_sets_6h = unique_genes_results$gene_sets_6h,      # Gene sets used for unique identification
    gene_sets_24h = unique_genes_results$gene_sets_24h,
 #  go_results = go_results_all,                          # commented out for now
    plots = list(upset = upset_results, venn = venn_plot)
  )
  
  # Generate cell line comparison plots
  message("\n=== Generating Cell Line Comparison Plots ===")
  comparison_results <- tryCatch({
    merged_file <- "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt"
    
    if (file.exists(merged_file)) {
      run_comparison_analysis(merged_file)
    } else {
      message("Merged file not found. Skipping comparison plots.")
      message("Please run the merge step first or check the file path.")
      NULL
    }
  }, error = function(e) {
    message(paste("Error generating comparison plots:", e$message))
    NULL
  })
  
  message("\n========================================")
  message("Comparative Analysis Complete!")
  message("========================================\n")
  
  return(results)
}

################################################################################
### SECTION 10: EXECUTE THE ANALYSIS
################################################################################

if (interactive()) {
  # Run the comparative analysis
  # Set debug = TRUE to see more detailed information
  # Optionally provide a mouse gene symbol conversion file
  comparative_results <- run_comparative_analysis(
    debug = TRUE,
    mouse_symbol_file = NULL,  # Set to path of JAX MGI file if available
    run_GO_analysis = TRUE    # set to false if you don't need to regen the whole GO analysis - it takes a long time to run.
  )
  
  # Run concordant/discordant analysis after main analysis completes
  message("\n=== Running Concordant/Discordant Analysis ===")
  concordant_discordant_results <- tryCatch({
    run_concordant_discordant_analysis()
  }, error = function(e) {
    message(paste("Error in concordant/discordant analysis:", e$message))
    NULL
  })
  
  # Print summary statistics for both time points
  message("\n=== Summary Statistics ===")
  
  # Summary for 6h
  message("\n--- 6h Time Point ---")
  for (cell_line in names(comparative_results$gene_sets_6h)) {
    n_changed_6h <- length(comparative_results$gene_sets_6h[[cell_line]])
    n_unique_6h <- 0
    if (cell_line %in% names(comparative_results$unique_genes_6h)) {
      n_unique_6h <- nrow(comparative_results$unique_genes_6h[[cell_line]])
    }
    message(paste(cell_line, ": ", n_changed_6h, " changed, ", n_unique_6h, " unique"))
  }
  
  # Summary for 24h
  message("\n--- 24h Time Point ---")
  for (cell_line in names(comparative_results$gene_sets_24h)) {
    n_changed_24h <- length(comparative_results$gene_sets_24h[[cell_line]])
    n_unique_24h <- 0
    if (cell_line %in% names(comparative_results$unique_genes_24h)) {
      n_unique_24h <- nrow(comparative_results$unique_genes_24h[[cell_line]])
    }
    message(paste(cell_line, ": ", n_changed_24h, " changed, ", n_unique_24h, " unique"))
  }
  
  # Combined summary (for reference with original pipeline)
  message("\n--- Combined (6h + 24h) ---")
  for (cell_line in names(comparative_results$changed_genes)) {
    n_changed <- length(comparative_results$changed_genes[[cell_line]])
    n_unique_6h <- 0
    n_unique_24h <- 0
    
    if (cell_line %in% names(comparative_results$unique_genes_6h)) {
      n_unique_6h <- nrow(comparative_results$unique_genes_6h[[cell_line]])
    }
    if (cell_line %in% names(comparative_results$unique_genes_24h)) {
      n_unique_24h <- nrow(comparative_results$unique_genes_24h[[cell_line]])
    }
    
    message(paste(cell_line, ": ", n_changed, " total changed, ",
                  n_unique_6h, " unique at 6h, ",
                  n_unique_24h, " unique at 24h"))
  }
  
  message("\n=== Session Information ===")
  writeLines(capture.output(sessionInfo()), "session_info/02_comparative_analysis_pipeline_sessioninfo.txt")
}