###############################################################################
#                     RNA-seq Analysis Pipeline Documentation                   #
#                                                                               #
# This pipeline performs differential expression analysis using DESeq2 with     #
# support for multiple group comparisons, visualizations, and gene annotation.  #
###############################################################################

# Required Packages:
# - dplyr: Data manipulation
# - DESeq2: Differential expression analysis
# - AnnotationDbi: Gene annotation tools
# - pheatmap: Heatmap visualization
# - tidyr: Data tidying
# - org.Mm.eg.db: Mouse genome database (for mouse data)
# - org.Hs.eg.db: Human genome database (for human data)
# Note: Choose the appropriate organism database for your data

# Load required libraries
library(dplyr)
library(DESeq2)
library(AnnotationDbi)
library(pheatmap)
library(tidyr)
library(org.Mm.eg.db)  # or org.Hs.eg.db for human data

#' Create condition table for multiple groups
#' @param sample_info Sample metadata dataframe
#' @param id_col Column name containing sample IDs
#' @param group_col Column name containing group information
#' @param groups Vector of group names to include
#' @return Dataframe with sample IDs and conditions
create_condition_table <- function(sample_info, id_col, group_col, groups) {
  if (!(id_col %in% colnames(sample_info))) stop(paste("Column", id_col, "not found in sample_info"))
  if (!(group_col %in% colnames(sample_info))) stop(paste("Column", group_col, "not found in sample_info"))
  
  sample_info_selected <- sample_info %>%
    filter(.data[[group_col]] %in% groups) %>%
    dplyr::select(all_of(id_col), all_of(group_col)) %>%
    mutate(condition = factor(.data[[group_col]], levels = groups)) %>%
    dplyr::select(all_of(id_col), condition)
  
  colnames(sample_info_selected)[1] <- "ID"
  return(sample_info_selected)
}

#' Validate sample groups and comparisons
#' @param condition_table Condition table with sample information
#' @param groups Vector of group names
#' @param comparisons List of pairwise comparisons
#' @return List with validation results and messages
validate_groups <- function(condition_table, groups, comparisons) {
  # Check sample counts for each group
  group_counts <- table(condition_table$condition)
  
  # Identify groups with no samples
  missing_groups <- groups[!groups %in% names(group_counts)]
  if (length(missing_groups) > 0) {
    stop("The following groups have no samples: ", 
         paste(missing_groups, collapse = ", "))
  }
  
  # Check sample counts
  low_sample_groups <- names(group_counts)[group_counts < 2]
  if (length(low_sample_groups) > 0) {
    warning("The following groups have less than 2 samples: ",
            paste(low_sample_groups, collapse = ", "))
  }
  
  # Validate comparisons
  invalid_comparisons <- list()
  valid_comparisons <- list()
  
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    
    # Check if both groups exist and have samples
    if (!all(c(group1, group2) %in% names(group_counts))) {
      invalid_comparisons[[length(invalid_comparisons) + 1]] <- comp
      next
    }
    
    # Check if both groups have enough samples
    if (group_counts[group1] < 2 || group_counts[group2] < 2) {
      warning("Comparison ", group1, " vs ", group2, 
              " has groups with less than 2 samples")
    }
    
    valid_comparisons[[length(valid_comparisons) + 1]] <- comp
  }
  
  if (length(invalid_comparisons) > 0) {
    stop("Invalid comparisons found: ", 
         paste(sapply(invalid_comparisons, paste, collapse=" vs "), collapse=", "))
  }
  
  # Print sample counts for each group
  message("Sample counts per group:")
  print(group_counts)
  
  return(valid_comparisons)
}


#' Prepare count matrix for DESeq2
#' @param count_matrix Raw count matrix
#' @param condition_table Condition table with sample IDs
#' @param min_count Minimum count threshold for filtering
#' @return Filtered count matrix
prepare_count_matrix <- function(count_matrix, condition_table, min_count = 10) {
  missing_samples <- setdiff(condition_table$ID, colnames(count_matrix))
  if (length(missing_samples) > 0) {
    stop("Missing samples in count matrix: ", paste(missing_samples, collapse = ", "))
  }
  
  filtered_matrix <- count_matrix %>%
    dplyr::select(all_of(condition_table$ID)) %>%
    mutate(across(everything(), as.integer))
  
  # Filter low-count genes
  keep_genes <- rowSums(filtered_matrix) > min_count
  filtered_matrix <- filtered_matrix[keep_genes, ]
  
  return(filtered_matrix)
}

#' Run DESeq2 analysis with multiple groups
#' @param count_matrix Filtered count matrix
#' @param condition_table Condition table
#' @param comparisons List of pairwise comparisons
#' @param independent_filtering Logical, whether to perform independent filtering
#' @return DESeq2 results object
run_deseq2_analysis <- function(count_matrix, condition_table, comparisons, 
                                independent_filtering = TRUE) {
  # Create and run DESeq2 object for all groups
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = condition_table,
    design = ~ condition
  )
  
  # Run DESeq2 once for global normalization
  dds <- DESeq(dds)
  
  # Get results for each comparison
  results_list <- list()
  for (comp in comparisons) {
    results_list[[paste(comp[1], "vs", comp[2])]] <- results(
      dds,
      contrast = c("condition", comp[1], comp[2]),
      independentFiltering = independent_filtering
    )
  }
  
  # Calculate normalized counts
  vsd <- vst(dds, blind = FALSE)
  
  return(list(
    dds = dds,
    results = results_list,
    vsd = vsd
  ))
}



#' Get available organism databases
#' @return Character vector of available organism databases
list_organism_dbs <- function() {
  return(c(
    "org.Hs.eg.db",  # Human
    "org.Mm.eg.db",  # Mouse
    "org.Rn.eg.db",  # Rat
    "org.Dm.eg.db",  # Fruit fly
    "org.Ce.eg.db",  # C. elegans
    "org.Sc.sgd.db", # Yeast
    "org.Dr.eg.db"   # Zebrafish
  ))
}

#' Load organism database
#' @param organism String specifying organism database or database object
#' @return Loaded organism database object
load_organism_db <- function(organism) {
  if (is.character(organism)) {
    if (!requireNamespace(organism, quietly = TRUE)) {
      stop("Package ", organism, " is not installed. Please install it using:\n",
           "BiocManager::install('", organism, "')")
    }
    db <- eval(parse(text = paste0(organism, "::", organism)))
    return(db)
  } else {
    return(organism)
  }
}




#' Annotate DESeq2 results for multiple comparisons
#' @param deseq_results Results from run_deseq2_analysis
#' @param organism String specifying organism database or database object
#' @param annotation_cols Columns to retrieve from annotation database
#' @param keytype Type of gene IDs in count matrix
#' @param remove_duplicates Logical, whether to remove duplicate gene symbols
#' @param remove_na_symbol Logical, whether to remove NA gene symbols
#' @param duplicate_agg Function for aggregating duplicates
#' @param report_stats Logical, whether to print annotation statistics
#' @param multiVals How to handle multiple values from annotation
#' @return List of annotated results dataframes
annotate_results <- function(deseq_results, 
                             organism = "org.Mm.eg.db",
                             annotation_cols = c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE"),
                             keytype = "ENSEMBL",
                             remove_duplicates = TRUE,
                             remove_na_symbol = TRUE,
                             duplicate_agg = NULL,
                             report_stats = TRUE,
                             multiVals = "first") {
  
  # Load organism database
  organism_db <- load_organism_db(organism)
  
  # Get normalized counts
  normalized_counts <- as.data.frame(assay(deseq_results$vsd))
  normalized_counts$ENSEMBL <- rownames(normalized_counts)
  
  # Get annotations
  if(report_stats) message("Fetching annotations from database...")
  gene_annotations <- AnnotationDbi::select(
    organism_db,
    keys = rownames(normalized_counts),
    columns = annotation_cols,
    keytype = keytype,
    multiVals = multiVals
  )
  
  # Report initial statistics
  if(report_stats) {
    message("\nInitial annotation statistics:")
    message("Total genes: ", length(unique(gene_annotations[[keytype]])))
    message("NA symbols: ", sum(is.na(gene_annotations$SYMBOL)))
    message("Duplicate symbols: ", sum(duplicated(gene_annotations$SYMBOL[!is.na(gene_annotations$SYMBOL)])))
  }
  
  # Handle NA symbols
  if(remove_na_symbol) {
    gene_annotations <- gene_annotations[!is.na(gene_annotations$SYMBOL), ]
    if(report_stats) message("Removed entries with NA symbols")
  }
  
  # Handle duplicates
  if(remove_duplicates) {
    if(!is.null(duplicate_agg)) {
      # Aggregate duplicates using provided function
      if(report_stats) message("Aggregating duplicate symbols using ", substitute(duplicate_agg))
      
      gene_annotations <- gene_annotations %>%
        group_by(SYMBOL) %>%
        summarise(across(where(is.numeric), duplicate_agg, na.rm = TRUE),
                  across(where(is.character), first),
                  .groups = "drop")
    } else {
      # Keep first occurrence of each symbol
      if(report_stats) message("Keeping first occurrence of duplicate symbols")
      gene_annotations <- gene_annotations[!duplicated(gene_annotations$SYMBOL), ]
    }
  }
  
  # Process each comparison
  annotated_results <- list()
  for (comparison_name in names(deseq_results$results)) {
    if(report_stats) message("\nProcessing comparison: ", comparison_name)
    
    # Convert DESeq2 results to data frame
    res_df <- as.data.frame(deseq_results$results[[comparison_name]])
    res_df$ENSEMBL <- rownames(res_df)
    
    # Merge annotations and normalized counts
    annotated_res <- res_df %>%
      merge(gene_annotations, by = keytype, all.x = TRUE) %>%
      merge(normalized_counts, by = keytype, all.x = TRUE)
    
    if(report_stats) {
      # Report annotation coverage
      total_genes <- nrow(res_df)
      annotated_genes <- sum(!is.na(annotated_res$SYMBOL))
      message("Annotation coverage: ", 
              round(annotated_genes/total_genes * 100, 2), 
              "% (", annotated_genes, "/", total_genes, " genes)")
    }
    
    annotated_results[[comparison_name]] <- annotated_res
  }
  
  return(annotated_results)
}


#' Generate PCA plot for a subset of samples
#' @param vst_counts VST-transformed count data
#' @param conditions Sample conditions
#' @param title Plot title
#' @return ggplot object
create_pca_plot <- function(vst_counts, conditions, title) {
  # Transpose the data for PCA
  pca_data <- t(vst_counts)
  
  # Calculate variances for each gene
  gene_vars <- apply(pca_data, 2, var)
  
  # Filter out genes with zero or very low variance
  # Keep genes with variance above a small threshold
  var_threshold <- 1e-10
  keep_genes <- gene_vars > var_threshold
  
  if (sum(keep_genes) < 2) {
    warning("Not enough variable genes for PCA in comparison: ", title)
    return(NULL)
  }
  
  # Subset the data to keep only variable genes
  pca_data <- pca_data[, keep_genes]
  
  # Perform PCA
  pca_result <- prcomp(pca_data, scale. = TRUE)
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
  
  # Create data frame for plotting
  plot_data <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2],
    Condition = conditions
  )
  
  # Create plot
  ggplot(plot_data, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3) +
    stat_ellipse(type = "t", level = 0.95) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
}



#' Generate visualization plots for multiple groups with subgroup analysis
#' @param deseq_results Results from run_deseq2_analysis
#' @param output_dir Directory for saving plots
#' @param comparisons List of pairwise comparisons
#' @param width Plot width in inches
#' @param height Plot height in inches
generate_plots <- function(deseq_results, output_dir = "plots", comparisons, width = 8, height = 6) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get all samples and conditions
  sample_data <- colData(deseq_results$dds)
  all_groups <- levels(sample_data$condition)
  
  # 1. MA plots for each comparison
  for (comparison_name in names(deseq_results$results)) {
    pdf(file.path(output_dir, paste0("MA_plot_", comparison_name, ".pdf")), 
        width = width, height = height)
    DESeq2::plotMA(deseq_results$results[[comparison_name]], ylim = c(-10, 10))
    dev.off()
  }
  
  # 2. Global plots
  # Select top variable genes for global plots
  vst_assay <- assay(deseq_results$vsd)
  gene_vars <- rowVars(vst_assay)
  top_var_genes <- order(gene_vars, decreasing = TRUE)[1:min(500, length(gene_vars))]
  
  # Global sample distance heatmap
  sample_distances <- dist(t(assay(deseq_results$vsd)))
  pdf(file.path(output_dir, "global_sample_distance_heatmap.pdf"), 
      width = width, height = height)
  pheatmap(as.matrix(sample_distances),
           clustering_distance_rows = sample_distances,
           clustering_distance_cols = sample_distances,
           main = "Global Sample Distance Heatmap",
           annotation_col = data.frame(Condition = sample_data$condition,
                                       row.names = rownames(sample_data)))
  dev.off()
  
  # Global PCA plot
  pdf(file.path(output_dir, "global_PCA_plot.pdf"), width = width, height = height)
  print(DESeq2::plotPCA(deseq_results$vsd, intgroup = "condition") +
          ggtitle("Global PCA Plot"))
  dev.off()
  
  # 3. Pairwise comparison plots
  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]
    
    # Subset VST data for these two groups
    idx <- sample_data$condition %in% c(group1, group2)
    if (sum(idx) < 2) next  # Skip if not enough samples
    
    # Get subset of VST data
    vst_counts <- vst_assay[, idx]
    subset_conditions <- sample_data$condition[idx]
    
    # Select top variable genes for this comparison
    gene_vars <- rowVars(vst_counts)
    top_var_genes <- order(gene_vars, decreasing = TRUE)[1:min(500, length(gene_vars))]
    vst_counts_filtered <- vst_counts[top_var_genes,]
    
    # Pairwise PCA plot
    pca_plot <- create_pca_plot(
      vst_counts_filtered,
      subset_conditions,
      paste("PCA Plot:", group1, "vs", group2)
    )
    
    if (!is.null(pca_plot)) {
      pdf(file.path(output_dir, paste0("PCA_", group1, "_vs_", group2, ".pdf")), 
          width = width, height = height)
      print(pca_plot)
      dev.off()
    }
    
    # Pairwise sample distance heatmap
    sample_distances_pair <- dist(t(vst_counts_filtered))
    pdf(file.path(output_dir, paste0("heatmap_", group1, "_vs_", group2, ".pdf")), 
        width = width, height = height)
    pheatmap(as.matrix(sample_distances_pair),
             clustering_distance_rows = sample_distances_pair,
             clustering_distance_cols = sample_distances_pair,
             main = paste("Sample Distance Heatmap:", group1, "vs", group2),
             annotation_col = data.frame(Condition = subset_conditions,
                                         row.names = colnames(vst_counts)))
    dev.off()
    
    # Pairwise correlation heatmap
    cor_matrix <- cor(vst_counts_filtered)
    pdf(file.path(output_dir, paste0("correlation_", group1, "_vs_", group2, ".pdf")), 
        width = width, height = height)
    pheatmap(cor_matrix,
             main = paste("Sample Correlation:", group1, "vs", group2),
             annotation_col = data.frame(Condition = subset_conditions,
                                         row.names = colnames(vst_counts)))
    dev.off()
  }
  
  # 4. Group-specific correlation heatmaps
  for (group in all_groups) {
    idx <- sample_data$condition == group
    if (sum(idx) < 2) next
    
    group_data <- vst_assay[top_var_genes, idx]
    cor_matrix <- cor(group_data, method = "pearson")
    
    pdf(file.path(output_dir, paste0("correlation_heatmap_", group, ".pdf")), 
        width = width, height = height)
    pheatmap(cor_matrix,
             main = paste("Sample Correlation Heatmap:", group),
             display_numbers = TRUE,
             number_format = "%.2f",
             fontsize_number = 8)
    dev.off()
  }
  
  # 5. Global sample correlation heatmap
  sample_cor <- cor(vst_assay[top_var_genes,])
  pdf(file.path(output_dir, "global_correlation_heatmap.pdf"), 
      width = width, height = height)
  pheatmap(sample_cor,
           main = "Global Sample Correlation Heatmap",
           annotation_col = data.frame(Condition = sample_data$condition,
                                       row.names = rownames(sample_data)),
           annotation_row = data.frame(Condition = sample_data$condition,
                                       row.names = rownames(sample_data)))
  dev.off()
}


#' Run DESeq2 analysis with multiple groups (updated with validation)
#' @param count_matrix Filtered count matrix
#' @param condition_table Condition table
#' @param comparisons List of pairwise comparisons
#' @param independent_filtering Logical, whether to perform independent filtering
#' @return DESeq2 results object
run_deseq2_analysis <- function(count_matrix, 
                                condition_table, 
                                comparisons, 
                                independent_filtering = TRUE) {
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = condition_table,
    design = ~ condition
  )
  
  # Run DESeq2
  tryCatch({
    dds <- DESeq(dds)
  }, error = function(e) {
    if (grepl("the model matrix is not full rank", e$message)) {
      # Get actual sample counts
      sample_counts <- table(condition_table$condition)
      stop("DESeq2 analysis failed due to insufficient samples.\n",
           "Sample counts per group:\n",
           paste(capture.output(print(sample_counts)), collapse="\n"), "\n",
           "Each group should have at least 2 samples.")
    } else {
      stop(e)
    }
  })
  
  # Get results for each comparison
  results_list <- list()
  for (comp in comparisons) {
    tryCatch({
      results_list[[paste(comp[1], "vs", comp[2])]] <- results(
        dds,
        contrast = c("condition", comp[1], comp[2]),
        independentFiltering = independent_filtering
      )
    }, error = function(e) {
      warning("Failed to get results for comparison ", comp[1], " vs ", comp[2], ": ", e$message)
    })
  }
  
  # Calculate normalized counts
  vsd <- vst(dds, blind = FALSE)
  
  return(list(
    dds = dds,
    results = results_list,
    vsd = vsd
  ))
}



#' Extract raw DESeq2 results and normalized counts
#' @param deseq_results Results from DESeq2 analysis
#' @param output_dir Directory for saving results
#' @return List of raw results and normalized counts
get_raw_results <- function(deseq_results, output_dir) {
  # Get normalized counts (both VST and size-factor normalized)
  vst_counts <- as.data.frame(assay(deseq_results$vsd))
  size_normalized <- as.data.frame(counts(deseq_results$dds, normalized=TRUE))
  
  # Process DESeq2 results for each comparison
  raw_results <- list()
  for (comparison_name in names(deseq_results$results)) {
    # Convert DESeq2 results to data frame
    res_df <- as.data.frame(deseq_results$results[[comparison_name]])
    res_df$gene_id <- rownames(res_df)
    
    # Reorder columns to put gene_id first
    res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
    
    # Save to file
    write.csv(res_df,
              file = file.path(output_dir, paste0("DESeq2_", comparison_name, "_results.csv")),
              row.names = FALSE)
    
    raw_results[[comparison_name]] <- res_df
  }
  
  # Save normalized counts
  write.csv(vst_counts,
            file = file.path(output_dir, "vst_normalized_counts.csv"),
            row.names = TRUE)
  write.csv(size_normalized,
            file = file.path(output_dir, "size_factor_normalized_counts.csv"),
            row.names = TRUE)
  
  return(list(
    raw_results = raw_results,
    vst_counts = vst_counts,
    size_normalized_counts = size_normalized
  ))
}


#' Main RNA-seq analysis function with annotation
#' @param count_matrix Pre-loaded count matrix
#' @param sample_info Pre-loaded sample information
#' @param groups Vector of group names to include
#' @param comparisons List of pairwise comparisons
#' @param id_col Column name containing sample IDs
#' @param group_col Column name containing group information
#' @param min_count Minimum count threshold for filtering
#' @param output_dir Directory for saving results
#' @param organism String specifying organism database or database object
#' @param keytype Type of gene IDs in count matrix
#' @param annotation_cols Columns to retrieve from annotation database
#' @param remove_duplicates Logical, whether to remove duplicate gene symbols
#' @param remove_na_symbol Logical, whether to remove NA gene symbols
#' @param duplicate_agg Function for aggregating duplicates
#' @param report_stats Logical, whether to print annotation statistics
#' @return List containing analysis results
analyze_rnaseq_multi <- function(count_matrix,
                                 sample_info,
                                 groups,
                                 comparisons,
                                 id_col = "ID",
                                 group_col = "Group",
                                 min_count = 10,
                                 output_dir = "results",
                                 organism = "org.Mm.eg.db",
                                 keytype = "ENSEMBL",
                                 annotation_cols = c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE"),
                                 remove_duplicates = TRUE,
                                 remove_na_symbol = TRUE,
                                 duplicate_agg = NULL,
                                 report_stats = TRUE) {
  
  # Create condition table for all groups
  condition_table <- create_condition_table(
    sample_info, id_col, group_col, groups
  )
  
  # Validate groups and comparisons
  message("Validating sample groups and comparisons...")
  valid_comparisons <- validate_groups(condition_table, groups, comparisons)
  
  # Print sample information
  message("\nSample information:")
  print(table(condition_table$condition))
  
  # Prepare count matrix
  message("\nPreparing count matrix...")
  filtered_counts <- prepare_count_matrix(
    count_matrix, condition_table, min_count
  )
  
  # Run DESeq2 analysis
  message("\nRunning DESeq2 analysis...")
  deseq_results <- run_deseq2_analysis(
    filtered_counts, condition_table, valid_comparisons
  )
  
  # Annotate results with enhanced options
  message("\nAnnotating results...")
  annotated_results <- annotate_results(
    deseq_results = deseq_results,
    organism = organism,
    annotation_cols = annotation_cols,
    keytype = keytype,
    remove_duplicates = remove_duplicates,
    remove_na_symbol = remove_na_symbol,
    duplicate_agg = duplicate_agg,
    report_stats = report_stats
  )
  
  # Generate plots
  message("\nGenerating plots...")
  generate_plots(
    deseq_results,
    output_dir = file.path(output_dir, "plots"),
    comparisons = valid_comparisons
  )
  
  # Save results
  message("\nSaving results...")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (comparison_name in names(annotated_results)) {
    write.csv(
      annotated_results[[comparison_name]],
      file = file.path(output_dir, paste0("DESeq2_", comparison_name, "_results.csv")),
      row.names = FALSE
    )
  }
  
  message("\nAnalysis complete!")
  return(list(
    deseq_results = deseq_results,
    annotated_results = annotated_results,
    condition_table = condition_table
  ))
}

#' Simple DESeq2 analysis function without annotation
#' @param count_matrix Pre-loaded count matrix
#' @param sample_info Pre-loaded sample information
#' @param groups Vector of group names to include
#' @param comparisons List of pairwise comparisons
#' @param id_col Column name containing sample IDs
#' @param group_col Column name containing group information
#' @param min_count Minimum count threshold for filtering
#' @param output_dir Directory for saving results
#' @param save_tables Logical, whether to save results to files
#' @return List containing analysis results
analyze_rnaseq_simple <- function(count_matrix,
                                  sample_info,
                                  groups,
                                  comparisons,
                                  id_col = "ID",
                                  group_col = "Group",
                                  min_count = 10,
                                  output_dir = "results",
                                  save_tables = TRUE) {
  
  # Create condition table for all groups
  condition_table <- create_condition_table(
    sample_info, id_col, group_col, groups
  )
  
  # Validate groups and comparisons
  message("Validating sample groups and comparisons...")
  valid_comparisons <- validate_groups(condition_table, groups, comparisons)
  
  # Print sample information
  message("\nSample information:")
  print(table(condition_table$condition))
  
  # Prepare count matrix
  message("\nPreparing count matrix...")
  filtered_counts <- prepare_count_matrix(
    count_matrix, condition_table, min_count
  )
  
  # Run DESeq2 analysis
  message("\nRunning DESeq2 analysis...")
  deseq_results <- run_deseq2_analysis(
    filtered_counts, condition_table, valid_comparisons
  )
  
  # Generate plots
  message("\nGenerating plots...")
  generate_plots(
    deseq_results,
    output_dir = file.path(output_dir, "plots"),
    comparisons = valid_comparisons
  )
  
  # Get and save raw results
  message("\nExtracting and saving results...")
  if (save_tables) {
    raw_results <- get_raw_results(deseq_results, output_dir)
  } else {
    raw_results <- list(
      raw_results = lapply(deseq_results$results, as.data.frame),
      vst_counts = as.data.frame(assay(deseq_results$vsd)),
      size_normalized_counts = as.data.frame(counts(deseq_results$dds, normalized=TRUE))
    )
  }
  
  message("\nAnalysis complete!")
  return(list(
    deseq_object = deseq_results$dds,
    vsd_object = deseq_results$vsd,
    raw_results = raw_results$raw_results,
    vst_counts = raw_results$vst_counts,
    size_normalized_counts = raw_results$size_normalized_counts,
    condition_table = condition_table
  ))
}
