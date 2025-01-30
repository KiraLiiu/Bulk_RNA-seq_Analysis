

# RNA-seq Analysis Pipeline Functions
library(dplyr)
library(DESeq2)
library(AnnotationDbi)
library(pheatmap)

#' Create condition table for differential expression analysis
#' @param sample_info Sample metadata dataframe
#' @param id_col Column name containing sample IDs
#' @param group_col Column name containing group information
#' @param group1 Experimental group name
#' @param group2 Control group name
#' @param exp_label Label for experimental group
#' @param ctrl_label Label for control group
#' @return Dataframe with sample IDs and conditions
create_condition_table <- function(sample_info, id_col, group_col, group1, group2, 
                                   exp_label = "Experimental", ctrl_label = "Control") {
  if (!(id_col %in% colnames(sample_info))) stop(paste("Column", id_col, "not found in sample_info"))
  if (!(group_col %in% colnames(sample_info))) stop(paste("Column", group_col, "not found in sample_info"))
  
  sample_info_selected <- sample_info %>%
    filter(.data[[group_col]] %in% c(group1, group2)) %>%
    dplyr::select(all_of(id_col), all_of(group_col)) %>%
    mutate(condition = ifelse(.data[[group_col]] == group1, exp_label, ctrl_label),
           condition = factor(condition, levels = c(ctrl_label, exp_label))) %>%
    dplyr::select(all_of(id_col), condition)
  
  colnames(sample_info_selected)[1] <- "ID"
  return(sample_info_selected)
}

#' Prepare count matrix for DESeq2
#' @param count_matrix Raw count matrix
#' @param condition_table Condition table with sample IDs
#' @param min_count Minimum count threshold for filtering
#' @return List containing filtered count matrix and matching condition table
prepare_count_matrix <- function(count_matrix, condition_table, min_count = 10) {
  # Ensure count matrix contains all samples from condition table
  missing_samples <- setdiff(condition_table$ID, colnames(count_matrix))
  if (length(missing_samples) > 0) {
    stop("Missing samples in count matrix: ", paste(missing_samples, collapse = ", "))
  }
  
  # Filter and convert to integer
  filtered_matrix <- count_matrix %>%
    dplyr::select(all_of(condition_table$ID)) %>%
    mutate(across(everything(), as.integer))
  
  # Filter low-count genes
  keep_genes <- rowSums(filtered_matrix) > min_count
  filtered_matrix <- filtered_matrix[keep_genes, ]
  
  return(filtered_matrix)
}

#' Run DESeq2 analysis
#' @param count_matrix Filtered count matrix
#' @param condition_table Condition table
#' @param independent_filtering Logical, whether to perform independent filtering
#' @return DESeq2 results object
run_deseq2_analysis <- function(count_matrix, condition_table, independent_filtering = TRUE) {
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = condition_table,
    design = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results with optional independent filtering
  res <- results(dds, independentFiltering = independent_filtering)
  
  # Calculate normalized counts
  vsd <- vst(dds, blind = FALSE)
  
  return(list(
    dds = dds,
    results = res,
    vsd = vsd
  ))
}

#' Annotate DESeq2 results
#' @param deseq_results Results from run_deseq2_analysis
#' @param organism_db Organism annotation database
#' @param annotation_cols Columns to retrieve from annotation database
#' @return Annotated results dataframe
annotate_results <- function(deseq_results, 
                             organism_db = org.Mm.eg.db,
                             annotation_cols = c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE")) {
  
  # Convert results to dataframe
  res_df <- as.data.frame(deseq_results$results)
  res_df$ENSEMBL <- rownames(res_df)
  
  # Get normalized counts
  normalized_counts <- as.data.frame(assay(deseq_results$vsd))
  normalized_counts$ENSEMBL <- rownames(normalized_counts)
  
  # Get annotations
  gene_annotations <- AnnotationDbi::select(
    organism_db,
    keys = res_df$ENSEMBL,
    columns = annotation_cols,
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Merge all data
  annotated_res <- res_df %>%
    merge(gene_annotations, by = "ENSEMBL", all.x = TRUE) %>%
    merge(normalized_counts, by = "ENSEMBL", all.x = TRUE)
  
  return(annotated_res)
}

#' Generate visualization plots
#' @param deseq_results Results from run_deseq2_analysis
#' @param output_dir Directory for saving plots
#' @param prefix Prefix for plot filenames
#' @param width Plot width in inches
#' @param height Plot height in inches
generate_plots <- function(deseq_results, 
                           output_dir = "plots", 
                           prefix = "", 
                           width = 8, 
                           height = 6) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # MA plot
  pdf(file.path(output_dir, paste0(prefix, "MA_plot.pdf")), width = width, height = height)
  DESeq2::plotMA(deseq_results$results, ylim = c(-10, 10))
  dev.off()
  
  # Sample distance heatmap
  sample_distances <- dist(t(assay(deseq_results$vsd)))
  pdf(file.path(output_dir, paste0(prefix, "sample_distance_heatmap.pdf")), 
      width = width, height = height)
  pheatmap(as.matrix(sample_distances),
           clustering_distance_rows = sample_distances,
           clustering_distance_cols = sample_distances)
  dev.off()
  
  # PCA plot
  pdf(file.path(output_dir, paste0(prefix, "PCA_plot.pdf")), width = width, height = height)
  print(DESeq2::plotPCA(deseq_results$vsd, intgroup = "condition"))
  dev.off()
}

#' Main RNA-seq analysis function
#' @param count_matrix Pre-loaded count matrix
#' @param sample_info Pre-loaded sample information
#' @param group1 Experimental group name
#' @param group2 Control group name
#' @param id_col Column name containing sample IDs
#' @param group_col Column name containing group information
#' @param min_count Minimum count threshold for filtering
#' @param output_dir Directory for saving results
#' @param organism_db Organism annotation database
#' @return List containing analysis results
analyze_rnaseq <- function(count_matrix,
                           sample_info,
                           group1,
                           group2,
                           id_col = "ID",
                           group_col = "Group",
                           min_count = 10,
                           output_dir = "results",
                           organism_db = org.Mm.eg.db) {
  
  # Create condition table
  condition_table <- create_condition_table(
    sample_info, id_col, group_col, group1, group2
  )
  
  # Prepare count matrix
  filtered_counts <- prepare_count_matrix(
    count_matrix, condition_table, min_count
  )
  
  # Run DESeq2 analysis
  deseq_results <- run_deseq2_analysis(
    filtered_counts, condition_table
  )
  
  # Annotate results
  annotated_results <- annotate_results(
    deseq_results, organism_db
  )
  
  # Generate plots
  generate_plots(
    deseq_results,
    output_dir = file.path(output_dir, "plots"),
    prefix = paste0(group1, ".vs.", group2, "_")
  )
  
  # Save results to file
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(
    annotated_results,
    file = file.path(output_dir, paste0("DESeq2_", group1, ".vs.", group2, "_results.csv")),
    row.names = FALSE
  )
  
  return(list(
    deseq_results = deseq_results,
    annotated_results = annotated_results,
    condition_table = condition_table
  ))
}

### run experiments
# Load your data (example)
count_matrix <- read.xlsx("./data/data_yoyo.xlsx", sheet = 1, rowNames = T)
sample_info <- read.xlsx("./data/data_yoyo.xlsx", sheet = 3)

# Run analysis
results <- analyze_rnaseq(
  count_matrix = count_matrix,
  sample_info = sample_info,
  group1 = "R1CC",
  group2 = "WD1AA",
  output_dir = "my_analysis"
)

# Access results
deseq_results <- results$deseq_results        # Raw DESeq2 results
annotated_results <- results$annotated_results # Annotated results table
condition_table <- results$condition_table     # Sample grouping information
