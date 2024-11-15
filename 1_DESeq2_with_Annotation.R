######################################
# RNA-seq Differential Expression Analysis with Annotation using DESeq2
######################################

# Author: Kira L
# Date: 2024-11-15
# Lab: Simcox Lab @ UW-Madison
# File Name: DESeq2_Analysis_with_Annotation.R

# Description:
# This script performs differential expression analysis on RNA-seq count data using the DESeq2 package.
# It includes data loading, sample filtering, DESeq2 analysis, gene annotation, and visualization of results.

# Workflow Outline:
# 1. Load and preprocess raw RNA-seq count data and sample metadata.
# 2. Filter low count genes to improve analysis robustness.
# 3. Set up DESeq2 dataset and perform differential expression analysis.
# 4. Extract DESeq2 results and apply variance-stabilizing transformation for visualization.
# 5. Annotate gene results with external databases.
# 6. Generate and save output files, including MA-plots, heatmaps, and PCA plots for sample clustering.

# Required Packages:
# openxlsx, dplyr, tidyverse, org.Mm.eg.db, DESeq2, AnnotationDbi, pheatmap

# Input Files:
# - data_RNAseq.xlsx: Contains raw RNA-seq counts (sheet 1) and sample information (sheet 3).

# Output Files:
# - DESeq2_results.xlsx: Differential expression results with gene annotations.
# - Visualizations (e.g., MA-plots, sample distance heatmap, PCA plot) saved in specified directories.

######################################

#### Loading Packages ####
# Load necessary libraries for analysis and data handling
library(openxlsx)         # For reading and writing Excel files
library(dplyr)            # For data manipulation
library(tidyverse)        # Collection of R packages for data science
library(org.Mm.eg.db)     # Mouse genome annotation database

# Step 1: Load RNA-seq raw count data and sample information
# Modify file paths as necessary
count_matrix <- read.xlsx("path/count_data.xlsx", sheet = 1, colNames = TRUE, rowNames = TRUE)
sample_info <- read.xlsx("path/sample_info.xlsx", sheet = 2)

# Step 2: Process sample metadata
# Create a `condition` column by combining experimental factors
# Example format:
# Sample    Condition
# A1        Treated
# A2        Treated
# B1        Untreated
sample_info <- sample_info %>%
  mutate(condition = paste(diet, temperature, sep = "_"))

# Step 3: Ensure count data is in integer format (required by DESeq2)
# Convert all count data columns to integer format
count_matrix <- count_matrix %>%
  mutate(across(everything(), as.integer))

### DESeq2 Analysis ###
library(DESeq2)  # For differential expression analysis

# Step 4: Create DESeq2 dataset
# Initialize DESeq2 dataset with count data and sample metadata
# 'condition' specifies the experimental groups for DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)

# Step 5: Filter out low-count genes to improve analysis robustness
# Retain genes with counts >10 in at least one sample to enhance statistical power
dds <- dds[rowSums(counts(dds)) > 10,]

# Step 6: Specify levels for 'condition' factor to set contrast direction
dds$condition <- factor(dds$condition, levels = c("Control_Group", "Treatment_Group"))

# Step 7: Run DESeq2 analysis
# Perform DESeq2 differential expression analysis based on the specified design
dds <- DESeq(dds)

# Step 8: Extract results for a specific contrast
# Retrieve results comparing experimental conditions with a specified significance threshold
contrast_names <- resultsNames(dds)
de_results <- results(dds, name = contrast_names[2], alpha = 0.05)  # Alpha sets the significance level
de_results_df <- as.data.frame(de_results)

# Step 9: Apply variance-stabilizing transformation (VST) for downstream analyses
# VST stabilizes variance across counts for improved visualization
vsd_transformed <- vst(dds, blind = FALSE)
normalized_counts <- as.data.frame(assay(vsd_transformed))

# Rename columns in `normalized_counts` using `condition` labels from `sample_info`
colnames(normalized_counts) <- sample_info$condition[match(colnames(normalized_counts), sample_info$sample)]
normalized_counts$ENSEMBL <- rownames(normalized_counts)  # Add gene ID as first column

### Gene Annotation ###
library(AnnotationDbi)  # For annotating gene information
library(org.Mm.eg.db)   # Mouse genome annotation database

# Retrieve relevant gene annotations using ENSEMBL IDs
de_results_df$ENSEMBL <- rownames(de_results_df)
gene_annotations <- AnnotationDbi::select(org.Mm.eg.db,
                                          keys = de_results_df$ENSEMBL,
                                          columns = c("SYMBOL", "ENTREZID", "GENENAME", "GENETYPE"),
                                          keytype = "ENSEMBL",
                                          multiVals = "first")

# Merge gene annotations with DESeq2 results and normalized data
annotated_de_results <- de_results_df %>%
  merge(gene_annotations, by = "ENSEMBL", all = TRUE) %>%
  merge(normalized_counts, by = "ENSEMBL", all = TRUE)

# Organize columns for improved readability
annotated_de_results <- annotated_de_results[, c("ENSEMBL", "SYMBOL", "ENTREZID", "GENENAME", "GENETYPE", 
                                                 "log2FoldChange", "pvalue", "padj", colnames(normalized_counts)[2:ncol(normalized_counts)])]

# Step 10: Export annotated DESeq2 results to an Excel file
write.xlsx(annotated_de_results, file = "path/DESeq2_results_annotated.xlsx", overwrite = TRUE)

# Step 11: Visualizations

# MA-plot to visualize differential expression
plotMA(de_results, ylim = c(-5, 5), main = "MA-plot of Differential Expression")

# Heatmap of sample distances using VST-transformed data
library(pheatmap)
sample_distances <- dist(t(assay(vsd_transformed)))
pheatmap(as.matrix(sample_distances), clustering_distance_rows = sample_distances, 
         clustering_distance_cols = sample_distances, main = "Sample Distance Heatmap")

# PCA plot for sample clustering based on condition
plotPCA(vsd_transformed, intgroup = "condition") + ggtitle("PCA Plot")

######################################
