#############################################################################
# RNA-seq Analysis Pipeline using DESeq2
# 
# This script performs differential expression analysis using DESeq2 with support
# for multiple group comparisons and visualization.
#
# Author: Kira Liu
# Date: 2025-02-10
#############################################################################

### Step 1: Load Required Packages ###########################################
#library(openxlsx)    # read.xlss For Excel file handling (recommended)
library(tidyverse)   # For data manipulation and visualization
library(DESeq2)      # For differential expression analysis
library(data.table)  # For efficient data reading (fread)


### Step 2: Load and Prepare Data ##########################################
# Read count matrix
count_matrix <- fread("./count_matrix.txt")

# format the count matrix
count_matrix <- data.frame(count_matrix)
rownames(count_matrix) <- count_matrix[,1]
count_matrix <- count_matrix[,-1]

# read.csv or read.xlsx are functions recommended to load files without formatting data frame
# count_matrix <- read.csv("./count_matrix.csv", row.names = TRUE)

### Step 3: Create Sample Information #####################################
# Read Sample Info
sample_info <- read.csv("./sample_info")

# Check sample counts before running analysis
unique(sample_info$Group)
table(sample_info$Group)


### Step 4: Define Analysis Parameters ###################################
# Define groups and comparisons
groups <- c("Control", "Experimental1", "Experimental2")

# Define comparisons
# Format: list of c("Experimental", "Control") pairs
comparisons <- list(
  c("Experimental1", "Control"),
  c("Experimental2", "Control")
  )


### Step 5: Run Analysis ################################################
# Option 1: Run analysis with gene annotation
results <- analyze_rnaseq_multi(
  count_matrix = count_matrix,
  sample_info = sample_info,
  groups = groups,
  comparisons = comparisons,
  output_dir = "RNAseq_DESeq2_results"
)

# Option 2: Run analysis without annotation
results <- analyze_rnaseq_simple(
  count_matrix = count_matrix,
  sample_info = sample_info,
  groups = groups,
  comparisons = comparisons,
  output_dir = "RNAseq_DESeq2_results"
)


### Step 6: Access Results #############################################
# The analysis creates several output files in your output directory:
# 1. DESeq2 results for each comparison (CSV files)
# 2. Normalized count matrices (VST and size-factor normalized)
# 3. Quality control plots in the 'plots' subdirectory:
#    - MA plots
#    - PCA plots (global and pairwise)
#    - Sample distance heatmaps
#    - Correlation heatmaps

# You can also access results programmatically:

## From annotated analysis (results_annotated):
# Access DESeq2 results with annotation
annotated_results <- results_annotated$annotated_results
# Access condition table
condition_table <- results_annotated$condition_table
# Access DESeq2 object for custom analysis
deseq_obj <- results_annotated$deseq_results$dds
# Access VST object for custom visualization
vst_obj <- results_annotated$deseq_results$vsd

## From simple analysis (results_simple):
# Access raw DESeq2 results
raw_results <- results_simple$raw_results
# Access normalized count data
vst_counts <- results_simple$vst_counts
size_norm_counts <- results_simple$size_normalized_counts
# Access DESeq2 object
deseq_obj_simple <- results_simple$deseq_object
# Access VST object
vst_obj_simple <- results_simple$vsd_object


### Optional: Example of Additional Analysis ###########################
# Example: Get significant genes for a specific comparison
sig_genes <- subset(raw_results[["Experimental1 vs Control"]], 
                    padj < 0.05 & abs(log2FoldChange) > 1)
