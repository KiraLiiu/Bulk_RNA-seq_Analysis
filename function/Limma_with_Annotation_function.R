library(openxlsx)
library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Mm.eg.db)

#### Define rowdata processing function using limma ####
perform_rowdata_processing <- function(exp, 
                                       group_info, 
                                       group_column,  # Custom column for group information
                                       group_levels = c("Control", "Experiment"),  # Custom group levels for comparison
                                       log_transform = TRUE, 
                                       gene_id_conversion = TRUE) {
  
  # Log2 Transformation
  if (log_transform) {
    exp <- log2(exp + 1)
  }
  
  # Remove Low-Expression Genes
  exp <- exp[rowMeans(exp) > 1, ]
  
  # Normalization
  boxplot(exp, outline = FALSE, notch = FALSE, las = 2, main = "Before Normalization")
  exp <- normalizeBetweenArrays(exp)
  boxplot(exp, outline = FALSE, notch = FALSE, las = 2, main = "After Normalization")
  
  # Prepare Group Information for the Model Matrix
  # Custom group column and levels provided
  group_list <- factor(group_info[[group_column]], levels = group_levels)
  design <- model.matrix(~group_list)
  
  # Differential Expression Analysis
  v <- voom(exp, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Extract DEGs
  df_RNAs <- topTable(fit, coef = 2, number = Inf)
  df_RNAs_cleans <- na.omit(df_RNAs)
  
  # Convert rownames to first column
  df_RNAs_cleans <- df_RNAs_cleans %>%
    rownames_to_column(var = "ENSEMBL")
  
  # Gene ID Conversion (ENSEMBL to SYMBOL and ENTREZID)
  if (gene_id_conversion) {
    gene_symbols <- df_RNAs_cleans$ENSEMBL %>%
      bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db) %>%
      dplyr::rename("Gene.Name" = "SYMBOL")
    
    gene_id <- df_RNAs_cleans$ENSEMBL %>%
      bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
      dplyr::rename("Gene.ID" = "ENTREZID")
    
    # Merge the results
    df_RNAs_cleans <- merge(gene_id, df_RNAs_cleans, by = "ENSEMBL") %>%
      merge(gene_symbols, by = "ENSEMBL") %>%
      dplyr::select("ENSEMBL", "Gene.Name", "Gene.ID", 
                    "logFC", "P.Value", "adj.P.Val")
  }
  
  # Return the final data frame
  return(df_RNAs_cleans)
}



#### Apply to raw data ####

### read raw data and group infomation
## for EXCEL files: preparing 1 xlsx file with 2 sheets:
# sheet1: raw data
# sheet2: groups information
rawdata_exp <- read.xlsx("{enter your file path}", sheet = 1, rowNames = TRUE)
group_info <- read.xlsx("{enter your file path}", sheet = 2)

## for CSV files: preparing 2 seperated files
# csv1: raw data
# csv2: group information
rawdata_exp <- read.csv("{enter your file path}",row.names = 1)
group_info <- read.csv("{enter your file path}",row.names = 1)

### input files to functions
## Reminder: modified the group information

result <- perform_rowdata_processing(exp = rawdata_exp, 
                                     group_info = group_info, 
                                     group_column = "Group",
                                     group_levels = c("Control", "Ile") # same to your group_info
                                     )

### output the xlsx or csv files
write.xlsx(result, "./Transcriptomics_Cleans.xlsx")
write.csv(result, "./Transcriptomics_Cleans.csv", row.names = FALSE)
