### library
library(openxlsx) # processing xlsx file
library(tidyverse) # Data_cleaning
library(limma) # alternative: library(DESeq2)
library(clusterProfiler) # converting geneid
library(org.Mm.eg.db)

### Read the Expression Matrix and Group Info
exp <- read.xlsx("./data/RawData_with_GroupInfo.xlsx", sheet = 1, rowNames = TRUE)
#group_info <- read.xlsx("./data/RawData_with_GroupInfo.xlsx", sheet = 2)
group_info <- read.xlsx("./data/RawData_with_GroupInfo.xlsx", sheet = 3)

### Data preprocessing
## log2 transformating
exp = log2(exp + 1) 
## remove low expression gene
exp = exp[rowMeans(exp) > 1, ] 

## Normalization
boxplot(exp,outline=FALSE, notch=F , las=2, main = "Before Normalization")
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=F , las=2, main = "After Normalization") # re-check the boxplot

### Differential Expression Analysis
## Prepare the Group Information for the Model and Design Matrix
group_list <- factor(group_info$Group,levels = c("Control","Ile"))
design <- model.matrix(~group_list) 
## Linear Regression and Empirical Bayes Moderation
v <- voom(exp, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
## Extract Differentially Expressed Genes (DEGs) and remove NAs
df_RNAs <- topTable(fit, coef = 2, number = Inf)
df_RNAs_cleans <- na.omit(df_RNAs) 

# Convert rownames to a new column named "Gene.ENSEMBL"
df_RNAs_cleans <- df_RNAs_cleans %>%
  rownames_to_column(var = "ENSEMBL")

### Convert to ENTREZID and Gene.Name
gene_symbols <- df_RNAs_cleans$ENSEMBL %>%
  bitr(fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db) %>%
  dplyr::rename("Gene.Name" = "SYMBOL")

gene_id <- df_RNAs_cleans$ENSEMBL %>%
  bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
  dplyr::rename("Gene.ID" = "ENTREZID")

df_RNAs_cleans <- merge(gene_id, df_RNAs_cleans, by = "ENSEMBL") %>%
  merge(gene_symbols, by = "ENSEMBL") %>%
  dplyr::select("ENSEMBL","Gene.Name","Gene.ID",
                "logFC", "P.Value", "adj.P.Val")

colnames(df_RNAs_cleans)

#Converting to 
write.xlsx(df_RNAs_cleans)

