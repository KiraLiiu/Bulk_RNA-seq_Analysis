### GO Enrichment: Function for integrating all categories with simple
integrate_GO_enrichment <- function(gene_ids, OrgDb, output_file, qvalue_cutoff = 0.05, pvalue_cutoff = 0.2) {
  library(clusterProfiler)
  library(dplyr)
  library(openxlsx)
  
  # Helper function to perform GO enrichment and simplify results
  enrich_and_simplify <- function(gene_ids, OrgDb, ont, qvalue_cutoff, pvalue_cutoff) {
    enrichment <- gene_ids %>%
      enrichGO(OrgDb = OrgDb, 
               ont = ont, 
               qvalueCutoff = qvalue_cutoff, 
               pvalueCutoff = pvalue_cutoff, 
               readable = TRUE, 
               pAdjustMethod = "BH") %>%
      simplify(cutoff = 0.5, by = "p.adjust", select_fun = min, measure = "Wang")
    
    # Add ontology category
    enrichment@result$Category <- ont
    return(enrichment@result)
  }
  
  # Perform GO enrichment for BP, MF, CC
  go_bp <- enrich_and_simplify(gene_ids, OrgDb, "BP", qvalue_cutoff, pvalue_cutoff)
  go_mf <- enrich_and_simplify(gene_ids, OrgDb, "MF", qvalue_cutoff, pvalue_cutoff)
  go_cc <- enrich_and_simplify(gene_ids, OrgDb, "CC", qvalue_cutoff, pvalue_cutoff)
  
  # Combine results
  go_combined <- rbind(go_bp, go_mf, go_cc)
  
  # Write to Excel
  write.xlsx(go_combined, output_file)
  
  # Return the result
  return(go_combined)
}

# Example usage:
# Perform the GO enrichment for a list of gene IDs and save the results to an Excel file
# result <- perform_GO_enrichment(gene_ids = union_deps_degs$Gene.ID, 
#                                 OrgDb = "org.Mm.eg.db", 
#                                 output_file = "./reports/GO_DEGs_DEPs_fc1.5_0.05_combined.xlsx")
