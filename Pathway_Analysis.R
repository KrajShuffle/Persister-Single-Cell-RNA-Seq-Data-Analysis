## ----setup, echo=FALSE, results='hide', warning=FALSE, error=FALSE, message=FALSE, cache=FALSE----
library(knitr)
opts_chunk$set(
  cache = FALSE,
  echo = TRUE,
  warning = FALSE,
  error = FALSE,
  message = FALSE
)

##--------- Loading in Dataset & Setting index as Ensembl Gene IDs
per_reps_counts = read.table('final_cleaned_all_reps.txt', sep = '\t', header = TRUE)
rep_names = c('par_rep1','par_rep2', 'per_rep1', 'per_rep2', 'per_rep3')
treat_factors = factor(c('control', 'control', 'treated', 'treated', 'treated'))
col_data = data.frame(row.names = rep_names, Treatment = treat_factors)
gene_expr = per_reps_counts[ , rep_names] # Selected all columns except gene_id column
per_reps_counts = per_reps_counts[rowSums(gene_expr) >= 1, ] # Removed empty genes

##---- String Function to Remove Version IDs from Ensembl Gene IDs 
clean_vers = function(my.string) {
  unlist(strsplit(my.string, split = "\\."))[1]
}
per_reps_counts$Gene_ID = sapply(per_reps_counts$Gene_ID,clean_vers)
##-------- RNA-seq Data Processing Workflow
##--- Converting Ensembl IDs to HGNC symbols
library(biomaRt)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
              values=per_reps_counts$Gene_ID, mart=mart)

matched = match(per_reps_counts$Gene_ID, genes$ensembl_gene_id)
per_reps_counts$HGNC = genes$hgnc_symbol[matched]
duplicated_genes = unique(per_reps_counts$HGNC[duplicated(per_reps_counts$HGNC)])














library(DESeq2)
# import data to DESeq2 and setup for data normalization
dset = DESeqDataSetFromMatrix(per_reps_counts,
                              colData=col_data, design=~Treatment)
dset_filter_step1 = dset[rowSums(counts(dset)) >= 1, ] # Filter out genes with 0 reads across all replicates
dset_filtered_counts = counts(dset_filter_step1)
dset_filter_step1 = DESeq(dset_filter_step1) # Normalization
gene_expr = getVarianceStabilizedData(dset_filter_step1) # Applied Variance Stabilizing Transformation
rownames(gene_expr) = sapply(rownames(gene_expr),clean_vers) # Stripped off gene version ID for Normalized Data
rownames(dset_filtered_counts) = sapply(rownames(dset_filtered_counts), clean_vers) # Stripped off gene version ID for Count Data
dset_filtered_counts = as.data.frame(dset_filtered_counts)


