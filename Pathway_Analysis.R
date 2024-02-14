## ----setup, echo=FALSE, results='hide', warning=FALSE, error=FALSE, message=FALSE, cache=FALSE----
library(knitr)
opts_chunk$set(
  cache = FALSE,
  echo = TRUE,
  warning = FALSE,
  error = FALSE,
  message = FALSE
)

##--------- Loading in Dataset & Setting index as Ensembl Gene IDs-----------
per_reps_counts = read.table('final_cleaned_all_reps.txt', sep = '\t', header = TRUE)
rep_names = c('par_rep1','par_rep2', 'per_rep1', 'per_rep2', 'per_rep3')
treat_factors = factor(c('control', 'control', 'treated', 'treated', 'treated'))
col_data = data.frame(row.names = rep_names, Treatment = treat_factors)
gene_expr_counts = per_reps_counts[ , rep_names] # Selected all columns except gene_id column
per_reps_counts = per_reps_counts[rowSums(gene_expr) >= 1, ] # Removed empty genes

##---- String Function to Remove Version IDs from Ensembl Gene IDs --------
clean_vers = function(my.string) {
  unlist(strsplit(my.string, split = "\\."))[1]
}
per_reps_counts$Gene_ID = sapply(per_reps_counts$Gene_ID,clean_vers)  # Stripped off gene version ID from Ensembl IDs

##----- Converting Ensembl IDs to HGNC symbols-------------
library(biomaRt)
library(magrittr)
library(dpylr)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
              values=per_reps_counts$Gene_ID, mart=mart)

matched = match(per_reps_counts$Gene_ID, genes$ensembl_gene_id)
per_reps_counts$HGNC = genes$hgnc_symbol[matched]
per_reps_counts$HGNC[per_reps_counts$HGNC == ""] <- NA # Replacing all empty strings as NA, since equivalent now

##---- Consolidating Genes with multiple HGNC symbols--------------
duplicated_genes = unique(per_reps_counts$HGNC[duplicated(per_reps_counts$HGNC)])
duplicated_genes
per_reps_empty = per_reps_counts[is.na(per_reps_counts$HGNC), ] # Genes with no HGNC symbol
per_reps_HGNC = per_reps_counts[is.na(per_reps_counts$HGNC) == FALSE, ] # Genes with HGNC symbol
# Combining read counts of genes with multiple Ensembl IDs
per_reps_uniqHGNC <- aggregate(. ~ HGNC, data = per_reps_HGNC[ , colnames(per_reps_HGNC)[2:7]], FUN = sum)
# Putting the dataframe with genes with no HGNC symbol and dataframe with genes with unique HGNC symbol together
cleaned_reps_counts <- rbind(per_reps_uniqHGNC, per_reps_empty[ , c("HGNC","par_rep1", "par_rep2", 
                                                                    "per_rep1", "per_rep2", "per_rep3")])

##---- RNA-seq Data Normalization----------------
library(DESeq2)
# import data to DESeq2 and setup for data normalization
dset = DESeqDataSetFromMatrix(cleaned_reps_counts[ , rep_names],
                              colData=col_data, design=~Treatment)
dset = DESeq(dset) # Normalization
gene_expr_VST = getVarianceStabilizedData(dset) # Applied Variance Stabilizing Transformation

