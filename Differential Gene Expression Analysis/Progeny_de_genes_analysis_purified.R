##--------- Reading in as dataframe -----
per_reps_counts = read.table('final_cleaned_all_reps.txt', sep = '\t', header = TRUE)
row_names = c('par_rep1','par_rep2', 'per_rep1', 'per_rep2', 'per_rep3')
treat_factors = factor(c('control', 'control', 'treated', 'treated', 'treated'))
col_data = data.frame(row.names = row_names, Treatment = treat_factors)
rownames(per_reps_counts) = per_reps_counts[ , 1]
per_reps_counts[ , 1] = NULL

##---- Remove version number from Ensembl Ids ----------------------------------
clean_vers = function(my.string) {
  unlist(strsplit(my.string, split = "\\."))[1]
}
rownames(per_reps_counts) = sapply(rownames(per_reps_counts),clean_vers)

## Pre-screening and loading into DESeq2----------------------------------------
library(DESeq2)
# Convert data into DESeq2 object
dset = DESeqDataSetFromMatrix(per_reps_counts,
    colData=col_data, design=~Treatment)

##--- Filtering out low read count genes----------------------------------------
dset_filter_step1 = dset[rowSums(counts(dset)) >= 1, ]
dset_filtered_counts = counts(dset_filter_step1)
check = counts(dset_filter_step1)
percent_filter =  1- (nrow(check) / nrow(per_reps_counts))
dset_filter_step1 = DESeq(dset_filter_step1)

# Checking results -------------------------------------------------------------
plotDispEsts(dset_filter_step1)
# Defining and Running Differential Gene Analysis
contrast_k = c("Treatment, treated", "control")
res_table_per = results(dset_filter_step1,contrast = c('Treatment', 'treated', 'control'), alpha = 0.05)
mcols(res_table_per, use.names=TRUE) # generates list of columns with info on each
sum(res_table_per$padj < 0.05, na.rm=TRUE ) # determines number of genes with p-value less than 0.05
summary(res_table_per)
resSig = res_table_per[which(res_table_per$padj < 0.05), ]
head( resSig[ order( resSig$log2FoldChange ), ] )
tail( resSig[ order( resSig$log2FoldChange ), ] )
# Adding in gene names
res_table_per$ensembl <- sapply( strsplit(rownames(res_table_per), split="\\+" ), "[", 1 )

# annotate matrix with HGNC symbols
# BiocManager::install("biomaRt")
library(biomaRt)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
              values=res_table_per$ensembl, mart=mart)
idx <- match( res_table_per$ensembl, genes$ensembl_gene_id )
res_table_per$hgnc_symbol = genes$hgnc_symbol[ idx ]
head(res_table_per,4)
#write.csv(as.data.frame(res_table_per), file = "C:\\Users\\Karthik\\Desktop\\CHECK_up_down_de_genes.csv")

# Removal of NaNs and Empty Values from HGNC column-----------------------------
df_de_analysis = as.data.frame(res_table_per)
non_na_tf = !is.na(df_de_analysis$hgnc_symbol) # Setting up non-NaNs slice
hgnc_wo_na = df_de_analysis[non_na_tf, ] # No NaNs in table
k  = !hgnc_wo_na$hgnc_symbol == "" # Setting up removal of empty values in Hgnc column
hgnc_clean = hgnc_wo_na[!(hgnc_wo_na$hgnc_symbol == ""), ] # Removed both nans and empty values

# Checking for duplicates and removing them-------------------------------------
value_counts_hgnc = table(hgnc_clean$hgnc_symbol) # Which are the duplicates
dupes = c('ITFG2-AS1', 'LINC01238', 'POLR2J4', 'RN7SL274P', 'TUBB7P')
t_f_dupes = hgnc_clean$hgnc_symbol %in% dupes
hgnc_dupes = hgnc_clean[t_f_dupes, ]
dupes_ensembl = c('ENSG00000263947', 'ENSG00000261186', 'ENSG00000272655', 'ENSG00000258325', 'ENSG00000173876') #looked through gene cards to determine which ones to remove
t_f_ensembl = !(hgnc_clean$ensembl %in% dupes_ensembl)
hgnc_all_clean = hgnc_clean[t_f_ensembl, ]
rownames(hgnc_all_clean) = hgnc_all_clean$hgnc_symbol

# In case BioMart versions differ, below code will result in dataframe without dupes -----
#tf_rem_dupes = !duplicated(hgnc_clean$hgnc_symbol)
#no_dupes_hgnc = hgnc_clean[tf_rem_dupes, ] 


# Purified DE genes table (Should work if uncommented code worked)--------------
noNA_padj_genes = hgnc_all_clean[!is.na(hgnc_all_clean$padj), ]
de_genes = noNA_padj_genes[noNA_padj_genes$padj < 0.05, ] # Differential Genes Dataframe
de_genes_sorted = de_genes[order(de_genes$padj), ]
View(de_genes_sorted) # Sorts the adjusted p-values in ascending order and enables viewing

# If BioMart versions differ, same cleaning as above -----------
# noNA_padj_genes = no_dupes_hgnc[!is.na(no_dupes_hgnc$padj), ]
# de_genes = noNA_padj_genes[noNA_padj_genes$padj < 0.05, ] # Differential Genes Dataframe
# de_genes_sorted = de_genes[order(de_genes$padj), ]
# View(de_genes_sorted) # Sorts the adjusted p-values in ascending order and enables viewing


# Volcano Plots -------------------------
BiocManager::install("EnhancedVolcano") # Installs EnhancedVolcano package
library(EnhancedVolcano) 

EnhancedVolcano(hgnc_all_clean, 
                lab = rownames(hgnc_all_clean),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = 'Persister Cell Dataset')
# In Case Biomart versions differ,--------------
# EnhancedVolcano(no_dupes_hgnc, lab = rownames(no_dupes_hgnc), x = 'log2FoldChange', y = 'padj',
# pCutoff = 0.05,
# FCcutoff = 1.5,
# title = 'Persister Cell Dataset')



# Export as Image and Size Note --------------
# To view plots, Export as image 
# with width: 750 and height: 950
## ----echo=FALSE---------------------------------------------------------------
sessionInfo()
