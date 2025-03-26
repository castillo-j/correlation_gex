library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)
library(ggpubr)

############ load object and name variables
se_obj<-readRDS('data_interim/luad_exp_maf_filtered.RDS')
# name of output file. will get saved in data_proocessed/
output<-'tcga_corr.csv'
# name the count assay in your SE object
count_assay<-'unstranded'
# name the column with patient id
patient_id<-'patient'
# name column with variable you want for correlation
variable<-'total_perMB'
###########################################

# Check if count_assay exists
if (!count_assay %in% names(assays(se_obj))) stop("Error: Count assay not found in se_obj")
# Check if metadata columns exist
if (!all(c(patient_id, variable) %in% colnames(colData(se_obj)))) stop("Error: One or more specified metadata columns are missing.")

### Create a normalized count assay to compare genes between samples
dds <- DESeqDataSetFromMatrix(countData = assay(se_obj, count_assay),
                              colData = as.data.frame(colData(se_obj)),
                              design = ~ 1)
dds<-estimateSizeFactors(dds)
assays(se_obj)$norm_count<-counts(dds, normalized=TRUE) # normalized counts
assays(se_obj)$vst_count<-assay(vst(dds,blind=F))       # for visualization

# Convert normalized counts to a dataframe
gex_table<-as.data.frame(t(assay(se_obj,'norm_count')))

# Extract relevant metadata
meta_table<-as.data.frame(colData(se_obj)[,c(patient_id, variable)])

# Merge metadata with expression data
gex_meta<-merge(meta_table, gex_table,by="row.names",all=TRUE)
rownames(gex_meta)<-gex_meta$Row.names
gex_meta[[variable]] <- as.numeric(gex_meta[[variable]])
gex_meta$Row.names<-NULL
rm(dds,gex_table,meta_table)

# Compute Spearman correlation coefficients
corr_table <- as.data.frame(sapply(gex_meta[-c(1:2)], function(x) cor(gex_meta[[variable]], x, method = "spearman")))
colnames(corr_table) <- "spearman"

# Compute p-values
sig_table <- as.data.frame(sapply(gex_meta[-c(1:2)], function(x) cor.test(gex_meta[[variable]], x, method = "spearman")$p.value))
colnames(sig_table) <- "p.value"

# Combine correlation and p-value tables
corr_table$p.value <- sig_table$p.value

# Adjust p-values using Benjamini-Hochberg correction
corr_table$p.adj <- p.adjust(corr_table$p.value, method = "BH")


# Add gene annotations
gene_info <- rowData(se_obj)[, c("gene_name", "gene_type", "gene_id")]
corr_table <- merge(gene_info, corr_table, by = "row.names", all.x = TRUE)
rownames(corr_table) <- corr_table$Row.names
corr_table$Row.names <- NULL

# Output dimensions and preview
head(corr_table)
write.csv(corr_table,paste0('data_processed/',output))






