############################################################
# Title: Differential Expression Analysis for E-MTAB-8031
# Author:Yadhukrishnan R
# Description:
#   This script downloads RNA-seq raw counts and metadata from
#   the E-MTAB-8031 dataset, processes the data, performs DESeq2
#   differential expression analysis, and extracts results for
#   selected genes.
############################################################

# Load required library
library(DESeq2)

# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-8031/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
View(counts)

# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-8031/resources/ExperimentDesignFile.RnaSeq/experiment-design")
View(metadata)

# Set row names and extract gene info for counts data
row.names(counts) = counts$Gene.ID
gene = counts[, c(1, 2)]
counts = counts[, -c(1,2)]

# Set row names for metadata
row.names(metadata) = metadata$Run

# Keep only necessary columns
metadata = metadata[,c("Sample.Characteristic.sampling.site.", "Sample.Characteristic.sex.", "Sample.Characteristic.individual.")]

#change column names
colnames(metadata) = c("sample", "sex", "patient")

# Siimplify factor names
metadata$sample[metadata$sample=="normal tissue adjacent to tumor"] = 'normal'
metadata$sample = factor(metadata$sample, levels = c("normal", "neoplasm"))
metadata$sex = factor(metadata$sex)
metadata$patient = factor(metadata$patient)

#Check if sample ordering matches between metadata and count matrix
all(rownames(metadata)==colnames(counts))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~sex+ sample)

#Filter low-count genes
dds <- dds[rowSums(counts(dds))>100, ]

#Run DESeq2
dds <- DESeq(dds)

#Get results for neoplasm vs normal
res = results(dds, contrast = c("sample", "neoplasm", "normal"), alpha = 0.05)

#Merge with gene info
res_df = as.data.frame(res)
res_df = merge(res_df, gene, by='row.names')

#Check selected genes
genes_check = c("GATA4", "HRH2", "CTLA4", "SUCNR1")
test = res_df[res_df$Gene.Name %in% genes_check, ]
#verify results
test
