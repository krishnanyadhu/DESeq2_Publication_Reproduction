# Load required library
library(DESeq2)

# Load metadata and raw count data from EBI Expression Atlas
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10041/resources/ExperimentDesignFile.RnaSeq/experiment-design")
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-10041/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")

# Set gene IDs as row names for count matrix
# Keep a copy of gene annotation columns separately
row.names(counts) = counts$Gene.ID
genes = counts[, c(1,2)]
counts = counts[, -c(1,2)]

# Wrangle metadata: keep only relevant columns and rename them
row.names(metadata) = metadata$Run
metadata = metadata[, c("Factor.Value.compound.","Factor.Value.phenotype.")]
colnames(metadata) = c("compound", "genotype")

# Simplify and standardize factor names in metadata
metadata$compound[metadata$compound == "Nutlin-3a 150 milligram per kilogram per day"] = "nutlin"
metadata$genotype[metadata$genotype == "p53 heterozygous mutant"] = "p53_hetero"
metadata$genotype[metadata$genotype == "p53R248Q heterozygous mutant"] = "R248Q"
metadata$genotype[metadata$genotype == "wild type phenotype"] = "wildtype"

# Remove specific runs from metadata (outliers or unwanted samples)
metadata <- metadata[!rownames(metadata) %in% c("ERR5216741", "ERR5216742", "ERR5216743", "ERR5216738"), ]

# Create a combined factor variable "group" for genotype + compound
metadata$group = factor(paste(metadata$genotype, metadata$compound))

# Check if sample ordering matches between metadata and count matrix
all(row.names(metadata)==colnames(counts))

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~group)

# Filter out genes with low counts across all samples
dds <- dds[rowSums(counts(dds))>10, ]

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Extract results for specific contrasts
# Contrast 1: p53_hetero nutlin vs p53_hetero none
res1 <- results(dds,contrast = c("group", "p53_hetero nutlin", "p53_hetero none") ,alpha = 0.05)

# Contrast 2: R248Q nutlin vs R248Q none
res2 <- results(dds,contrast = c("group", "R248Q nutlin", "R248Q none") ,alpha = 0.05)

#Convert results to data frames and merge with gene annotations
res1_df = as.data.frame(res1)
res1_df = merge(res1_df, genes, by="row.names")

res2_df = as.data.frame(res2)
res2_df = merge(res2_df, genes, by="row.names")

# Retrieve results for specific genes of interest to verify with published results
res2_df[res2_df$Gene.Name == "Pitx2", ]
res1_df[res1_df$Gene.Name == "Gzmk", ]
