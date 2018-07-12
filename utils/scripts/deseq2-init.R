
library("DESeq2")
args = commandArgs(trailingOnly = TRUE)

threads = args[1]
counts = args[2]
params = args[3]
output = args[4]
dds_design = args[5]
row_names = args[6]

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(counts, header=TRUE, row.names=1, sep="\t")
coldata <- read.table(params, header=TRUE, row.names=row_names,sep="\t")

coldata <- coldata[colnames(cts),]
cts <- cts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design= as.formula(paste('~',dds_design)))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=output)
