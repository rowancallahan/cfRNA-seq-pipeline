library("DESeq2")

args = commandArgs(trailingOnly = TRUE)
rds = args[1]
output = args[2]
rds = args[3]

# load deseq2 data
dds <- readRDS(rds)

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)

svg(output)
plotPCA(counts, intgroup=snakemake@params[["pca_labels"]])
dev.off()

saveRDS(counts, file=rds)
