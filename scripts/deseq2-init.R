
library("DESeq2")

counts = snakemake@input[['counts']]

params = snakemake@params[['samples']]

output = snakemake@output[['rds']]

dds_design = snakemake@params[['design']]

row_names = snakemake@params[['row_names']]

out_table = snakemake@output[['normed_counts']]

rld_out = snakemake@output[['rld_out']]

parallel <- FALSE
if (snakemake@threads > 1) {
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

normed_counts <-counts(dds,normalized=TRUE)
write.table(normed_counts,quote=F,sep='\t',file=out_table)

# obtain normalized counts
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file=rld_out)
