library("DESeq2")

args = commandArgs(trailingOnly = TRUE)

threads = args[1]
rds = args[2]
contrast1 = args[3]
contrast2 = args[4]
condition = args[5]
ma_plot = args[6]
out_table = args[7]

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(threads))
    parallel <- TRUE
}

dds <- readRDS(rds)

contrast <- c(condition, contrast1, contrast2)
res <- results(dds, contrast=contrast, parallel=parallel)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# MA plot
svg(ma_plot)
plotMA(res, ylim=c(-2,2))
dev.off()



write.table(as.data.frame(res), file=out_table)
