library("DESeq2")

args = commandArgs(trailingOnly = TRUE)

print('Setting parameters')

threads = args[1]
print(c('threads: ',threads))

rds = args[2]
print(c('RDS object: ', rds))

contrast1 = args[3]
print(c('baseline: ', contrast1))

contrast2 = args[4]
print(c('contrast: ',contrast2))

condition = args[5]
print(c('Condition: ',condition))

ma_plot = args[6]
print(c('MA plot', ma_plot))

out_table = args[7]
print(c('Summary results table', out_table))

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
pdf(ma_plot)
plotMA(res, ylim=c(-2,2))
dev.off()


write.table(as.data.frame(res), file=out_table)
