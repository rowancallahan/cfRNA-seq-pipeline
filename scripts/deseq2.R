library("DESeq2")


print('Setting parameters')

rds = snakemake@input[['rds']]
cat(sprintf(c('RDS object: ',rds,'\n')))

condition = snakemake@params[['condition']]
cat(sprintf(c('Condition: ',condition,'\n')))

ma_plot = snakemake@output[['ma_plot']]
cat(sprintf(c('MA plot', ma_plot,'\n')))

out_table = snakemake@output[['table']]
cat(sprintf(c('Summary results table', out_table,'\n')))

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(rds)

contrast <- c(condition, snakemake@params[["contrast"]])
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

p_hist = snakemake@output[['p_hist']]
pdf(p_hist)
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white", main='P values for genes with mean normalized count larger than 1',xlab='pvalue')
dev.off()

write.table(as.data.frame(res), file=out_table, quote=FALSE, sep='\t')
