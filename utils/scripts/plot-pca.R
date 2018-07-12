library("DESeq2")
library("ggplot2")
library("pheatmap")
library("vsn")
library("RColorBrewer")

rds = snakemake@input[['rds_object']]
pca_plot = snakemake@output[['pca']]
rld_out = snakemake@output[['rld_out']]
labels = snakemake@params[['pca_labels']]
sd_mean_plot =snakemake@output[['sd_mean_plot']]
heatmap_plot =snakemake@output[['heatmap_plot']]
distance_plot =snakemake@output[['distance_plot']]
ggplot_pca =snakemake@output[['ggplot_pca']]
panel_ma =snakemake@output[['panel_ma']]

# load deseq2 data
dds <- readRDS(rds)

# obtain normalized counts
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file=rld_out)


svg(pca_plot)
plotPCA(rld, intgroup=labels)
dev.off()

# SD mean plot
svg(sd_mean_plot)
meanSdPlot(assay(rld))
dev.off()

# Heatmap top 20 DE genes

res <-results(dds)
topGenes <- head(order(res$padj),20)

df <- as.data.frame(colData(dds))


svg(heatmap_plot)
pheatmap(assay(rld)[topGenes,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

# Heatmap of distances

svg(distance_plot)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Type, rld$RNA_extracted_by, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# ggplot PCA
svg(ggplot_pca)
pcaData <- plotPCA(rld, intgroup=c("Type", "RNA_extracted_by"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Type, shape=RNA_extracted_by)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()

#panel ma plot
svg(panel_ma)
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
dev.off()

