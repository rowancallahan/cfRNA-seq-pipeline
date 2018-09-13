library("DESeq2")
library("ggplot2")
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("genefilter")

cat(sprintf(c('Working directory',getwd())))

cat(sprintf('Setting parameters'))

rds <- snakemake@input[['rds_object']]
cat(sprintf(c('RDS: ',rds)))

dds <- snakemake@input[['dds_object']]
cat(sprintf(c('DDS: ',rds)))

pca_plot <- snakemake@output[['pca']]
cat(sprintf(c('PCA plot: ',pca_plot)))

labels <- snakemake@params[['pca_labels']]
cat(sprintf(c('PCA Labels: ',labels)))

sd_mean_plot <- snakemake@output[['sd_mean_plot']]
cat(sprintf(c('SD Mean plot: ',sd_mean_plot,'\n')))

heatmap_plot <- snakemake@output[['heatmap_plot']]
cat(sprintf(c('Heatmap plot: ',heatmap_plot,'\n')))

distance_plot <- snakemake@output[['distance_plot']]
cat(sprintf(c('Distance plot: ',distance_plot,'\n')))

panel_ma <- snakemake@output[['panel_ma']]
cat(sprintf(c('panel MA plot: ',panel_ma,'\n')))

var_heat <- snakemake@output[['var_heat']]
cat(sprintf(c('Top Variance Genes Heatmap: ',var_heat,'\n')))

ggplot_pca_factor <- snakemake@output[['ggplot_pca_factor']]
cat(sprintf(c('ggplot2 PCA factor plot: ', ggplot_pca_factor,'\n')))

cat(sprintf('Load rlog DESeqTransform object'))
rld <-readRDS(rds)

cat(sprintf('Load dds DESeqTransform object'))
dds <-readRDS(dds)

pdf(pca_plot)
plotPCA(rld, intgroup=labels)
dev.off()

# SD mean plot
pdf(sd_mean_plot)
meanSdPlot(assay(rld))
dev.off()

# Heatmap top 20 DE genes
# TODO Add flexible/extensible filtering criteria to parameter configuration file for plotting utilities
res <-results(dds)
topGenes <- head(order(res$padj),20)

df <- as.data.frame(colData(rld))


pdf(heatmap_plot)
pheatmap(assay(rld)[topGenes,], cluster_rows=T, fontsize=6,fontsize_row=6,fontsize_col=6,show_rownames=T,
         cluster_cols=T, annotation_col=df[,c('Type','RNA_extracted_by','Lib_prep_date')],labels_col=as.character(df$SampleID))
dev.off()

# Heatmap of distances

pdf(distance_plot)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

# TODO allow for dynamic row naming conventions
rownames(sampleDistMatrix) <- paste(rld$Type, rld$RNA_extracted_by, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, fontsize=5,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# ggplot PCA

pcaData <- plotPCA(rld, intgroup=c("Type", "RNA_extracted_by"), returnData=TRUE)

pdf(ggplot_pca_factor)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# TODO allow for dynamic color and shape naming conventions
g<-ggplot(pcaData, aes(PC1, PC2, color=Type, shape=RNA_extracted_by)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
g + theme_gray()
dev.off()

#panel ma plot
pdf(panel_ma)
par(mfrow=c(2,2),mar=c(2,2,1,1) +0.1)
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
mtext(resG@elementMetadata$description[[2]], outer=T, cex=.6,line=-1)
dev.off()


pdf(var_heat)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld))
pheatmap(mat, annotation_col = anno[,c('Type','RNA_extracted_by','Lib_prep_date')],fontsize=6,)
dev.off()