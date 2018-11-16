args <- commandArgs()

help <- function(){
    cat("RNAseq_makeVolcano.R :
- For the deseq2 output in the pipeline, make a volcano plot.
- Currently it calculates the adjusted p-val using the BH method, as there are many NAs in the table.
- Color options can be hex followed by saturation ex. #FE117A60 or rcolors
- The plot will have the same name as the degFile but with a .pdf extension.\n")
    cat("Usage: \n")
    cat("--degFile : deseq2 table with log2FoldChange and pvalue [ required ]\n")
    cat("--adjp    : FDR adjusted p-value cutoff                 [ default = 0.01 ]\n")
    cat("--FC      : fold change cutoff (not log2 transformed)   [ default = 2 ]\n")
    cat("--upCol   : up regulated genes color                    [ default = red ]\n")
    cat("--downCol : down regulated genes color                  [ default = blue ]\n")    
    cat("--ncCol   : non-changing genes color                    [ default = grey ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args))){
    help()
} else {
    degFile  <-sub('--degFile=', '', args[grep('--degFile=', args)])
    adjp     <-sub('--adjp=', '', args[grep('--adjp=', args)])
    FC       <-sub('--FC=', '', args[grep('--FC=', args)])
    upCol    <- sub('--upCol=', '', args[grep('--upCol=', args)])
    downCol  <- sub('--downCol=', '', args[grep('--downCol=', args)])
    ncCol    <- sub('--ncCol=', '', args[grep('--ncCol=', args)])
}

## set defaults if options are not provided
if (identical(adjp,character(0))){
   adjp<-0.01
}else{
    adjp <- as.numeric(adjp)
}

if (identical(FC,character(0))){
   FC <- 2
}else{
    FC <- as.numeric(FC)
}

if (identical(downCol,character(0))){
   downCol <- "blue"
}
if (identical(ncCol,character(0))){
   ncCol <- "grey"
}
if (identical(upCol,character(0))){
   upCol <- "red"
}

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)

## check if an rda file or tab sep
if(grepl('rda|RData',degFile)){
    deg <- get(load(file=degFile))
}
if(grepl('txt|tsv',degFile)){
    deg <- read.delim(file=degFile)
}

head(deg)
dim(deg)

## set all NA missing p-values to 1 (NA is DESeq2 default)
deg[is.na(deg$padj), "padj"] <- 1

## select up regulated genes
up <- deg$padj < adjp & deg$log2FoldChange > log2(FC)
sum(up)

## select down regulated genes
down <- deg$padj < adjp & deg$log2FoldChange < -log2(FC)
sum(down)

## set labels for pdf
adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.txt$|\\.rda|\\.tsv","",degFile)
pdfFile <- paste(comparison,adjplabel,"VolcanoPlot.pdf",sep=".")
print(pdfFile)

## remove directory to use in plot
comparison <- gsub("^.*/","",comparison)

Dir <- sub("$", "/Volcano", dirname(comparison))
if(!(file.exists(Dir))) {
      dir.create(Dir,FALSE,TRUE)
}

## calculate the -log10(adjp) for the plot
deg$log10padj <- -log10(deg$padj)

pdf(pdfFile,height=6,width=6)
plot(deg$log2FoldChange
    ,deg$log10padj
    ,col=ncCol
    ,pch=19
    ,cex=1
    ,main=comparison,xlab="log2(Fold Change)"
    ,ylab="-log10(padj)"
    ,cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
points(deg[up, "log2FoldChange"],deg[up,"log10padj"], col=upCol, pch=19,cex=1)
points(deg[down, "log2FoldChange"],deg[down, "log10padj"], col=downCol, pch=19,cex=1)
legend("topright", c(paste(sum(up),"genes p.adj-val < ",adjp),paste(sum(down),"genes p.adj-val < ",adjp)),
       pch=c(19,19,1), col=c(upCol,downCol,ncCol),cex=0.9,)
abline(v = -log2(FC), lty="dashed", col="black")
abline(v = log2(FC), lty="dashed", col="black")
abline(h = -log10(adjp), lty="dashed", col="black")
dev.off()
