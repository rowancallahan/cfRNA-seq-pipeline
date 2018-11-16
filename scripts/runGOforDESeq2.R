args <- commandArgs()

help <- function(){
    cat("runGOforDESeq2.R :
- For the deseq2 output in the pipeline, run GO analysis on significant genes.
- Currently this is compatible with genome assemblies hg19 (Ens75), hg38 (Ens89) and hg38 (Ens90)
- Input DESEq2 table must be in .txt format
- Color options can be hex followed by saturation ex. #FE117A60 or rcolors
- The plot will have the same name as the degFile but with a .pdf extension.\n")
    cat("Usage: \n")
    cat("--degFile  : deseq2 table with log2FoldChange and pvalue [ required ]\n")
    cat("--adjp     : FDR adjusted p-value cutoff                 [ default = 0.01 ]\n")
    cat("--assembly : genome assembly                             [ requires ]\n")
    cat("--FC       : fold change cutoff (not log2 transformed)   [ default = 2 ]\n")
    cat("--printTree: option to print GOterm tree (0/1)           [ default = 0 ]\n")
    cat("\n")
    q()
}

## Save values of each argument
if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args))){
    help()
} else {
    degFile   <-sub('--degFile=', '', args[grep('--degFile=', args)])
    adjp      <-sub('--adjp=', '', args[grep('--adjp=', args)])
    assembly  <-sub('--assembly=', '', args[grep('--assembly=', args)])
    FC        <-sub('--FC=', '', args[grep('--FC=', args)])
    printTree <-sub('--printTree=', '', args[grep('--printTree=', args)])
}

if (identical(adjp,character(0))){
   adjp<-0.01
}else{
    adjp <- as.numeric(adjp)
}

if (identical(FC,character(0))){
    FC <- 2
} else{
    FC <- as.numeric(FC)
}

if (identical(printTree,character(0))){
    printTree <- 0
} else{
    printTree <- as.numeric(printTree)
}

library(GO.db)
library(topGO)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(GenomicFeatures)
library(Rgraphviz)

##----------load differentially expressed genes --------#
print("Loading differential expressed gene table")
print(degFile)

if(grepl('rda$|RData$',degFile)){
   deg <- get(load(file=degFile))
}
if(grepl('txt$|tsv$',degFile)){
    deg <- read.delim(file=degFile,header=TRUE,sep="\t")
}
rownames(deg) <- sub("\\.[0-9]*","",rownames(deg))

##---------load correct Biomart------------------------#
if (assembly == "hg19") {
    organismStr <- "hsapiens"
    geneID2GO <- get(load("/home/groups/CEDAR/anno/biomaRt/hg19.Ens_75.biomaRt.GO.geneID2GO.RData"))
    xx <- get(load("/home/groups/CEDAR/anno/biomaRt/GO.db.Term.list.rda"))
}
if (assembly == "hg38.89") {
    organismStr <- "hsapiens"
    ### to get to hg38 mappings ensembl 89!
    geneID2GO <- get(load("/home/groups/CEDAR/anno/biomaRt/hg38.Ens_89.biomaRt.GO.geneID2GO.RData"))
    xx <- get(load("/home/groups/CEDAR/anno/biomaRt/GO.db.Term.list.rda"))
}
if (assembly == "hg38.90") {
    organismStr <- "hsapiens"
    ### to get to hg38 mappings ensembl 90!
    geneID2GO <- get(load("/home/groups/CEDAR/anno/biomaRt/hg38.Ens_90.biomaRt.GO.geneID2GO.RData"))
    xx <- get(load("/home/groups/CEDAR/anno/biomaRt/GO.db.Term.list.rda"))
}

##-----------------------------------Functions--------------------------------------#
runGO <- function(geneList,xx=xx,otype,setName){
    setLength       <- sum(as.numeric(levels(geneList))[geneList]) 
    fname           <- paste(Dir, paste(setName, otype, "GO.txt", sep="_"), sep="/")
    GOData          <- new("topGOdata", ontology=otype, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher    <- runTest(GOData, algorithm = "classic", statistic = "fisher")## statistical test for topGO
    x               <- GenTable(GOData, classicFisher=resultFisher, topNodes=length(names(resultFisher@score)))## make go table for all terms
    x               <- data.frame(x)
    pVal            <- data.frame(pval=signif(resultFisher@score, 6)) ## get unrounded pvalue
    x$enrich        <- x$Significant/x$Expected ## calculate enrichment based on what you expect by chance
    x$p.unround     <- pVal[x$GO.ID,"pval"]## put unrounded pvalue in the table
    x$p.adj         <- signif(p.adjust(x$p.unround, method="BH"), 6)## calculate the adjusted pvalue with Benjamini & Hochberg correction
    x$log.p.adj     <- -log10(x$p.adj) ## convert adjusted p value to -log10 for plot magnitude
    #x$Term.full     <- sapply(x$GO.ID, FUN=function(n){Term(xx[[n]])}) ## get the full term name
    x <- x[order(x$GO.ID),]
    write.table(x, file=fname, sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE) ## save the table
    ## you can print the tree if you want, but since I keep the list of all of them skip
    if(printTree>0){
        printGraph(GOData,## make the tree for the go data
                   resultFisher,
                   firstSigNodes = 5,
                   fn.prefix = sub("_GO.txt$", "", fname),
                   useInfo = "all",
                   pdfSW = TRUE
                   )
    }
    return(x)  
}

## function to make barplot of -log10 adjusted pvalues colored by enrichment

drawBarplot <- function(go, ontology, setName){
    go <- go[!go$p.adj > 0.01,]
    if(nrow(go)>1){
        #go$Term.full <- make.unique(paste(sapply(strsplit(as.character(substring(go$Term.full,1,50)), "\\,"), `[`, 1)))
        go$Term <- make.unique(paste(sapply(strsplit(as.character(substring(go$Term,1,50)), "\\,"), `[`, 1)))
        print(setName)
        go <- go[with(go, order(p.adj, -enrich)),]
        ## Currently there is a discrepency between xx and x, so we only use Term right now, not Term.full
        #go$Term.full <-factor(paste(go$Term.full), levels=rev(paste(go$Term.full))) ## sort table by adjusted p-value
        go$Term <-factor(paste(go$Term), levels=rev(paste(go$Term))) ## sort table by adjusted p-value
        ptitle <- paste(ontology, setName) ## plot title
        ptitle <- gsub("^.*/","",ptitle)
        pfname <- paste(setName,ontology,"pdf",sep=".")## name of png file
        if(nrow(go) < 20 ){
            toprange <- 1:nrow(go)
        }else{
            toprange <- 1:20
        }
        top <- go[toprange,]
        col <- colorRampPalette(c("white","navy"))(16)   
        pdf(file=paste(Dir, pfname, sep="/"),height=5,width=7)
        print({
           p <- ggplot(top, aes(y=log.p.adj, x=Term, fill=enrich)) + ## ggplot barplot function
               geom_bar(stat="identity",colour="black") +
               ggtitle(ptitle) +
               xlab("") + ylab("-log10(fdr)") +
               scale_fill_gradient(low=col[2], high=col[15], name="enrichment", limits=c(0,ceiling(max(top$enrich))))+
               coord_flip()+
               theme(panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
               theme(text = element_text(size=8),
                     axis.text.x = element_text(vjust=1,color="black",size=8),
                     axis.text.y = element_text(color="black",size=8),
                     plot.title=element_text(size=10))   
       })
       dev.off()
    }
}

print("get up genes and make geneList")
up <- deg$padj < adjp & deg$log2FoldChange >= log2(FC)
up <- unique(rownames(deg[up,]))
all <-unique(names(geneID2GO))
up.geneList <-  factor(as.integer(all %in% up))
names(up.geneList) <- all

up.setsize <- sum(as.numeric(levels(up.geneList))[up.geneList])
print("setsize for significant genes") 
up.setsize

adjplabel <- gsub("^0\\.","",adjp)
comparison <- gsub("\\.tsv$|\\.txt$|\\.rda$|\\.RData$","",degFile)

Dir <- sub("$", "/GOterms", dirname(comparison))
if(!(file.exists(Dir))) {
      dir.create(Dir,FALSE,TRUE)
}

print("make GO table for the up genes")
#################################
go.UP.BP <- runGO(geneList=up.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"upFC",FC, "adjp", adjp, sep="."))
#go.UP.MF <- runGO(geneList=up.geneList,xx=xx,otype="MF",setName=paste(comparison,"up",sep="."))

print("make the png for the up genes")
drawBarplot(go=go.UP.BP,ontology="BP",setName=paste(basename(comparison),"upFC",FC, "adjp", adjp, sep="."))
    
print("get down genes and make geneList")
dn <- deg$padj < adjp & deg$log2FoldChange <= -log2(FC)
dn <- unique(rownames(deg[dn,]))
all <-unique(names(geneID2GO))
dn.geneList <-  factor(as.integer(all %in% dn))
names(dn.geneList) <- all

dn.setsize <- sum(as.numeric(levels(dn.geneList))[dn.geneList])
print("setsize for significant genes") 
dn.setsize

print("make GO table for down genes")
go.DN.BP <- runGO(geneList=dn.geneList,xx=xx,otype="BP",setName=paste(basename(comparison),"downFC",FC, "adjp", adjp, sep="."))

print("make barplot for down genes")
drawBarplot(go=go.DN.BP,ontology="BP",setName=paste(basename(comparison),"downFC",FC, "adjp", adjp, sep="."))
