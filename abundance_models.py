
qc_model = """
require(data.table)
require(NMF)
require(affy)
require(limma)
require(biomaRt)
library(stringr)
source("{code_dir}AbundanceFunctions/BiasReduce.R")
source("{code_dir}AbundanceFunctions/ExtractTransformLoad.R")
source("{code_dir}AbundanceFunctions/DifferentialAnalysis.R")
source("{code_dir}AbundanceFunctions/NonVisualOutput.R")
source("{code_dir}GenomicsFunctions/ReadAndParse.R")
source("{code_dir}AssociationFunctions/gs.wrapper.R")
source("{code_dir}AssociationFunctions/PathwayAnalysis.R")
source("{code_dir}BcorePlotting/SummaryPlots.R")
source("{code_dir}BcorePlotting/MultipleTestingCorrection.R")
source("{code_dir}BcorePlotting/ClusteringPlots.R")


setwd("{results}")

dir.create(file.path(getwd(),'{project_title}'),showWarnings=FALSE)
oneclass = {oneclass}
# constants
# max data rows for hclust
hclust.limit = 2^16
annCol.lm_by = {annCollm_by}
# quantile of data distribution reqd in one group's worth of data if too many rows for hclust()
hc.frac.cut = 0.75;
SJ.counts.na.frac = 0.25;
# max fraction of samples not having detected a splice junction for the splice
# junction to be retained in raw data


# regression parameters
na.lim = 0 # max NAs per row tolerated by lm() at least in some cases
do.not.regress = "alograw" # control norm not to be used for regression stats
# plotting colors
colors.rgb = c(rgb(0,0,0),rgb(0.1,0.1,1),rgb(0,.7,.7),rgb(0,.7,0),rgb(.7,1,0),rgb(.7,0,.7))
md.file = "{meta_file}"
md.orientation = "byRow" # sample_ids are in @ row. alt:byCol (IDs in @ col)
md.IDcol = "{sample_id}" # reqd if md.orientation is byRow; byCol==headers are IDs

# gene annotation
taxID = {tax_id}
gene2ENSfile = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/ncbi/gene2ensembl.gz"
gene2ENS.col = c("taxID","EntrezID","Gene","RefSeqTranscript","EnsemblTranscript","RefSeqProtein","EnsemblProtein")
gtf_file = "{gtf_file}"
gtf.feature = "{gtf_feature}"
gtf.orig.col = c("gene_id","gene_name","gene_biotype")
gtf.col = c("Gene","Symbol","biotype")
proj.title = "{project_title}"
gtf.Rdir = "{gtf_read_dir}"
genome.func = "{code_dir}GenomicsFunctions"


readir = "{read_dir}"
readpattern = "{read_pattern}"
useme.cols = "{useme_cols}"
label.from.colname = "{label_from_colname}"
samps = dir(path=readir, pattern=readpattern)
samp.labels = gsub(label.from.colname,'\\\\1', samps)

annCol.names = "group"
annCol.label = "{label_from_colname}"

# read in STAR alignments

print("Read in Abundance data")

STAR.data = read_STAR(readir=readir,useme.cols=useme.cols,label.from.colname=label.from.colname,annCol.label=annCol.label,annCol.names=annCol.names,annCol.normBy=NULL,annCol.lm_by=annCol.lm_by, readpattern=readpattern,outypes=list(gene.counts='Gene.out.tab'), unstranded.col=list(gene.counts=c(1:4)),outkey=list(gene.counts=1),outcol=list(gene.counts=c("Gene","Reads","FwdReads","RevReads")),outsum=list(gene.counts = c("N_ambiguous","N_multimapping","N_noFeature","N_unmapped")))

print("Read in Abundance data : Complete")

print("Filter SJ.counts data")
# filter out SJs with too many NAs
{{
    if(any( names(STAR.data$LoM.raw)=="SJ.counts" & exists("SJ.counts.na.frac")))
    {{
      STAR.data$SJ.counts.orig = STAR.data$LoM.raw$SJ.counts
      for(tag in names(STAR.data$LoM.raw$SJ.counts))
      {{
        STAR.data$LoM.raw$SJ.counts[[tag]] = STAR.data$LoM.raw$SJ.counts[[tag]][rowSums(is.na(STAR.data$LoM.raw$SJ.counts[[tag]]))<= (SJ.counts.na.frac*ncol(STAR.data$LoM.raw$SJ.counts[[tag]])), ]
      }}
    }}
}}

print("Filtering SJ.counts : Complete")
save.image("./{project_title}.RData")

gtf = readENSgtf(filename=gtf_file)
genes.gtf = gtf[feature==gtf.feature, mget(gtf.orig.col)]
names(genes.gtf) = gtf.col
setkeyv(genes.gtf,gtf.col[1])

##### parse of NCBI's EntrezID to Ensembl translation
Entrez2Ensembl = fread(paste("zgrep",paste0("-E '^",taxID,"'"),gene2ENSfile))
names(Entrez2Ensembl) = gene2ENS.col
setkeyv(Entrez2Ensembl,gene2ENS.col[3])

##### add EntrezIDs to genes.gtf
tmp = Entrez2Ensembl[,mget(gene2ENS.col[2:3])]; tmp=tmp[!duplicated(tmp),]
genes.gtf = merge(genes.gtf,tmp,all.x=TRUE)
rm(tmp)

myreads = STAR.data$myreads

print("Filter gene.counts data")

{{
    if(any( names(STAR.data$LoM.raw)=="gene.counts"))
    {{
      STAR.data$LoM.orig.raw = STAR.data$LoM.raw$gene.counts
      for(tag in names(STAR.data$LoM.raw$gene.counts))
      {{
        STAR.data$LoM.raw$gene.counts[[tag]] = STAR.data$LoM.raw$gene.counts[[tag]][rowSums(STAR.data$LoM.raw$gene.counts[[tag]])>0, ]
      }}
    }}
}}

print("Filtering gene.counts : Complete")

save.image("./{project_title}.RData")

LoM.norms = vector(mode='list',length=length(STAR.data$LoM.raw))
names(LoM.norms) = names(STAR.data$LoM.raw)


print("Bias reduction with normMatrix")

for(tag in names(STAR.data$LoM.raw))
{{
  if( length(STAR.data$LoM.raw[[tag]])>1 )
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[myreads]])
  }}
  else
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[1]])
  }}
}}

print("Bia reduction with normMatrix : Complete")

save.image("./{project_title}.RData")
## Iterate through bias reduced mats and write to working directory

#for (n in names(LoM.norms[['gene.counts']])){{
#    out_table = paste(n ,'.txt',sep='')
#    write.table(LoM.norms[['gene.counts']][[n]],paste(getwd(),proj.title,out_table,sep='/'),sep='\t',quote=F)
#}}


## custom: reannotate from metadata
## run if necessary metadata is in file and not also in FASTQ file names
# map annotation to read matrix
md.dt = fread("{meta_file}")

print("Evaluate meta-data file for sample name consistency")

{{
    if( md.orientation == "byCol" ){{ # samples are one per column
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), names(md.dt) )
    }}
    else
    {{ # samples are one per row
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), md.dt[,get(md.IDcol)] )
    }}
}}

idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
{{
    if( sum(idx1==1:ncol(STAR.data$LoM.raw[[1]][[myreads]]))==ncol(STAR.data$LoM.raw[[1]][[myreads]]))
    {{
      if( md.orientation == "byCol" )
      {{ # samples are one per column
        annCol = NULL
        namecol = setdiff( 1:ncol(md.dt), idx2 )
        md.factors = as.character(md.dt[,namecol,with=F])
        for(k in 1:nrow(md.dt) )
        {{
          annCol[[ md.dt[k,get(names(md.dt)[namecol])] ]][idx1] = as.vector( md.dt[ k, mget(names(md.dt)[idx2]) ] )
        }}
      }}
      else
      {{ # samples are one per row
        annCol = NULL
        namecol = setdiff(names(md.dt), md.IDcol)
        md.factors = as.character(namecol)
        for( k in setdiff(names(md.dt),md.IDcol))
        {{
          annCol[[ k ]][idx1] = as.vector( md.dt[ idx2, get(k) ] )
        }}
      }}
    }}
    else
    {{
      stop(paste("Some samples have no annotation in",md.file))
    }}
}}

print("Evaluation of meta-data : Complete")

mytypes = names(LoM.norms)
mynorms = names(LoM.norms[[1]])

#grpBy = annCol[[annCol.lm_by]]
grpBy = annCol[[oneclass]]
annCol.plotme = {ann_colplotme}
clim.pct=0.96
histbins=20
dir.create(file.path(getwd(),'{project_title}/summary.plotsPlots'),showWarnings=FALSE)
dir.create(file.path(getwd(),'{project_title}/qc.clustersPlots'),showWarnings=FALSE)

print("Generate summary and qc.cluster plots")

#### Loop through mynorms ######

for(i in 1:length(mytypes))
{{
  for(j in 1:length(mynorms))
  {{
    # set up matrices and config for plotting
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}

    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels
    plotdata = list(plotdir='./{project_title}/summary.plotsPlots',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)))

    if(sum(rowmask)>hclust.limit)
    {{
      rowmask = rowSums(rawmat>=quantile(rawmat,probs=hc.frac.cut,na.rm=T),na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)) )
    }}

    ans = summary.plots(rawmat=log2(rawmat +1), normmat=normmat, mynorm=mynorms[j], samp.labels=samp.labels, samp.classes=grpBy, plotdata=plotdata,plot2file=TRUE,histbins=histbins, colorspec=colors.rgb)
    plotdata = list(plotdir='./{project_title}/qc.clustersPlots',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    ans = qc.clusters(rawmat=log2(rawmat[rowmask,] +1), normmat=normmat[rowmask,], attribs=annCol[annCol.plotme], oneclass=oneclass, colorspec=colors.rgb, plotdata=plotdata, plot2file=TRUE, clim.pct=clim.pct)
  }}
}}

print("Generate summary and qc.cluster plots : Complete")

save.image("./{project_title}.RData")

regress_lsls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(regress_lsls) = names(STAR.data$LoM.raw)

contr_ls = {contr_ls}

# set baseline for regression in parameter(s) of interest
# contr.treatment generates regression coefficients that are like (adjusted mean) ratios of other groups to the baseline group.
# contr.sum generates coefficients that are like (adjusted mean) ratios to average all for all but the mandadory ommitted treatment group (because there is always one fewer independent pairwise comparison than there are pairs).


lm_expr = "{lm_expr}"

rowmask_ls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(rowmask_ls) = names(STAR.data$LoM.raw)
dir.create(file.path(getwd(),'{project_title}/regressMatrixPlots'),showWarnings=FALSE)

print("RegressMatrix with specified contrasts")

for(i in 1:length(mytypes))
{{
  for(j in which(!mynorms %in% do.not.regress))
  {{
    # set up data matrices
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}
    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp=paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./{project_title}/regressMatrixPlots/',plotbase=paste(mynorms[j],mytypes[i],'vs',tmp,sep='.'),plottitle=proj.title)
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels

    # prepare rowmask for heatmap/MDS (remove non-expr or low expr>hclustlim)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) <= na.lim
    rowmask_ls[[ mytypes[i] ]][[ mynorms[j] ]] = rowmask # save for later

    # regression

    regress_lsls[[mytypes[i]]][[mynorms[j]]] = regressMatrix(normmat[rowmask,], expt.design=annCol[annCol.lm_by], lm_expression=lm_expr, contr_list = contr_ls, plot2file = TRUE, plotdata = plotdata)
  }}
}}

print("RegressMatrix with specified contrasts : Complete")

save.image("./{project_title}.RData")

topn=500 # number with which to play

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress) )
  {{
    ans = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    cat(mytypes[i],mynorms[j],"1-p0:",signif(unlist(lapply(2:length(ans),function(x){{1-ans[[x]]$pi0}})),2),"\n")
    cat(mytypes[i],mynorms[j],"qcut:",signif(unlist(lapply(2:length(ans),function(x){{quantile(ans[[x]]$qvalues,probs=topn/length(ans[[x]]$qvalues))}})),3),
    "returns ~",topn,"out of",length(ans[[2]]$qvalues),"\n")
  }}
}}

# Plots for each design factor (default) or for factors specified in facSel
#   1) Histogram of p values that were included in the design
#      Look for a peak on the left, and no peaks in the middle or on the right
#   2) qvalue's default plots, with full qvalue range c(0,1) plotted
#      Look for descent to pi0 in the top left tuning plot with good asymptote
#      The slow/steep rise in q-values in remaining plots depends on resolving power of data


lm_by = annCol.lm_by
histbins=20
select_lsmk = vector(mode='list',length=length(STAR.data$LoM.raw))
names(select_lsmk) = names(STAR.data$LoM.raw)
dir.create(file.path(getwd(),'{project_title}/qcQvaluePlots'),showWarnings=FALSE)

print("Perform qcQvalue & Generate masks")

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(normmat)=samp.labels
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]
    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    # pull out regression design and set up plot config

    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp=paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./{project_title}/qcQvaluePlots',plotbase=paste(mynorms[j],mytypes[i],tmp,sep='.'),plottitle=proj.title)

    ##for(fac in names(reg_ls$q_list)[grep(lm_by,names(reg_ls$q_list))])

    for(fac in names(reg_ls$q_list))
    {{
      if(!fac %in% '(Intercept)')
      {{
        pdata = plotdata
        pdata$plotbase = paste(plotdata$plotbase,make.names(fac),sep='.')
        ans = qcQvalues(norm_x=normmat[mymask,], pvalue_v=reg_ls$p_mat[,fac], obj_qvalue=reg_ls$q_list[[fac]], attribs=annCol[lm_by],
        oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, histbins=histbins, plot2file=TRUE)
        select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac)] = ans['rowmask']
      }}
    }}
  }}
}}

print("Perform qcQvalue & Generate masks : Complete")

save.image("./{project_title}.RData")

# Build of ratios to baselines for desired design factors
# Also build masks selecting genes based on q-values per factor
# and follow-up considerations such as expression level and fold change
# Loop over q-value cuts to assist with final cut selection.
# Plot heatmaps and MDS plots based on selected genes and designed ratios


ngene_v = c(200,500,1000) # q-value cuts by number; can also cut by q-value
ratioby_ls = list("{lm_by}"=contr_ls${lm_by}$baseline)
ratio_fold = 1.3
intensity_fold = 2

cut_ls = list(q_combine="OR", rcut_fold=ratio_fold, icut_fold=intensity_fold)

# settings for heatmaps
#annCol.plotme = annCol.lm_by # heatmap tracks

clustrowmin = 10 # min data rows for heatmap and MDS plots

# save ratios and selections

ratio_lsmat = vector(mode='list',length=length(STAR.data$LoM.raw))
names(ratio_lsmat) = names(STAR.data$LoM.raw)

print("Perform PlotRatios")

dir.create(file.path(getwd(),'{project_title}/plotRatiosPlots'),showWarnings=FALSE)
for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(normmat)=samp.labels

    if( grepl('SJ',mytypes[i]))
    {{#not the most robust way to find SJs
      # map SJ positions to overlapping genes
      # specifying defaults to fix "Gene" and Pos as colnames
      gtf.col = c("Gene","Symbol","Chr","start","stop")
      pos.col = "Pos"
      mySJ_dt = mapSJ2feature(STAR.data$LoM.raw[[ mytypes[i] ]][[1]], pos.col=pos.col, gtf.col=gtf.col, gtf.file=gtf_file, gtf.Rdir=genome.func)
      # min data rows for heatmap and MDS plots
      # map gene IDs back to data matrix using mapping built earlier
      # this enables include_ID to select rows

      idx2 = match( rownames(STAR.data$LoM.raw[[ mytypes[i] ]][[1]]), mySJ_dt[,get(pos.col)] )
      idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
      rownames(normmat)[idx1] = mySJ_dt[idx2,get(gtf.col[1])]
    }}

    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    reg_ls = reg_ls[!grepl('Intercept',names(reg_ls))]
    plotdata = list(plotdir='./{project_title}/plotRatiosPlots',plotbase=paste(mynorms[j],sub('counts','ratios',mytypes[i]),'minratio',paste0(ratio_fold,'x'),'minexpr',round(min(normmat,na.rm=T)+log2(intensity_fold),1),sep='.'),plottitle=proj.title)

    #for(fac in names(reg_ls)[grep(lm_by,names(reg_ls))])
    for(fac in names(reg_ls))
    {{
        for( ngenes in ngene_v )
        {{
          # calculate and plot ratios
          cut_ls$qcut = ngenes
          pdata = plotdata
          pdata$plotbase = paste(plotdata$plotbase,'ratio.vs',ratioby_ls${lm_by},ngenes,sep='.')

          # this function makes the ratios and cuts
          ans = designRatios(normmat[mymask,], q_list=reg_ls[[fac]], attribs=annCol, ratioby_ls=ratioby_ls, cut_ls=cut_ls)

          # this function creates the heatmap and MDS plot
          if( sum(ans$rowmask)>clustrowmin )
          {{
            ans2= plotRatios( ratiomat=ans$ratiomat, attribs=annCol[annCol.plotme], oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, rowmask=ans$rowmask, plot2file=TRUE)
          }}
          # save selections
          select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac,ngenes,sep="_")] = ans['rowmask']
        }}
        # save ratiomat (not dependent on ngene cut) after ngene loop
        ratio_lsmat[[mytypes[i]]][[mynorms[j]]] = ans$ratiomat
    }}
  }}
}}

print("Perform PlotRatios : Complete")

save.image("./{project_title}.RData")

save.image("./{project_title}.RData")

print("Perform biomaRt gene annotation")

ratio_mat = ratio_lsmat[["{path_type}"]][["{path_norms}"]]
reg_ls = regress_lsls[["{path_type}"]][["{path_norms}"]]
ID = rownames(ratio_mat)


project_species = useMart("ENSEMBL_MART_ENSEMBL", dataset = "{mart_dataset}", host = "dec2016.archive.ensembl.org")
human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "dec2016.archive.ensembl.org")

print("biomaRt gene annotations : Loaded")

hsa_entrezID = getLDS(attributes = "ensembl_gene_id",mart = project_species,  attributesL = "entrezgene", martL = human)
Entrez.ID = character(length(ID))
Entrez.ID[] = NA
idx2 = match(ID,hsa_entrezID$Gene.ID)
idx1 = which(!is.na(idx2)); idx2= idx2[idx1]
Entrez.ID[idx1] = hsa_entrezID$EntrezGene.ID[idx2];
Entrez.ID = as.numeric(Entrez.ID)

###Human gene symbols

hsa_symbol = getLDS(attributes = "ensembl_gene_id",mart = project_species,  attributesL = "hgnc_symbol", martL = human)
Anno.Symbol = character(length(ID))
Anno.Symbol = NA
idx2 = match(ID,hsa_symbol$Gene.ID)
idx1 = which(!is.na(idx2)); idx2= idx2[idx1]
Anno.Symbol[idx1] = unlist(hsa_symbol$HGNC.symbol[idx2]);
# assemble IDs for annotation

print("Species event ID's mapped to Hsapiens_gene_ensembl")


backgroundset = as.data.table(cbind(ID, Entrez.ID, Anno.Symbol))

# assemble signatures

reg_ls = regress_lsls$gene.counts${path_norms}

masks = names(select_lsmk[["{path_type}"]][["{path_norms}"]])
dir.create(file.path(getwd(),'{project_title}/tables'),showWarnings=FALSE)

print("Perform pathway enrichment for all factor levels relative to baseline")

for (k in masks)
{{
    sig_gmt = NULL
    rowmask = select_lsmk[["{path_type}"]][["{path_norms}"]][[k]]
    sig_gmt[['all']] =  ID[rowmask]
    fileSettings = list(directory=file.path(getwd(),'{project_title}/tables'),baseFilename=paste(k,'with_Uni',sep='_'))
    path_ls = ags.wrapper(setlist=sig_gmt, backgroundset, include.identifiers=TRUE, anno.uni=TRUE,
        fileSettings = fileSettings,
        functiondir = "{code_dir}AssociationFunctions/",
        resourcedir = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/old_association_resources",
        return.OM=TRUE, ecut=0.05, ocut=5)

    fileSettings = list(directory=file.path(getwd(),'{project_title}/tables'),baseFilename=paste(k,'with_out_Uni',sep='_'))
    path_ls = ags.wrapper(setlist=sig_gmt, backgroundset, include.identifiers=TRUE, anno.uni=FALSE,
        fileSettings = fileSettings,
        functiondir = "{code_dir}AssociationFunctions/",
        resourcedir = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/old_association_resources",
        return.OM=TRUE, ecut=0.05, ocut=5)
}}

print("Perform pathway enrichment for all factor levels relative to baseline : Complete")

print("Generate Abundance and Ratio table with associated q-value and p-values")

gtf.Rdir = "{gtf_read_dir}"

out_norm_mat = LoM.norms$gene.counts$loess[rowmask_ls${path_type}${path_norms},]
out_table = outputTable(normmat= out_norm_mat, gtf.file = "{gtf_file}",ratiomat = ratio_lsmat${path_type}${path_norms}, q_list=reg_ls$q_list, gtf.Rdir=genome.func, gtf.key='transcript')
out_table = out_table[!duplicated(out_table$Gene),]

write.table(out_table, row.names = FALSE, file=file.path(getwd(), './{project_title}/tables', paste("{project_title}","{path_norms}","Normed_with_Ratio_and_Abundance.txt", sep="_")),quote=FALSE,sep='\t')

print("Generate Abundance and Ratio table with associated q-value and p-values : Complete")

save.image("./{project_title}.RData")
"""


qc_matrix_model = """
require(data.table)
require(NMF)
require(RColorBrewer)
require(affy)
require(limma)
require(biomaRt)
library(stringr)
source("{code_dir}AbundanceFunctions/BiasReduce.R")
source("{code_dir}AbundanceFunctions/ExtractTransformLoad.R")
source("{code_dir}AbundanceFunctions/DifferentialAnalysis.R")
source("{code_dir}AbundanceFunctions/NonVisualOutput.R")
source("{code_dir}GenomicsFunctions/ReadAndParse.R")
source("{code_dir}AssociationFunctions/gs.wrapper.R")
source("{code_dir}AssociationFunctions/PathwayAnalysis.R")
source("{code_dir}BcorePlotting/SummaryPlots.R")
source("{code_dir}BcorePlotting/MultipleTestingCorrection.R")
source("{code_dir}BcorePlotting/ClusteringPlots.R")


setwd("{results}")

oneclass = {oneclass}
# constants
# max data rows for hclust
hclust.limit = 2^16
annCol.lm_by = {annCollm_by}

# quantile of data distribution reqd in one group's worth of data if too many rows for hclust()
hc.frac.cut = 0.75;
SJ.counts.na.frac = 0.25;
# max fraction of samples not having detected a splice junction for the splice
# junction to be retained in raw data


# regression parameters
na.lim = 0 # max NAs per row tolerated by lm() at least in some cases
do.not.regress = "alograw" # control norm not to be used for regression stats
# plotting colors
colors.rgb = c(rgb(0,0,0),rgb(0.1,0.1,1),rgb(0,.7,.7),rgb(0,.7,0),rgb(.7,1,0),rgb(.7,0,.7))
md.file = "{meta_file}"
md.orientation = "byRow" # sample_ids are in @ row. alt:byCol (IDs in @ col)
md.IDcol = "{sample_id}" # reqd if md.orientation is byRow; byCol==headers are IDs

# gene annotation
taxID = {tax_id}
gene2ENSfile = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/ncbi/gene2ensembl.gz"
gene2ENS.col = c("taxID","EntrezID","Gene","RefSeqTranscript","EnsemblTranscript","RefSeqProtein","EnsemblProtein")
gtf_file = "{gtf_file}"
gtf.feature = "{gtf_feature}"
gtf.orig.col = c("gene_id","gene_name","gene_biotype")
gtf.col = c("Gene","Symbol","biotype")
proj.title = "{project_title}"
gtf.Rdir = "{gtf_read_dir}"
genome.func = "{code_dir}GenomicsFunctions"


readir = "{read_dir}"
readpattern = "{read_pattern}"
useme.cols = "{useme_cols}"
label.from.colname = "{label_from_colname}"
samps = dir(path=readir, pattern=readpattern)
samp.labels = gsub(label.from.colname,'\\\\1', samps)

annCol.names = "group"
annCol.label = "{label_from_colname}"

# read in STAR alignments

print("Read in Abundance data")


filt_pattern = str_replace_all(readpattern,'[.*]','')
tep_mat = as.matrix(read.table("../{dataset}",check.names=F))
filt_dt = data.table(tep_mat,check.names=F)
samps = sort(grep(filt_pattern, names(filt_dt),value=T))
samp.labels = samps
filt_mat = as.matrix(filt_dt[,samps,with=F])
row.names(filt_mat) = row.names(tep_mat)
STAR.data = list(LoM.raw=list(gene.counts=list(Reads=filt_mat)),expt.design=list(group=samps), myreads = 'Reads')
attr(STAR.data$expt.design,'names')='group'
attr(STAR.data$expt.design,'lm_by')=annCol.lm_by

myreads = STAR.data$myreads

print("Filter gene.counts data")

{{
    if(any( names(STAR.data$LoM.raw)=="gene.counts"))
    {{
      STAR.data$LoM.orig.raw = STAR.data$LoM.raw$gene.counts
      for(tag in names(STAR.data$LoM.raw$gene.counts))
      {{
        STAR.data$LoM.raw$gene.counts[[tag]] = STAR.data$LoM.raw$gene.counts[[tag]][rowSums(STAR.data$LoM.raw$gene.counts[[tag]])>0, ]
      }}
    }}
}}

print("Filtering gene.counts : Complete")

save.image("./{project_title}.RData")

LoM.norms = vector(mode='list',length=length(STAR.data$LoM.raw))
names(LoM.norms) = names(STAR.data$LoM.raw)


print("Bias reduction with normMatrix")
dir.create(file.path(getwd(),'tables'),showWarnings=FALSE)

for(tag in names(STAR.data$LoM.raw))
{{
  if( length(STAR.data$LoM.raw[[tag]])>1 )
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[myreads]],fileData_ls=list(fileDir='./tables/',fileBase='{project_title}'))
  }}
  else
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[1]],fileData_ls=list(fileDir='./tables/',fileBase='{project_title}'))
  }}
}}

print("Bia reduction with normMatrix : Complete")

save.image("./{project_title}.RData")
## Iterate through bias reduced mats and write to working directory


## custom: reannotate from metadata
## run if necessary metadata is in file and not also in FASTQ file names
# map annotation to read matrix
md.dt = fread("{meta_file}")

print("Evaluate meta-data file for sample name consistency")

{{
    if( md.orientation == "byCol" ){{ # samples are one per column
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), names(md.dt) )
    }}
    else
    {{ # samples are one per row
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), md.dt[,get(md.IDcol)] )
    }}
}}

idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
{{
    if( sum(idx1==1:ncol(STAR.data$LoM.raw[[1]][[myreads]]))==ncol(STAR.data$LoM.raw[[1]][[myreads]]))
    {{
      if( md.orientation == "byCol" )
      {{ # samples are one per column
        annCol = NULL
        namecol = setdiff( 1:ncol(md.dt), idx2 )
        md.factors = as.character(md.dt[,namecol,with=F])
        for(k in 1:nrow(md.dt) )
        {{
          annCol[[ md.dt[k,get(names(md.dt)[namecol])] ]][idx1] = as.vector( md.dt[ k, mget(names(md.dt)[idx2]) ] )
        }}
      }}
      else
      {{ # samples are one per row
        annCol = NULL
        namecol = setdiff(names(md.dt), md.IDcol)
        md.factors = as.character(namecol)
        for( k in setdiff(names(md.dt),md.IDcol))
        {{
          annCol[[ k ]][idx1] = as.vector( md.dt[ idx2, get(k) ] )
        }}
      }}
    }}
    else
    {{
      stop(paste("Some samples have no annotation in",md.file))
    }}
}}

print("Evaluation of meta-data : Complete")

mytypes = names(LoM.norms)
mynorms = names(LoM.norms[[1]])

grpBy = annCol[[oneclass]]
annCol.plotme = {ann_colplotme}
clim.pct=0.96
histbins=40
dir.create(file.path(getwd(),'1'),showWarnings=FALSE)
dir.create(file.path(getwd(),'2'),showWarnings=FALSE)

print("Generate summary and qc.cluster plots")

#### Loop through mynorms ######

for(i in 1:length(mytypes))
{{
  for(j in 1:length(mynorms))
  {{
    # set up matrices and config for plotting
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}

    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    normmat[is.na(normmat)] <-0
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels
    plotdata = list(plotdir='./1',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) )
    # Too stringent of a rowmask for CFRNA project
    # rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy))) 

    if(sum(rowmask)>hclust.limit)
    {{
      rowmask = rowSums(rawmat>=quantile(rawmat,probs=hc.frac.cut,na.rm=T),na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)) )
    }}

    ans = summary.plots(rawmat=log2(rawmat +1), normmat=normmat, mynorm=mynorms[j], samp.labels=samp.labels, samp.classes=grpBy, plotdata=plotdata,plot2file=TRUE,histbins=histbins, colorspec=colors.rgb)
    plotdata = list(plotdir='./2',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    ans = qc.clusters(rawmat=log2(rawmat[rowmask,] +1), normmat=normmat[rowmask,], attribs=annCol[annCol.plotme], oneclass=oneclass, colorspec=colors.rgb, plotdata=plotdata, plot2file=TRUE, clim.pct=clim.pct)
  }}
}}

print("Generate summary and qc.cluster plots : Complete")

save.image("./{project_title}.RData")

regress_lsls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(regress_lsls) = names(STAR.data$LoM.raw)

contr_ls = {contr_ls}

# set baseline for regression in parameter(s) of interest
# contr.treatment generates regression coefficients that are like (adjusted mean) ratios of other groups to the baseline group.
# contr.sum generates coefficients that are like (adjusted mean) ratios to average all for all but the mandadory ommitted treatment group (because there is always one fewer independent pairwise comparison than there are pairs).


lm_expr = "{lm_expr}"

rowmask_ls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(rowmask_ls) = names(STAR.data$LoM.raw)

print("RegressMatrix with specified contrasts")

for(i in 1:length(mytypes))
{{
  for(j in which(!mynorms %in% do.not.regress))
  {{
    # set up data matrices
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}
    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp=paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./3/',plotbase=paste(mynorms[j],mytypes[i],'vs',tmp,sep='.'),plottitle=proj.title)
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    normmat[is.na(normmat)] <-0
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels

    # prepare rowmask for heatmap/MDS (remove non-expr or low expr>hclustlim)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) )
    # Too stringent of a rowmask for CFRNA project

    rowmask_ls[[ mytypes[i] ]][[ mynorms[j] ]] = rowmask # save for later

    # regression

    regress_lsls[[mytypes[i]]][[mynorms[j]]] = regressMatrix(normmat[rowmask,], expt.design=annCol[annCol.lm_by], lm_expression=lm_expr, contr_list = contr_ls, plot2file = FALSE, plotdata = plotdata)
  }}
}}

print("RegressMatrix with specified contrasts : Complete")

save.image("./{project_title}.RData")

topn=50 # number with which to play

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress) )
  {{
    ans = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    cat(mytypes[i],mynorms[j],"1-p0:",signif(unlist(lapply(2:length(ans),function(x){{1-ans[[x]]$pi0}})),2),"\n")
    cat(mytypes[i],mynorms[j],"qcut:",signif(unlist(lapply(2:length(ans),function(x){{quantile(ans[[x]]$qvalues,probs=topn/length(ans[[x]]$qvalues))}})),3),
    "returns ~",topn,"out of",length(ans[[2]]$qvalues),"\n")
  }}
}}

# Plots for each design factor (default) or for factors specified in facSel
#   1) Histogram of p values that were included in the design
#      Look for a peak on the left, and no peaks in the middle or on the right
#   2) qvalue's default plots, with full qvalue range c(0,1) plotted
#      Look for descent to pi0 in the top left tuning plot with good asymptote
#      The slow/steep rise in q-values in remaining plots depends on resolving power of data


lm_by = annCol.lm_by
histbins=40
select_lsmk = vector(mode='list',length=length(STAR.data$LoM.raw))
names(select_lsmk) = names(STAR.data$LoM.raw)
print("Perform qcQvalue & Generate masks")

dir.create(file.path(getwd(),'3'),showWarnings=FALSE)

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    normmat[is.na(normmat)] <-0
    colnames(normmat)=samp.labels
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]
    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    # pull out regression design and set up plot config

    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp = paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./3',plotbase=paste(mynorms[j],mytypes[i],tmp,sep='.'),plottitle=proj.title)

    ##for(fac in names(reg_ls$q_list)[grep(lm_by,names(reg_ls$q_list))])

    for(fac in names(reg_ls$q_list))
    {{
      if(!fac %in% '(Intercept)')
      {{
        pdata = plotdata
        pdata$plotbase = paste(plotdata$plotbase,make.names(fac),sep='.')
        ans = qcQvalues(norm_x=normmat[mymask,], pvalue_v=reg_ls$p_mat[,fac], obj_qvalue=reg_ls$q_list[[fac]], attribs=annCol[lm_by],
        oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, histbins=histbins, plot2file=TRUE)
        select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac)] = ans['rowmask']
      }}
    }}
  }}
}}

print("Perform qcQvalue & Generate masks : Complete")

save.image("./{project_title}.RData")

# Build of ratios to baselines for desired design factors
# Also build masks selecting genes based on q-values per factor
# and follow-up considerations such as expression level and fold change
# Loop over q-value cuts to assist with final cut selection.
# Plot heatmaps and MDS plots based on selected genes and designed ratios


ngene_v = c(50,100,200) # q-value cuts by number; can also cut by q-value
ratioby_ls = list("{lm_by}"=contr_ls${lm_by}$baseline)
ratio_fold = 1.3
intensity_fold = 2

cut_ls = list(q_combine="OR", rcut_fold=ratio_fold, icut_fold=intensity_fold)

clustrowmin = 10 # min data rows for heatmap and MDS plots

# save ratios and selections

ratio_lsmat = vector(mode='list',length=length(STAR.data$LoM.raw))
names(ratio_lsmat) = names(STAR.data$LoM.raw)

print("Perform PlotRatios")

dir.create(file.path(getwd(),'4'),showWarnings=FALSE)
for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    normmat[is.na(normmat)] <-0
    colnames(normmat)=samp.labels

    if( grepl('SJ',mytypes[i]))
    {{#not the most robust way to find SJs
      # map SJ positions to overlapping genes
      # specifying defaults to fix "Gene" and Pos as colnames
      gtf.col = c("Gene","Symbol","Chr","start","stop")
      pos.col = "Pos"
      mySJ_dt = mapSJ2feature(STAR.data$LoM.raw[[ mytypes[i] ]][[1]], pos.col=pos.col, gtf.col=gtf.col, gtf.file=gtf_file, gtf.Rdir=genome.func)
      # min data rows for heatmap and MDS plots
      # map gene IDs back to data matrix using mapping built earlier
      # this enables include_ID to select rows

      idx2 = match( rownames(STAR.data$LoM.raw[[ mytypes[i] ]][[1]]), mySJ_dt[,get(pos.col)] )
      idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
      rownames(normmat)[idx1] = mySJ_dt[idx2,get(gtf.col[1])]
    }}

    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    reg_ls = reg_ls[!grepl('Intercept',names(reg_ls))]
    plotdata = list(plotdir='./4',plotbase=paste(mynorms[j],sub('counts','ratios',mytypes[i]),'minratio',paste0(ratio_fold,'x'),'minexpr',round(min(normmat,na.rm=T)+log2(intensity_fold),1),sep='.'),plottitle=proj.title)

    for(fac in names(reg_ls))
    {{
        for( ngenes in ngene_v )
        {{
          # calculate and plot ratios
          cut_ls$qcut = ngenes
          pdata = plotdata
          pdata$plotbase = paste(plotdata$plotbase,'ratio.vs',ratioby_ls${lm_by},ngenes,sep='.')

          # this function makes the ratios and cuts
          ans = designRatios(normmat[mymask,], q_list=reg_ls[[fac]], attribs=annCol, ratioby_ls=ratioby_ls, cut_ls=cut_ls)

          # this function creates the heatmap and MDS plot
          if( sum(ans$rowmask)>clustrowmin )
          {{
            ans2= plotRatios( ratiomat=ans$ratiomat, attribs=annCol[annCol.plotme], oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, rowmask=ans$rowmask, plot2file=TRUE)
          }}
          # save selections
          select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac,ngenes,sep="_")] = ans['rowmask']
        }}
        # save ratiomat (not dependent on ngene cut) after ngene loop
        ratio_lsmat[[mytypes[i]]][[mynorms[j]]] = ans$ratiomat
    }}
  }}
}}

print("Perform PlotRatios : Complete")

save.image("./{project_title}.RData")

print("Perform biomaRt gene annotation")

ratio_mat = ratio_lsmat[["{path_type}"]][["{path_norms}"]]
reg_ls = regress_lsls[["{path_type}"]][["{path_norms}"]]
ID = rownames(ratio_mat)

# assemble signatures

reg_ls = regress_lsls$gene.counts${path_norms}

gtf.Rdir = "{gtf_read_dir}"

out_norm_mat = LoM.norms$gene.counts$loess[rowmask_ls${path_type}${path_norms},]
out_table = outputTable(normmat= out_norm_mat, gtf.file = "{gtf_file}",ratiomat = ratio_lsmat${path_type}${path_norms}, q_list=reg_ls$q_list, gtf.Rdir=genome.func, gtf.key='transcript')
out_table = out_table[!duplicated(out_table$Gene),]

write.table(out_table, row.names = FALSE, file=file.path(getwd(), './tables', paste("{project_title}","Normed_with_Ratio_and_Abundance.txt", sep="_")),quote=FALSE,sep='\t')

print("Generate Abundance and Ratio table with associated q-value and p-values : Complete")

save.image("./{project_title}.RData")
"""


qc_matrix_w_deseq_model = """
#### Running DESeq2 first ####
print('RUNNING DESEq2')
setwd("{results}")

library(DESeq2)
library(data.table)
library(stringr)
readpattern = "{read_pattern}"
filt_pattern = str_replace_all(readpattern,'[.*]','')
tep_mat = round(as.matrix(read.table("{dataset}",check.names=F)))
filt_dt = data.table(tep_mat,check.names=F)
samps = sort(grep(filt_pattern, names(filt_dt),value=T))
cts = as.matrix(filt_dt[,samps,with=F])
row.names(cts) = row.names(tep_mat)

md.dt = read.csv("{meta_file}",row.names=1,sep="\t")
md.dt = md.dt[,{ann_colplotme}]
coldata = md.dt[colnames(cts),]

dds = DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ {lm_by})

dds${lm_by} = relevel(dds${lm_by}, ref = "{baseline}")
dds = DESeq(dds)
vsd <- vst(dds, blind=FALSE)

write.table(assay(vsd), file="{results}/{project_title}_deseq.txt",quote=F,sep="\t")

print('DESEQ2 finished')
#### Done with DESeq2 ####
require(data.table)
require(NMF)
require(affy)
require(limma)
require(biomaRt)
library(stringr)
source("{code_dir}AbundanceFunctions/BiasReduce.R")
source("{code_dir}AbundanceFunctions/ExtractTransformLoad.R")
source("{code_dir}AbundanceFunctions/DifferentialAnalysis.R")
source("{code_dir}AbundanceFunctions/NonVisualOutput.R")
source("{code_dir}GenomicsFunctions/ReadAndParse.R")
source("{code_dir}AssociationFunctions/gs.wrapper.R")
source("{code_dir}AssociationFunctions/PathwayAnalysis.R")
source("{code_dir}BcorePlotting/SummaryPlots.R")
source("{code_dir}BcorePlotting/MultipleTestingCorrection.R")
source("{code_dir}BcorePlotting/ClusteringPlots.R")
library(ggplot2)

dir.create(file.path(getwd(),'{project_title}'),showWarnings=FALSE)
oneclass = {oneclass}
# constants
# max data rows for hclust
hclust.limit = 2^16
annCol.lm_by = {annCollm_by}

# quantile of data distribution reqd in one group's worth of data if too many rows for hclust()
hc.frac.cut = 0.75;
SJ.counts.na.frac = 0.25;
# max fraction of samples not having detected a splice junction for the splice
# junction to be retained in raw data


# regression parameters
na.lim = 0 # max NAs per row tolerated by lm() at least in some cases
do.not.regress = "alograw" # control norm not to be used for regression stats
# plotting colors
colors.rgb = c(rgb(0,0,0),rgb(0.1,0.1,1),rgb(0,.7,.7),rgb(0,.7,0),rgb(.7,1,0),rgb(.7,0,.7))
md.file = "{meta_file}"
md.orientation = "byRow" # sample_ids are in @ row. alt:byCol (IDs in @ col)
md.IDcol = "{sample_id}" # reqd if md.orientation is byRow; byCol==headers are IDs

# gene annotation
taxID = {tax_id}
gene2ENSfile = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/ncbi/gene2ensembl.gz"
gene2ENS.col = c("taxID","EntrezID","Gene","RefSeqTranscript","EnsemblTranscript","RefSeqProtein","EnsemblProtein")
gtf_file = "{gtf_file}"
gtf.feature = "{gtf_feature}"
gtf.orig.col = c("gene_id","gene_name","gene_biotype")
gtf.col = c("Gene","Symbol","biotype")
proj.title = "{project_title}"
gtf.Rdir = "{gtf_read_dir}"
genome.func = "{code_dir}GenomicsFunctions"


readir = "{read_dir}"
readpattern = "{read_pattern}"
useme.cols = "{useme_cols}"
label.from.colname = "{label_from_colname}"
samps = dir(path=readir, pattern=readpattern)
samp.labels = gsub(label.from.colname,'\\\\1', samps)

annCol.names = "group"
annCol.label = "{label_from_colname}"

# read in STAR alignments

print("Read in Abundance data")

logged_B = {logged_B}

filt_pattern = str_replace_all(readpattern,'[.*]','')
tep_mat = as.matrix(read.table("{dataset}",check.names=F))
filt_dt = data.table(tep_mat,check.names=F)
samps = sort(grep(filt_pattern, names(filt_dt),value=T))
samp.labels = samps
filt_mat = as.matrix(filt_dt[,samps,with=F])
row.names(filt_mat) = row.names(tep_mat)
STAR.data = list(LoM.raw=list(gene.counts=list(Reads=filt_mat)),expt.design=list(group=samps), myreads = 'Reads')
attr(STAR.data$expt.design,'names')='group'
attr(STAR.data$expt.design,'lm_by')=annCol.lm_by

myreads = STAR.data$myreads

print("Filter gene.counts data")

{{
    if(any( names(STAR.data$LoM.raw)=="gene.counts"))
    {{
      STAR.data$LoM.orig.raw = STAR.data$LoM.raw$gene.counts
      for(tag in names(STAR.data$LoM.raw$gene.counts))
      {{
        STAR.data$LoM.raw$gene.counts[[tag]] = STAR.data$LoM.raw$gene.counts[[tag]][rowSums(STAR.data$LoM.raw$gene.counts[[tag]])>0, ]
      }}
    }}
}}

print("Filtering gene.counts : Complete")

save.image("./{project_title}.RData")

LoM.norms = vector(mode='list',length=length(STAR.data$LoM.raw))
names(LoM.norms) = names(STAR.data$LoM.raw)


print("Bias reduction with normMatrix")
dir.create(file.path(getwd(),'{project_title}/tables'),showWarnings=FALSE)

for(tag in names(STAR.data$LoM.raw))
{{
  if( length(STAR.data$LoM.raw[[tag]])>1 )
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[myreads]],fileData_ls=list(fileDir='./tables/',fileBase='{project_title}'),logged_B=logged_B)
  }}
  else
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[1]],fileData_ls=list(fileDir='./tables/',fileBase='{project_title}'), logged_B=logged_B)
  }}
}}

des_idx = row.names(LoM.norms$gene.counts$loess)
de_mat = as.matrix(read.table("{results}/{project_title}_deseq.txt",check.names=F))
de_mat = de_mat[des_idx,]

LoM.norms$gene.counts$Deseq  = de_mat

print("Bia reduction with normMatrix : Complete")

save.image("./{project_title}.RData")
## Iterate through bias reduced mats and write to working directory


## custom: reannotate from metadata
## run if necessary metadata is in file and not also in FASTQ file names
# map annotation to read matrix
md.dt = fread("{meta_file}")

print("Evaluate meta-data file for sample name consistency")

{{
    if( md.orientation == "byCol" ){{ # samples are one per column
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), names(md.dt) )
    }}
    else
    {{ # samples are one per row
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), md.dt[,get(md.IDcol)] )
    }}
}}

idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
{{
    if( sum(idx1==1:ncol(STAR.data$LoM.raw[[1]][[myreads]]))==ncol(STAR.data$LoM.raw[[1]][[myreads]]))
    {{
      if( md.orientation == "byCol" )
      {{ # samples are one per column
        annCol = NULL
        namecol = setdiff( 1:ncol(md.dt), idx2 )
        md.factors = as.character(md.dt[,namecol,with=F])
        for(k in 1:nrow(md.dt) )
        {{
          annCol[[ md.dt[k,get(names(md.dt)[namecol])] ]][idx1] = as.vector( md.dt[ k, mget(names(md.dt)[idx2]) ] )
        }}
      }}
      else
      {{ # samples are one per row
        annCol = NULL
        namecol = setdiff(names(md.dt), md.IDcol)
        md.factors = as.character(namecol)
        for( k in setdiff(names(md.dt),md.IDcol))
        {{
          annCol[[ k ]][idx1] = as.vector( md.dt[ idx2, get(k) ] )
        }}
      }}
    }}
    else
    {{
      stop(paste("Some samples have no annotation in",md.file))
    }}
}}

print("Evaluation of meta-data : Complete")

mytypes = names(LoM.norms)
mynorms = names(LoM.norms[[1]])

grpBy = annCol[[oneclass]]
annCol.plotme = {ann_colplotme}
clim.pct=0.96
histbins=40
dir.create(file.path(getwd(),'{project_title}/1'),showWarnings=FALSE)
dir.create(file.path(getwd(),'{project_title}/2'),showWarnings=FALSE)

print("Generate summary and qc.cluster plots")

#### Loop through mynorms ######

for(i in 1:length(mytypes))
{{
  for(j in 1:length(mynorms))
  {{
    # set up matrices and config for plotting
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}

    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels
    plotdata = list(plotdir='./1',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)))

    if(sum(rowmask)>hclust.limit)
    {{
      rowmask = rowSums(rawmat>=quantile(rawmat,probs=hc.frac.cut,na.rm=T),na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)) )
    }}

    ans = summary.plots(rawmat=log2(rawmat +1), normmat=normmat, mynorm=mynorms[j], samp.labels=samp.labels, samp.classes=grpBy, plotdata=plotdata,plot2file=TRUE,histbins=histbins, colorspec=colors.rgb)
    plotdata = list(plotdir='./2',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    ans = qc.clusters(rawmat=log2(rawmat[rowmask,] +1), normmat=normmat[rowmask,], attribs=annCol[annCol.plotme], oneclass=oneclass, colorspec=colors.rgb, plotdata=plotdata, plot2file=TRUE, clim.pct=clim.pct)
  }}
}}

print("Generate summary and qc.cluster plots : Complete")

save.image("./{project_title}.RData")

regress_lsls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(regress_lsls) = names(STAR.data$LoM.raw)

contr_ls = {contr_ls}

# set baseline for regression in parameter(s) of interest
# contr.treatment generates regression coefficients that are like (adjusted mean) ratios of other groups to the baseline group.
# contr.sum generates coefficients that are like (adjusted mean) ratios to average all for all but the mandadory ommitted treatment group (because there is always one fewer independent pairwise comparison than there are pairs).


lm_expr = "{lm_expr}"

rowmask_ls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(rowmask_ls) = names(STAR.data$LoM.raw)

print("RegressMatrix with specified contrasts")

for(i in 1:length(mytypes))
{{
  for(j in which(!mynorms %in% do.not.regress))
  {{
    # set up data matrices
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}
    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp=paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./3/',plotbase=paste(mynorms[j],mytypes[i],'vs',tmp,sep='.'),plottitle=proj.title)
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels

    # prepare rowmask for heatmap/MDS (remove non-expr or low expr>hclustlim)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) <= na.lim
    rowmask_ls[[ mytypes[i] ]][[ mynorms[j] ]] = rowmask # save for later

    # regression

    regress_lsls[[mytypes[i]]][[mynorms[j]]] = regressMatrix(normmat[rowmask,], expt.design=annCol[annCol.lm_by], lm_expression=lm_expr, contr_list = contr_ls, plot2file = FALSE, plotdata = plotdata)
  }}
}}

print("RegressMatrix with specified contrasts : Complete")

save.image("./{project_title}.RData")

topn=500 # number with which to play

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress) )
  {{
    ans = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    cat(mytypes[i],mynorms[j],"1-p0:",signif(unlist(lapply(2:length(ans),function(x){{1-ans[[x]]$pi0}})),2),"\n")
    cat(mytypes[i],mynorms[j],"qcut:",signif(unlist(lapply(2:length(ans),function(x){{quantile(ans[[x]]$qvalues,probs=topn/length(ans[[x]]$qvalues))}})),3),
    "returns ~",topn,"out of",length(ans[[2]]$qvalues),"\n")
  }}
}}

# Plots for each design factor (default) or for factors specified in facSel
#   1) Histogram of p values that were included in the design
#      Look for a peak on the left, and no peaks in the middle or on the right
#   2) qvalue's default plots, with full qvalue range c(0,1) plotted
#      Look for descent to pi0 in the top left tuning plot with good asymptote
#      The slow/steep rise in q-values in remaining plots depends on resolving power of data


lm_by = annCol.lm_by
histbins=40
select_lsmk = vector(mode='list',length=length(STAR.data$LoM.raw))
names(select_lsmk) = names(STAR.data$LoM.raw)
print("Perform qcQvalue & Generate masks")

dir.create(file.path(getwd(),'{project_title}/3'),showWarnings=FALSE)

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    normmat[is.na(normmat)] <- 0
    colnames(normmat)=samp.labels
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]
    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    # pull out regression design and set up plot config

    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp = paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./3',plotbase=paste(mynorms[j],mytypes[i],tmp,sep='.'),plottitle=proj.title)

    ##for(fac in names(reg_ls$q_list)[grep(lm_by,names(reg_ls$q_list))])

    for(fac in names(reg_ls$q_list))
    {{
      if(!fac %in% '(Intercept)')
      {{
        pdata = plotdata
        pdata$plotbase = paste(plotdata$plotbase,make.names(fac),sep='.')
        ans = qcQvalues(norm_x=normmat[mymask,], pvalue_v=reg_ls$p_mat[,fac], obj_qvalue=reg_ls$q_list[[fac]], attribs=annCol[lm_by],
        oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, histbins=histbins, plot2file=TRUE)
        select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac)] = ans['rowmask']
      }}
    }}
  }}
}}

print("Perform qcQvalue & Generate masks : Complete")

save.image("./{project_title}.RData")

# Build of ratios to baselines for desired design factors
# Also build masks selecting genes based on q-values per factor
# and follow-up considerations such as expression level and fold change
# Loop over q-value cuts to assist with final cut selection.
# Plot heatmaps and MDS plots based on selected genes and designed ratios


ngene_v = c(200,500,1000) # q-value cuts by number; can also cut by q-value
ratioby_ls = list("{lm_by}"=contr_ls${lm_by}$baseline)
ratio_fold = 1.0
intensity_fold = 1

cut_ls = list(q_combine="OR", rcut_fold=ratio_fold, icut_fold=intensity_fold)

clustrowmin = 10 # min data rows for heatmap and MDS plots

# save ratios and selections

ratio_lsmat = vector(mode='list',length=length(STAR.data$LoM.raw))
names(ratio_lsmat) = names(STAR.data$LoM.raw)

print("Perform PlotRatios")

dir.create(file.path(getwd(),'{project_title}/4'),showWarnings=FALSE)
for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    normmat[is.na(normmat)] <-0
    colnames(normmat)=samp.labels
    groupdat = data.frame(annCol)
    sample_id = groupdat$Sample
    colorfactor= "{lm_by}"
    #png(filename=paste("./{project_title}/4/", mynorms[j], "_PCA_Plot.png", sep=""))
    gene_pcaplot(normmat, sample_id, groupdat=groupdat, colorfactor=colorfactor, shapefactor=NULL, plot_sample_ids=TRUE, pcnum=1:2, plottitle = paste(mynorms[j]," PCA Plot",sep=""))
    ggsave(paste("./{project_title}/4/", mynorms[j], "_PCA_Plot.png", sep=""))
    if( grepl('SJ',mytypes[i]))
    {{#not the most robust way to find SJs
      # map SJ positions to overlapping genes
      # specifying defaults to fix "Gene" and Pos as colnames
      gtf.col = c("Gene","Symbol","Chr","start","stop")
      pos.col = "Pos"
      mySJ_dt = mapSJ2feature(STAR.data$LoM.raw[[ mytypes[i] ]][[1]], pos.col=pos.col, gtf.col=gtf.col, gtf.file=gtf_file, gtf.Rdir=genome.func)
      # min data rows for heatmap and MDS plots
      # map gene IDs back to data matrix using mapping built earlier
      # this enables include_ID to select rows

      idx2 = match( rownames(STAR.data$LoM.raw[[ mytypes[i] ]][[1]]), mySJ_dt[,get(pos.col)] )
      idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
      rownames(normmat)[idx1] = mySJ_dt[idx2,get(gtf.col[1])]
    }}

    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    reg_ls = reg_ls[!grepl('Intercept',names(reg_ls))]
    plotdata = list(plotdir='./4',plotbase=paste(mynorms[j],sub('counts','ratios',mytypes[i]),'minratio',paste0(ratio_fold,'x'),'minexpr',round(min(normmat,na.rm=T)+log2(intensity_fold),1),sep='.'),plottitle=proj.title)

    for(fac in names(reg_ls))
    {{
        for( ngenes in ngene_v )
        {{
          # calculate and plot ratios
          cut_ls$qcut = ngenes
          pdata = plotdata
          pdata$plotbase = paste(plotdata$plotbase,'ratio.vs',ratioby_ls${lm_by},ngenes,sep='.')

          # this function makes the ratios and cuts
          ans = designRatios(normmat[mymask,], q_list=reg_ls[[fac]], attribs=annCol, ratioby_ls=ratioby_ls, cut_ls=cut_ls)

          # this function creates the heatmap and MDS plot
          if( sum(ans$rowmask)>clustrowmin )
          {{
            ans2= plotRatios( ratiomat=ans$ratiomat, attribs=annCol[annCol.plotme], oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, rowmask=ans$rowmask, plot2file=TRUE)
          }}
          # save selections
          select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac,ngenes,sep="_")] = ans['rowmask']
        }}
        # save ratiomat (not dependent on ngene cut) after ngene loop
        ratio_lsmat[[mytypes[i]]][[mynorms[j]]] = ans$ratiomat
    }}
  }}
}}

print("Perform PlotRatios : Complete")

save.image("./{project_title}.RData")

print("Perform biomaRt gene annotation")

ratio_mat = ratio_lsmat[["{path_type}"]][["{path_norms}"]]
reg_ls = regress_lsls[["{path_type}"]][["{path_norms}"]]
ID = rownames(ratio_mat)

# assemble signatures

reg_ls = regress_lsls$gene.counts${path_norms}

gtf.Rdir = "{gtf_read_dir}"

out_norm_mat = LoM.norms$gene.counts$loess[rowmask_ls${path_type}${path_norms},]
out_table = outputTable(normmat= out_norm_mat, gtf.file = "{gtf_file}",ratiomat = ratio_lsmat${path_type}${path_norms}, q_list=reg_ls$q_list, gtf.Rdir=genome.func, gtf.key='transcript')
out_table = out_table[!duplicated(out_table$Gene),]

write.table(out_table, row.names = FALSE, file=file.path(getwd(), './tables', paste("{project_title}","Normed_with_Ratio_and_Abundance.txt", sep="_")),quote=FALSE,sep='\t')

print("Generate Abundance and Ratio table with associated q-value and p-values : Complete")

save.image("./{project_title}.RData")
"""



deseq_model = """
library(DESeq2)
library(ggplot2)
readpattern = "{read_pattern}"
filt_pattern = str_replace_all(readpattern,'[.*]','')
tep_mat = round(as.matrix(read.table("{dataset}",check.names=F)))
filt_dt = data.table(tep_mat,check.names=F)
samps = sort(grep(filt_pattern, names(filt_dt),value=T))
cts = as.matrix(filt_dt[,samps,with=F])

md.dt = read.csv("{meta_file}",row.names=1,sep="    ")
md.dt = md.dt[,{ann_colplotme}]
coldata = md.dt[colnames(cts),]

dds = DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ {lm_by})
                              
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

dds${lm_by} = relevel(dds${lm_by}, ref = "{baseline}")
dds = DESeq(dds)
res = results(dds)

dds_cts = counts(dds, normalized=TRUE)
write.table(dds_cts, file="{results}/{project_title}_deseq.txt",quote=F,sep="   ")
write.table(res, file="{results}/{project_title}_deseq_results.txt",quote=F,sep="   ")
"""
