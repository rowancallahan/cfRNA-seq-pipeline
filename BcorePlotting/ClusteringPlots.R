#####################################################################################################
# Set of functions useful for plotting clustering-based visual data
#
# gene_pcaplot     --- Generates a PCA plot with labeled axes
# plotRatios       --- Plot heatmap and MDS for ratio data
# makeMDSplot      --- MDS plot child function called by qc functions
# makeHeatmap      --- heatmap child function called by qc & design functions
# getCidFromHeatmap--- Extract dendrogram info from heatmap return values
#
# Author: Julja Burchard, Theresa Lusardi, Jessica Minnier, Mark Fisher, and Wes Horton
# Started - July, 2016
#####################################################################################################
library(ggplot2)
library(RColorBrewer)

gene_pcaplot <-
function(exprdat,sampleid,groupdat=NULL,colorfactor=NULL,shapefactor=NULL, plot_sampleids=TRUE, pcnum=1:2, plottitle = "PCA Plot") {
  #borrowed from Jessica Minnier, 2016
  #adapted from DESeq2:::plotPCA.DESeqTransform
  pca <- prcomp(t(exprdat))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if(is.null(groupdat)){ groupdat = data.frame("group"=rep(1,ncol(exprdat)))}
  intgroup = colnames(groupdat)
  allgroup <- if (length(intgroup) > 1) {
    factor(apply(groupdat, 1, paste, collapse = ":"))
  }else{allgroup <- intgroup}
  d <- data.frame(PC1 = pca$x[, pcnum[1]], PC2 = pca$x[, pcnum[2]], uniquegroup = allgroup,groupdat, name = sampleid)
  percentVar <- round(100 * percentVar)
  if(is.null(colorfactor)) {d$color=as.factor(1)}else{
    colnames(d)[colnames(d)==colorfactor] <- "color"}
  if(is.null(shapefactor)) {d$shape=as.factor(1)}else{
    colnames(d)[colnames(d)==shapefactor] <- "shape"
  }
  if(identical(shapefactor,colorfactor)) {d$shape = d$color}
  p <- ggplot(d, aes(PC1, PC2, color=color, shape=shape, size=3)) 
  if(plot_sampleids) {
    p <- p + geom_text(aes(label=name,size=10))
  }else{p <- p + geom_point()}
  
  if(!is.null(colorfactor)) {
    p <- p + guides(color=guide_legend(title=colorfactor))
  }else {
    p <- p + guides(color = "none")
  }
  if(!is.null(shapefactor)) {
    p <- p + guides(shape=guide_legend(title=shapefactor))
  }else{
    p <- p + guides(shape = "none")
  }
  p <- p + guides(size= "none") + theme_bw() + 
    xlab(paste0("PC",pcnum[1],": ",percentVar[pcnum[1]],"% variance")) +
    ylab(paste0("PC",pcnum[2],": ",percentVar[pcnum[2]],"% variance")) + ggtitle(plottitle)
  
  return(p)
}
plotRatios <-
function (ratiomat, attribs, oneclass, plotdata, colorspec, 
                       rowmask=NULL, tag="selected",
                       NAtol=ncol(ratiomat)/2, SDtol=0.01, 
                       heatmap_plot=TRUE, MDS_plot=TRUE, 
                       clim.pct=0.99, clim_fix=NULL, 
                       plot2file = FALSE, filesep='/', 
                       annColors = NA, annRow = NA,
                       cexRow=0.00001, png_res=1200) {
  #' Differential expression visualizations. Heatmap and MDS plot of ratios
  #' @description 
  #' This function takes a matrix of ratios and displays them graphically in a heatmap and/or MDS plot
  #' @param ratiomat matrix of ratio data in columns (with headers!)
  #' @param rowmask optional mask of rows to plot
  #'   default is to create an all-TRUE rowmask of length nrow(ratiomat)
  #' @param NAtol max number of NAs tolerated; rows with more are masked out
  #' @param SDtol min SD within row tolerated; rows with less are masked out
  #' @param attribs list of sample classifications to be tracked in clustering
  #'   each list element contains a string vector with one label per sample
  #' @param tag word or phrase to indicate what's special about _these_ ratios
  #' @param oneclass string name of attribs element to be used in MDS plot
  #' @param heatmap_plot boolean to plot (TRUE) or not (FALSE) heatmap
  #' @param MDS_plot boolean to plot (TRUE) or not (FALSE) MDS plot
  #' @param plotdata list of info relevant to labeling and saving the plot
  #'   plotdir  plot destination directory
  #'   plotbase  base filename for the plot
  #'   plottitle  title for all plots
  #' @param colorspec is a vector of color specifiers for colorRampPalette  
  #' @param clim.pct  0:1 - fraction of data to limit max color
  #' @param clim_fix if set, max abs ratio to show; data>clim_fix are shown as clim_fix
  #' @param annRow optional named lists of row annotations; can be single list
  #' @param annColors  optional named lists of annotation colors;
  #'    names should match names in annRow and attribs
  #' @param cexRow rowlabel size for heatmap, set to ~invisible by default
  #' @param png_res resolution of saved png files in dpi; default 300
  #' @export

  # imports
  require(NMF)

  # arguments

  # test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }

  # check for rowmask; create if NULL
  if( is.null(rowmask) ){
    rowmask = logical(nrow(ratiomat)); rowmask[]= TRUE
  }
  # adjust rowmask if needed to keep problem rows out of plots
  NAmask = rowSums(is.na(ratiomat)) <= NAtol
  SDmask = sapply(1:nrow(ratiomat), 
                  function(x){sd(ratiomat[x,],na.rm=T)}) >= SDtol
  SDmask[ is.na(SDmask) ] = FALSE # in case of all-NA rows
  rowmask = rowmask & NAmask & SDmask
  # report selection size
  message(sprintf('selected rows = %s, excess NA rows = %s, low SD rows = %s',
                   sum(rowmask), sum(!NAmask), sum(!SDmask) ))

  ah_ls = NULL
  if( heatmap_plot ){
    plotID = '5q1'
    plotDesc = paste('Heatmap', tag, sep="_")
    if(plot2file) {
    png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc), width=5,height=7,units="in",res=png_res)
    }
    
    # plot heatmap
    ah_ls = makeHeatmap( ratiomat=ratiomat[rowmask,], attribs=attribs, plottitle = plotdata$plottitle, clim.pct=clim.pct, clim_fix=clim_fix, cexRow=cexRow, annColors = annColors, annRow = annRow)
    
    if(plot2file) dev.off()
  }


  obj_MDS = NULL
  if( MDS_plot ){
    plotID = '6q1'
    plotDesc = paste('MDS', tag, sep="_") 
    if(plot2file) {
    png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc), width=5.4,height=5.4,units="in",res=png_res)
    }
  
    # plot MDS
    obj_MDS = makeMDSplot(normmat=ratiomat[rowmask,], attribs=attribs, oneclass=oneclass, colorspec=colorspec, plottitle=plotdata$plottitle, ngenes=sum(rowmask))
  
    if(plot2file) dev.off()
  }


  # return processed data
  invisible( list( ah_ls=ah_ls, obj_MDS=obj_MDS) )
  
}

# Map colors to plotting factor and samples
assignMappingSpecs <- function(attribs, oneclass, colorspec, varPoints_v=NULL){
    #' Map colors to plotting factor and samples
    #' @description
    #' This function is called from within makeMDSplot, appropriately colors vectors and changes point size based on input
    #' @param attribs list of sample classifications to be tracked in clustering. Each list element contains a string vector with one lable per sample.
    #' @param oneclass string name of attribs element to be used in MDS plot
    #' @param colorspec vector of color specifiers for colorRampPalette. colors are generated per sample by their clustering variable values
    #' @param  varPoints_v vector of values to scale points by. Default is no scaling
    #' @return list of vectors of plotting paramters
    #' @export
    
  # Get factor to plot
  sampClasses_v = attribs[[oneclass]]
  
  # Set up plotting colors based on classes
  uSampClasses_v = unique(sampClasses_v)
  colorMap = colorRampPalette(colorspec)
  plotColors_v = colorMap(length(uSampClasses_v))
  # Map factor-based colors back to individual samples
  plotColors_v = plotColors_v[as.numeric(as.factor(sampClasses_v))]
  # R doesn't re-sort on unique: saves order of occurrence
  uPlotColors_v = unique(plotColors_v)
  # Set point size otuput. Either all points are 0.6, or points are scaled.
  if (is.null(varPoints_v)){
    plotCex_v = 0.6
  } else {
    # Size range of values
    sizeRange_v = range(varPoints_v)
    # Get step size
    stepSize_v = (sizeRange_v[2] - sizeRange_v[1]) / 14
    # Scale
    scaledVarPoints_v = sapply(varPoints_v, function(x) round((x-sizeRange_v[1])/stepSize_v)) + 1
    # Transform to cex values
    plotCex_v = ((scaledVarPoints_v / 10) + 0.4)
  } # fi
    
  # Return
  return(list("sampClasses_v" = sampClasses_v, 
              "uSampClasses_v" = uSampClasses_v, 
              "plotColors_v" = plotColors_v, 
              "uPlotColors_v" = uPlotColors_v,
              "plotCex_v" = plotCex_v))
  
} # assignMappingSpecs

# Replicate color mapping with point types and legend point types
assignLabelSpecs <- function(extraParams_ls, sampClasses_v, uSampClasses_v, normmat){
  #' Replicate color mapping with point types and assign legend point types
  #' @description 
  #' This function is called from within makeMDSplot, appropriately assigns point size for plot and legend.
  #' @param extraParams_ls list of extra plotting parameters. Name of element must be exact name of default plotting argument
  #' @param sampClasses_v vector of length == ncol(normmat) containing grouping identifier for each sample
  #' @param uSampClasses_v vector of each unique grouping identifier
  #' @param normmat data matrix, with unique & informative names
  #' @return list of vectors of plotting parameters
  #' @export

  # Assign colnames to labels (default), or assign shapes to sample classes (will mirror colors)
  if ("pch" %in% names(extraParams_ls)){
    # pch is ignored in plots if plot labels is non-null
    plotLabels_v = NULL
    # Get pch portion of extraParams and map to samples
    pointType_v = extraParams_ls$pch
    names(pointType_v) = uSampClasses_v
    extraParams_ls$pch = sapply(sampClasses_v, function(x) x = pointType_v[x], USE.NAMES = F)
    # Legend needs to be updated also
    pchLegend_v = pointType_v
  } else {
    # Set to default column names and legend point
    plotLabels_v = colnames(normmat)
    extraParams_ls$pch = NULL # ignored due to plotLabels_v assignment, but need for consistency
                              # changed from pch = NULL to extraParams_ls$pch = NULL. I don't think
                              # it's necessary, but may be more informative.
    pchLegend_v = 16
  } # fi
  # Return
  return(list("pch" = extraParams_ls$pch,
              "plotLabels_v" = plotLabels_v,
              "pchLegend_v" = pchLegend_v))
} # assignLabelSpecs

createLegend <- function(legendPos_v, uSampClasses_v, uPlotColors_v, pchLegend_v){
  #' Create legend to display on plots
  #' @description 
  #' This function is currently called from within makeMDSplot, although it may be suitable to other
  #' functions as well, as-is or with slight tweaking.
  #' @param legendPos_v "auto" means legend in default configuration at upper limit of user space, 
  #' can specify standard options (see ?legend)
  #' @param uSampClasses_v vector of each unique grouping identifier
  #' @param uPlotcolors_v: vector of colors, one for each unique entry in uSampClasses_v
  #' (should be created using colorRampPalette - see assignMappingSpecs)
  #' @param pchLegend_v: vector of integers specifying point types. Should be either length of 1 or
  #' length(pchLegend_v) == length(uSampClasses_v)
  #' @return legend object
  #' @export

  standardArgs_v = list(legend=uSampClasses_v, col=uPlotColors_v, pch=pchLegend_v, cex=.6)
  # If position not specified, use default
  if (legendPos_v == "auto"){
    axl = par("usr") # c(x1,x2,y1,y2) == extremes of user coords in plot region
    obj_legend = do.call(legend, c(x = axl[2]-.025*abs(diff(axl[1:2])), y = axl[4], standardArgs_v))
  } else {
    # Move legend to specified area
    obj_legend = do.call(legend, c(x=legendPos_v, bg = "transparent", bty = "n", standardArgs_v))
  } # fi
  # Return
  return(obj_legend)
} # createLegend

makeMDSplot <-
function (normmat, attribs, oneclass, colorspec, plottitle, 
                        subtitle=NULL, ngenes=NULL, legendPos_v="auto", varPoints_v=NULL, ...) {
#' Make multi-dimensional scaling plot
#' @description 
#' This function considers only a single clustering variable in coloring MDS plot... Will need to be gneralized
#' Uses the matrix values to cluster the data
#' @param  normmat  data matrix, with unique & informative colnames 
#' @param  attribs  list of sample classifications to be tracked in clustering
#'  each list element contains a string vector with one label per sample
#' @param  oneclass string name of attribs element to be used in MDS plot
#' @param  colorspec vector of color specifiers for colorRampPalette  
#' colors are generated per sample by their clustering variable values
#' @param  plottitle  title for all plots
#' @param  subtitle optional subtitle to add below title on this plot
#' @param  ngenes number of "top" genes for plotMDS to select for distance calculation
#' @param  legendPos_v "auto" means legend in default configuration at upper limit of user space, can specify
#' standard options (see ?legend)
#' @param varPoints_v vector of values to scale points by. Default is no scaling. 
#' If used, values are placed into 15 bins and assigned cex values from 0.5 to 1.9
#' @param  ... any extra plotting parameters to pass to plotMDS. (change plotting points, axis labels, etc.)
#' @return mds plot
#' @export

  # imports
  require(limma)

  # number of "top" genes used for sample-sample distance calculation
  if( is.null(ngenes) ) { ngenes = nrow(normmat) }
  
  # Grab extra arguments specified in ...
  extraParams_ls = list(...)
  
  # Map colors to plotting factor and samples
  mappingSpecs_lsv = assignMappingSpecs(attribs, 
                                        oneclass,
                                        colorspec,
                                        varPoints_v)
  
  # Replicate color mapping with point types and legend point types
  labelSpecs_lsv <- assignLabelSpecs(extraParams_ls, 
                                     mappingSpecs_lsv$sampClasses_v, 
                                     mappingSpecs_lsv$uSampClasses_v, 
                                     normmat)
  
  # Update original input with assignLabelSpecs output
  extraParams_ls$pch = labelSpecs_lsv$pch

  # plot MDS and capture numeric results
  par( mar=c(5,4,4,4)+0.1 ) #default mar=c(5,4,4,2)+0.1; margin in lines
  # Concatenate all the standard arguments into a vector
  standardArgs_v = list(normmat, 
                        col=mappingSpecs_lsv$plotColors_v,
                        labels=labelSpecs_lsv$plotLabels_v,
                        top=ngenes,
                        main=plottitle,
                        cex=mappingSpecs_lsv$plotCex_v)
  # Create MDS object with all standard arguments and any extra arguments, if exist.
  obj_MDS = do.call(plotMDS, c(standardArgs_v, extraParams_ls))
  
  # Add subtitle, if specified
  if( !is.null(subtitle) ){ mtext(subtitle, cex=.8) }

  # allow plotting anywhere on device, and widen right margin
  par(xpd=NA) 
  
  # Create and apply legend
  createLegend(legendPos_v, mappingSpecs_lsv$uSampClasses_v, mappingSpecs_lsv$uPlotColors_v, labelSpecs_lsv$pchLegend_v)
  
  return(obj_MDS)
} # makeMDSplot

setColSpecs <- function(ratiomat, attribs, setCol_v=NULL, colOrder_v=NULL, labcoltype=c("colnames","colnums")){
  #' Set column specs for heatmap
  #' @description 
  #' This function is called in makeHeatmap and determines column ordering and labelling
  #' @param ratiomat  ratio data matrix, optimally derived from normmat,
  #'    with unique & informative colnames
  #' @param attribs  list of sample classifications to be tracked in clustering
  #'    each list element contains a string vector with one label per sample
  #'    set to NA to omit
  #' @param setCol_v name of attrib to sort matrix columns by (will turn off column-clustering in heatmap output)
  #'      Defaults to NULL to keep original clustering (hierarchical)
  #' @param colOrder_v Vector containing all unique values of attribs[[setColv]] in desired output order. Defaults
  #'      to NULL so that if setColv is specified alone, output order will be alphabetical.
  #' @param labcoltype colnames to show ratiomat column names, colnums to show col #s
  #' @return list of column parameters
  #' @export

  # TODO: incorporate normmat as well.
  
  # Set column labels
  if (any(grepl("colnames", labcoltype,ignore.case = T))){
    labCol_v = colnames(ratiomat)
  } else {
    labCol_v = 1:ncol(ratiomat)
    colnames(ratiomat) <- labCol_v
  } # fi
  
  # Exit here if attribs is NA
  if (unique(is.na(attribs))){
    setCol_v=NULL
    return(list("ratiomat" = ratiomat, "attribs" = attribs, "setCol_v" = setCol_v, "labCol_v" = labCol_v))
  } else {
    # If specified sort ratiomat by attrib name found in setCol_v
    if (is.null(setCol_v)){
      setCol_v=NULL # Default is to use clustering
    } else {
      # Combine column indeces with ordering attribute in a data.table TODO: multiple attribs
      colIndex_v = 1:ncol(ratiomat)
      temp_dt = as.data.table(cbind("Index" = colIndex_v, "sortBy" = attribs[[setCol_v]]))
      
      # Sort
      if (is.null(colOrder_v)){
        setkey(temp_dt, `sortBy`) # Default is alphabetical
      } else {
        temp_dt <- temp_dt[order(match(`sortBy`, colOrder_v))] # Sort by given vector (Note: can be same length of ncol, or can be 1 for each unique attrib)
      } # fi
      
      # Order everything by sortBy
      colOut_v = as.numeric(temp_dt[,`Index`]) # Extract from sorted data.table
      attribs[[setCol_v]] <- temp_dt[,`sortBy`] # Reorder attribs
      ratiomat = ratiomat[,colnames(ratiomat)[colOut_v]] # Reorder ratiomat
      labCol_v = labCol_v[colOut_v] # Reorder column labels
      setCol_v = NA # Turns off ordering in aheatmap call, which is what we want b/c we ordered mat appropriately
    } # fi
    # Output
    return(list("ratiomat" = ratiomat, "attribs" = attribs, "setCol_v" = setCol_v, "labCol_v" = labCol_v))
  } # fi
} # setColSpecs

uni_qual_colors <- function(attribs, pal = 'Set1'){
    #' @description
	#' This function expands qualitative ColorBrewer palettes to accommodate datasets which large
	#' number of factors that cannot be represented by a single ColorBrewer palette.
    #' @param attribs list of sample classifications to be tracked in clustering
    #'    each list element contains a string vector with one label per sample
    #'    set to NA to omit
    require(RColorBrewer)
	colourCount = 0

	`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

	for(i in names(annCol[annCol.plotme])){
	colourCount %+=% length(unique(annCol[i][[1]]))[[1]]; colourCount}

	getPalette = colorRampPalette(brewer.pal(9, pal))

	ann_color_palette = getPalette(colourCount)
	subcolors = ann_color_palette
	annColors = list()

	for(i in names(attribs)){
		icolors = subcolors[1:length(unique(attribs[i][[1]]))[[1]]]
		names(icolors) = unique(attribs[i][[1]])
		annColors[i]=list(icolors)
		subcolors=subcolors[-(1:length(unique(attribs[i][[1]])))]
		}
	return(annColors)
}


makeHeatmap <-
function (ratiomat, attribs, plottitle, subtitle=NULL, normmat=NULL,
                        clim.pct=.99, clim_fix=NULL, colorbrew="-PiYG:64", 
                        cexRow=0.00001, annRow = NA, annColors = NA,
                        cexCol=min(0.2 + 1/log10(ncol(ratiomat)), 1.2),
                        labcoltype=c("colnames","colnums") , labRowType = NULL,
                        setCol_v=NULL, setRowv=NULL, colOrder_v=NULL, inclLegend = TRUE) {
  #' Make heatmap of ratio data
  #' @description 
  #' This function makes a heatmap of ratio data, displaying experimental design values as tracks
  #' Uses the matrix values to cluster the data
  #' @param normmat  abundance data matrix, optionally used to set color limits
  #'   set to null to use ratiomat for color limits
  #' @param ratiomat  ratio data matrix, optimally derived from normmat,
  #'    with unique & informative colnames 
  #' @param attribs list of sample classifications to be tracked in clustering
  #'    each list element contains a string vector with one label per sample
  #'    set to NA to omit
  #' @param colorbrew is a colorBrewer string specifying a heatmap color scale 
  #'    colors are generated per sample by their clustering variable values
  #' @param plottitle  title for all plots
  #' @param clim.pct  0:1 - fraction of data to limit max color
  #' @param clim_fix if set, max abs ratio to show; data>clim_fix are shown as clim_fix
  #' @param cexRow rowlabel size for heatmap, set to ~invisible by default
  #' @param cexCol collabel size for heatmap, set to aheatmap default by default
  #' @param annRow named lists of row annotations; can be single list
  #' @param annColors  named lists of annotation colors; names should match names in annRow and attribs
  #' @param labcoltype colnames to show ratiomat column names, colnums to show col #s
  #' @param labRowType  labels for rows; default NULL doesn't label
  #' @param setCol_v name of attrib to sort matrix columns by (will turn off column-clustering in heatmap output)
  #'      Defaults to NULL to keep original clustering (hierarchical)
  #' @param setRowv pass through to Rowv argument of aheatmap. Default to NULL, hierarchical clustering
  #' @param colOrder_v Vector containing all unique values of attribs[[setColv]] in desired output order. Defaults
  #'      to NULL so that if setColv is specified alone, output order will be alphabetical.
  #' @param inclLegend logical, should legend be included in plot.

  # imports
  require(NMF)
  require(RColorBrewer)
  colSpecs_lsv <- setColSpecs(ratiomat = ratiomat, attribs = attribs, setCol_v = setCol_v, colOrder_v = colOrder_v)
  ratiomat = colSpecs_lsv$ratiomat

  annColors = uni_qual_colors(attribs)

  # set row labels
  if(!is.null(labRowType)) {
    if (length(labRowType) == nrow(ratiomat)) {
      labRow = labRowType
      cexRow = 1
    } else {
      labRow = rownames(ratiomat)
    }
  } else { labRow = NULL }

  # calculate the number of colors in each half vector
  if(length(colorbrew)>1){ # color vector given
    halfbreak = length(colorbrew)/2
  } else if(any(grepl('^[^:]+:([0-9]+)$',colorbrew)) ){
    halfbreak = as.numeric(sub('^[^:]+:([0-9]+)$','\\1',colorbrew))/2
  } else { # aheatmap also takes a single colorspec and makes 2-color map
    halfbreak = 1
  }
  if(halfbreak>1){halfbreak=halfbreak-1} # make room for outer limits!

  # calculate the color limits 
  # if normmat is provided: clim.pct percentile of data, divided in half (~dist to med) gets color scale
  # if no normmat: clim.pct percentile of ratio data gets color scale
  #    values outside this range get max color
  # inner limit
  if( !is.null(normmat) ){
    c.lim = quantile( sapply(1:nrow(normmat),
                           function(x){ diff(range( normmat[x,], na.rm=T )) }),
                    probs=clim.pct, na.rm=T)/2
  } else { 
    c.lim = quantile( ratiomat, probs=clim.pct, na.rm=T )
  }
  # outer limits
  if( !is.null(clim_fix) ){
    ratiomat[ratiomat>clim_fix] = clim_fix
    ratiomat[ratiomat<=-clim_fix] = -clim_fix
  }
  c.lim0 = max(abs(ratiomat),na.rm=T)
  c.min0 = min(ratiomat,na.rm=T) # true up lower limit -- make this optional?
  
  # make vector of colorbreaks
  if(halfbreak>1) { # use inner and outer limits
    colorbreaks=c(-c.lim0, seq(from=-c.lim, to=c.lim, by=c.lim/halfbreak), c.lim0)
  } else { # use outer limits only
    colorbreaks=c(c.min0,(c.min0+c.lim0)/2, c.lim0)
  }
  # If quantiles are the same for multiple columns, will have non-unique breaks which will fail   
  colorbreaks <- unique(colorbreaks)

  # Plot
  ah_ls = aheatmap(ratiomat, cexRow=cexRow, cexCol = cexCol, 
           color=colorbrew, breaks=colorbreaks,
           annCol=colSpecs_lsv$attribs, labCol=colSpecs_lsv$labCol_v, labRow=labRow,
           main=plottitle, annRow=annRow, annColors=annColors,
           sub=subtitle, Colv=colSpecs_lsv$setCol_v, Rowv = setRowv, annLegend = inclLegend)

  return(ah_ls)
}

getCidFromHeatmap <-
function(aheatmapObj, numSubTrees, cutByRowLogic=T){
  #' Return a list of cluster IDs and elements of dendrogram in left-right order
  #' @description 
  #' Takes a heatmap object (see makeHeatmap), and splits rows (or columns) into a dendrogram with specified number of trees. Cluster IDs and elemenets returned.
  #' @param aheatmapObj object returned by the aheatmap() function from the NMF package
  #' @param numSubTrees numerical vector of length 1 specifying number of subclades/subtrees to split the dendrogram row (or column) into.
  #' @param cutByRowLogic boolean specifying to return cluster IDs and orders for rows (TRUE) or columns (FALSE)
  #' @examples
  #' @export
  # 
  # (old comments and partial example)
  # Return a list of length 2 containing 1) a named vector of cluster IDs when cutting the dendrogram by numSubTrees and 2) a vector of elements in the order as they appear in the dendrogram (left to right for columns and top to bottom for rows).
  # aheatmapObj: an object returned by the aheatmap() function from the NMF package
  # numSubTrees: a numerical argument specifying the number of subclades/subtrees to split the dendrogram row or column into.
  # cutByRowLogic: a boolean specifying whether you want to return cluster IDs and orders for rows (TRUE) or columns (FALSE).
  # returns a list of length 2 containing 1) a named vector of cluster IDs when cutting the dendrogram by numSubTrees and 2) a vector of elements in the order as they appear in the dendrogram (left to right for columns and top to bottom for rows).
  # examples
  # normmat = NULL
  # ratiomat = structure(c(-2.96...
  # attribs = list(ARPEcellStatus=c(rep("No",1),rep("ARPE",1),rep("No",1)))
  # plottitle = "CommonCellTypesInOurImportantTSSs"
  # clim.pct = .99
  # clim_fix = NULL
  # colorbrew = "-PiYG:64"
  # cexRow=0.01
  # labcoltype=c("colnames","colnums")
  # source("abundance_functions.R")
  # aheatmapObj = makeHeatmap (normmat=NULL, ratiomat=ratiomat, attribs=attribs, plottitle=plottitle, clim.pct=clim.pct, colorbrew=colorbrew, cexRow=cexRow, labcoltype=labcoltype)
  # numSubTrees = 4
  # cutByRowLogic = T
  # test_ls = getCidFromHeatmap(aheatmapObj, numSubTrees, cutByRowLogic)
  # test_ls$clusterIDByElement
  #  #1  2  3  4  5  6  7  8  9 10 
  #  #1  2  2  3  3  3  2  3  3  4 
  # test_ls$elementOrderByDendrogram
  # #[1]  3  2  7  4  9  8  6  5  1 10

  if(cutByRowLogic == TRUE){
    clusterIDByElement = cutree(as.hclust(aheatmapObj$Rowv), k = numSubTrees)
    elementOrderByDendrogram = rev(as.hclust(aheatmapObj$Rowv)$order)
    return_ls = list(clusterIDByElement, elementOrderByDendrogram)
    names(return_ls) = c("clusterIDByElement", "elementOrderByDendrogram")
  }else{
    clusterIDByElement = cutree(as.hclust(aheatmapObj$Colv), k = numSubTrees)
    elementOrderByDendrogram = as.hclust(aheatmapObj$Colv)$order
    return_ls = list(clusterIDByElement, elementOrderByDendrogram)
    names(return_ls) = c("clusterIDByElement", "elementOrderByDendrogram")
  }
  return(return_ls)
}

