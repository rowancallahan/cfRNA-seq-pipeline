#####################################################################################################
# Set of functions useful for running and plotting results of multiple testing correction (e.g., q-value)
#
# qcQ_wrapper     --- Plots p-value and q matrix for all factors of a linear regression 
#                     Optionally plots only selected factors
# qcQvalues       --- Plot p-value histogram, q-value auto-qc & MDS w q-cut
# qc.clusters     --- Plot heatmap and MDS for abundance data
#
# Author: Theresa Lusardi, Julja Burchard
# Started - July, 2016
#####################################################################################################



qcQ_wrapper <-
function(norm_x, p_mat, q_list, metadata, plotdata, plot2file = TRUE,
                       bonusMDS = FALSE, oneclass = NULL, qcut = 0.1, 
                       histbins=40, colorspec = c("#7b3294", "#008837"),
                       facSel = NULL, filesep="/") {

# Wrapper for qcQvalue function to verify q-value analysis terms
# Plots for each design factor (default) or for factors specified in facSel
#   1) Histogram of p values that were included in the design
#   2) qvalue's default plots, with full qvalue range c(0,1) plotted
#      Look for even descent to pi0 in the top left tuning plot
#      The slow/steep rise in q-values in remaining plots depends on resolving power of data
#   3) Optional MDS plot (bonusMDS=TRUE) at q-value cut (default 0.1, see below)
#
# Arguments
#   norm_x: abundance data input into regression (eg, regressMatrix()). nrow(norm_x)==length(p_mat)
#   p_mat: matrix of p-values for each factor. Returned by regressMatrix()
#   q_list: S3 object; list of p-value lists for each factor. Returned by regressMatrix() #Feedback supposed to be q-value lists?
#   metadata: list of factors; each list element is a factor (eg, sex = c("m", "f"...)) #Feedback perhaps as.factor(c("m","f"))?
#             Vector order same as norm_x column order
#             was attribs
#   plotdata is a list of info relevant to labeling and saving the plot #Feedback with the following elements:
#     plotdir:  plot destination directory
#     plotbase:  base filename for the plot. Suggest: bias reduction method
#     plottitle:  title for all plots
#     plotSubtitle:  optional subtitle; suggest lm_expression
#   plot2file if TRUE, otherwise to studio...
#   bonusMDS: if TRUE plots an MDS of data at qcut (see below)
#   oneclass: metadata factor for labeling bonusMDS data; #Feedback is this actually a factor, or is it a vector?
#             if NULL or not specified, defaults to the first element of metadata
#   qcut: either a number between 0 and 1 used as an upper qvalue limit for MDS,
#     OR a number >1, assumed to be the top n qvalues to use in MDS
#   histbins: # bins for p-value histogram; may want to change if there is suspicious behavior
#   colorspec: optional range of colors for MDS plot. Defaults to hotpink/green!
#   facSel: optional vector of experimental factors to plot (rather than all factors) #Feedback to plot for the MDS plot? What if you want to plot based on second element of metadata list?
#   #Feedback what about filesep?
# Returns... Error message or... #Feeback LOL! Thank you for the good humor!!
#
  
  # Define factors to plot
  if(is.null(facSel)) {
    facUse = colnames(p_mat)
  } else {
    # Confirm elements of facSel in the model
    facCheck = facSel %in% colnames(p_mat)
    if (sum(facCheck) == 0) {
      message(paste0("facSel value ", facSel, " not in model.\n"))
      message("Plotting defaults.")
      facUse = colnames(p_mat)
    } else
    if (sum(facCheck) != length(facSel)) {
      facBad = facSel[which(!facCheck)]
      message(paste0("facSel value ", facBad, " not in model.\n"))
      facUse = facSel[which(facCheck)]
      message(paste0("Plotting ", facUse, "\n"))
    }
  }

  # test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }

  # Establish the number of elements in the normalized data
  if(is.null(dim(norm_x))) {  # TRUE - norm_x is a vector
    norm_len = length(norm_x)
  } else { norm_len = nrow(norm_x) }

  ##### Plot data for each factor in facUse
  message("facUse length = ", length(facUse))
  orig_plotbase = plotdata$plotbase
  orig_plottitle = plotdata$plottitle
  for (factor in facUse) {
    skipFac = FALSE
    plotdata$plotbase = paste(orig_plotbase, factor, sep = "_")
    plotdata$plottitle = paste(orig_plottitle, factor, sep = ", ")
    # Check length of p-value list against norm_x
    if(nrow(p_mat) != norm_len) {
      message(sprintf("WARNING: Length of %s p_vals does not match norm_x",factor))
      skipFac = TRUE
    }
    # Check length of qvalues
    if(length(q_list[[factor]]$qvalues) != norm_len) {
      message(sprintf("WARNING: Length of %s q_list$qvalues does not match norm_x", factor))
      skipFac = TRUE
    }

    if (!skipFac) {
      message(sprintf("Factor %s to qcQvalues!", factor))
      qcQvalues(norm_x = norm_x, pvalue_v = p_mat[,factor], 
                obj_qvalue = q_list[[factor]], qcut = qcut, attribs = metadata,
                oneclass = oneclass, plotdata = plotdata, colorspec = colorspec,
                histbins=histbins, plot2file = plot2file, filesep=filesep,
                p_hist=TRUE, q_plots=TRUE, MDS_plot=bonusMDS)
    }
  }
  return("qcQ_wrapper done!")
}

qcQvalues <-
function (norm_x, pvalue_v, obj_qvalue, attribs, oneclass, plotdata, 
                      colorspec, qcut=1, histbins=40, plot2file = FALSE, 
                      p_hist=TRUE, q_plots=TRUE, MDS_plot=FALSE, 
                      filesep='/') {
# This is a check on the proper running of q-value analysis
# plots:
#  1) histogram of p-values input to q-value. 
#     these must NOT have a right peak or qvalue() will return spurious results
#  2) qvalue's default plots, with full qvalue range c(0,1) plotted
#     these should show an even descent to pi0 in the top left tuning plot,
#     and a slow or steep rise in q-values in remaining plots depending on 
#     resolving power of data
#  3) MDS plot restricted by q-value cut.
#  norm_x:  abundance data, vector or matrix. nrow(norm_x)==length(pvalue_v) 
#  pvalue_v: vector of p-values previously used as input to qvalue()
#  obj_qvalue: S3 object of qvalue class (a list!) returned by qvalue()
#  qcut: either a number between 0 and 1 used as an upper qvalue limit for MDS,
#     OR a number >1, assumed to be the top n qvalues to use in MDS
#  attribs:  list of sample classifications to be tracked in clustering
#     each list element contains a string vector with one label per sample
#  oneclass: string name of attribs element to be used in MDS plot
#  p_hist, q_plots, MDS_plot: flags to plot (TRUE) or skip (FALSE) plot types
#     MDS_plot now defaults to FALSE -- please use plotRatios() instead
#     after selecting a q-value cut based on these QC plots
#  plotdata is a list of info relevant to labeling and saving the plot
#     plotdir:  plot destination directory
#     plotbase:  base filename for the plot. Suggest: bias reduction method
#     plottitle:  title for all plots
#     plotSubtitle:  optional subtitle; suggest lm_expression
#  colorspec is a vector of color specifiers for colorRampPalette  
#  histbins: number of bins for p-value histogram
#  plot2file: if TRUE save plots to location defined in plotdata$plotdir
#  fac_sel:  NULL defaults to plotting all factors; can specify subset of factors in a vector
#  lm_expression is the linear model (eg, "y ~ factor1 + factor2"), included as a subtitle and filename value
#  p_hist, q_plots, MDS_plot FALSE allow suppression of a plot type
#  filesep: filepath separator

  # constants
  fudgefac = 2 # fold excess permitted in rows recovered by top n qcut when specified as a quantile of qvalues

  # imports
  require(qvalue)
  require(limma)

  # argument tests
  if( is.null(dim(norm_x)) ){ # not matrix
    norm_len = length(norm_x)
  } else { norm_len = nrow(norm_x) 
  }
  if( !exists("qvalues", obj_qvalue) ){
    stop("qvalue object does not contain qvalues element")
  }
  if( norm_len != length(pvalue_v) | norm_len != length(obj_qvalue$qvalues) ){
    stop("Data matrix has ",norm_len," rows, p-value vector has ",length(pvalue_v)," elements and qvalue object has ",length(obj_qvalue$qvalues)," q-values, but all must be the same size")
  }

  # test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }

  if( p_hist ){
    plotID = '8p'
    pD = 'p.value_histogram' 
    subTrue = is.vector(plotdata$plotSubtitle)
    plotDesc = ifelse(subTrue, paste(make.names(plotdata$plotSubtitle),pD, sep = '_'), pD)
    if(plot2file) {
    png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
        width=5,height=5.4,units="in",res=300)
    }
    hist(pvalue_v, nclass=histbins, main=plotdata$plottitle)
    if(is.vector(plotdata$plotSubtitle)) {
      mtext(plotdata$plotSubtitle)
    }
    if(plot2file) dev.off()
  }


  if( q_plots ){
    plotID = '8q'
    pD = 'q.value_QC.plot.array' 
    plotDesc = ifelse(subTrue, paste(make.names(plotdata$plotSubtitle),pD, sep = '_'), pD)
    if(plot2file) {
    png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
        width=7,height=5.4,units="in",res=300)
    }
    plot(obj_qvalue, rng=c(0,1), cex.axis=.6)
  
    if(plot2file) dev.off()
  }
  

  # make first-draft MDS plot if abundance data are a matrix 
  #  (multiple samples to compare)
  obj_MDS = NULL
  
  if( MDS_plot ){
    if( !is.null(dim(norm_x)) ){
      plotID = '6q0'
      plotDesc = 'MDS_q.value_QC' 
      if(plot2file) {
      png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
          width=5.4,height=5.4,units="in",res=300)
      }
  
      # data to plot
      if( !is.numeric(qcut) | qcut<0 | qcut>length(obj_qvalue$qvalues) ){
        stop(paste(qcut," must be a number > 0 and <= nrow(data)"))
      } else if( qcut <=1 ) { # assume this is a q-value on which to cut
        mymask = obj_qvalue$qvalues < qcut & !is.na(obj_qvalue$qvalues)
        titlestr= paste("qvalue <",signif(qcut,2) )
      } else { # assume this is a number of top genes by qvalue on which to cut
        q_quan = quantile(obj_qvalue$qvalues, probs=qcut/length(obj_qvalue$qvalues), na.rm=T)
        mymask = obj_qvalue$qvalues <= q_quan & !is.na(obj_qvalue$qvalues)
        titlestr= paste("top",sum(mymask),"q-values <=",signif(q_quan,3) )
        # test for lack of resolution in q-values; if so use p-values instead
        if( sum(mymask)>(qcut*fudgefac) ) { #identical qvalues at cut increase sum
          p_quan = quantile( pvalue_v, probs=qcut/length(pvalue_v), na.rm=T)  
          mymask = pvalue_v <= p_quan & !is.na(obj_qvalue$qvalues)
          titlestr= paste("top",sum(mymask),"p-values <=",signif(p_quan,3) )
        }
      }
  
      obj_MDS = makeMDSplot(normmat=norm_x[mymask,], attribs=attribs, oneclass=oneclass, colorspec=colorspec, plottitle=plotdata$plottitle, subtitle=titlestr, ngenes=sum(mymask))
  
      if(plot2file) dev.off()
    }
  }

  invisible( list(rowmask=mymask,obj_MDS=obj_MDS) )

}
qc.clusters <-
function (normmat, rawmat, attribs, oneclass, plotdata,
                        colorspec, center = 'norm', clim.pct, mask = NULL, 
                        plot2file = FALSE, filesep='/') {
# This function considers only a single clustering variable... Will need to be gneralized
# Data will be centered around the average of all data in the normalized data set
# As a QC plot, want to filter out non-changing data, otherwise the alg gets slow
# Uses the attributes to cluster the data; could make this a longer list...
#  normmat:  data matrix 
#  rawmat:  raw data matrix (optionally used to center normmat)
#  attribs:  list of sample classifications to be tracked in clustering
  # each list element contains a string vector with one label per sample
#  oneclass: string name of attribs element to be used in MDS plot
#  plotdata is a list of info relevant to labeling and saving the plot
#    plotdir:  plot destination directory
#    plotbase:  base filename for the plot
#    plottitle:  title for all plots
#  colorspec is a vector of color specifiers for colorRampPalette  
#  center:  center data on 'norm' for normalized average
#           none indicates no-centering needed for zero-centered data
#           any other value will trigger centering on the average of rawmat
#  clim.pct:  0:1 - fraction of data to limit max color
#  mask:  NULL default masks rows with SD > SD of entire data set; 
#         Numeric vector of length 1 masks rows with SD > value
#         Vector of length nrows(normmat) is a custom mask   ***** TBImplemented
#  plot2file:  if true, plot to file as specified in plotdata

  # imports
  require(NMF)
  require(limma)

# test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }


  plotID = 6
  plotDesc = 'Heatmap' 
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5,height=7,units="in",res=300)
  }

  # Center data
  if(center == 'norm') {
    ratiomat = sweep(normmat, 1, rowMeans(normmat,na.rm=T), '-')
    message(sprintf('centered on norm, max=%1.3f, min=%1.3f', max(ratiomat,na.rm=T), min(ratiomat,na.rm=T)))
  } else if (center == 'none') {
    ratiomat = normmat
  } else ratiomat = sweep(normmat, 1, rowMeans(rawmat,na.rm=T), '-')
  
  # Determine a lower threshold for data of interest
  # If there is too much non-changing data, the cluster calculation will stall
  # As this is a QC plot, limit to genes more variable than average
  mySD = sd(ratiomat, na.rm = T)
  mymask = sapply(1:nrow(ratiomat),function(x){sd(ratiomat[x,],na.rm=T)}) > mySD
  message(sprintf('mySD = %1.3f, mymask rows = %s, !mymask rows = %s',
                   mySD, sum(mymask), sum(!mymask)))

  ah_ls = makeHeatmap(normmat=normmat[mymask,], ratiomat=ratiomat[mymask,], attribs=attribs, plottitle = plotdata$plottitle, clim.pct=clim.pct)

  if(plot2file) dev.off()


  plotID = 7
  plotDesc = 'MDS' 
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5.4,height=5.4,units="in",res=300)

  obj_MDS = makeMDSplot(normmat=normmat, attribs=attribs, oneclass=oneclass, colorspec=colorspec, plottitle=plotdata$plottitle, ngenes=sum(mymask))

  if(plot2file) dev.off()

  invisible(list(ah_ls=ah_ls, obj_MDS=obj_MDS))
  }
}
