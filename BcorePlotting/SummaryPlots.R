#####################################################################################################
# Functions involved in generating many summary plots for rapid QC
#
# summary.plots    --- Creates box, scatter, density, SD/meanm, and (optionally) MA plots
# scatterplot      --- X-Y Scatter plot
# plot_raw_v_norm  --- Creates raw vs. norm plot; needs improvements
# maplot           --- A plain vanilla MA plot with loess fit line
# makePairwisePlotMatrix --- A wrapper for GGally::ggpairs for a matrix of MA plots. Incomplete and not incorporated into summary.plots
#
# Author: Theresa Lusardi, Julja Burchard, and Mark Fisher
# Started - July, 2016
#####################################################################################################

##Please be on the lookout for #Attn. and #Feedback tags in lines that warrant them.

summary.plots <-
    function (rawmat, normmat, mynorm, samp.labels, samp.classes, colorspec, plotdata, plot2file = FALSE, histbins = 40, expand.2D = 5, SDrange=7, filesep="/", plotIDOffset = 0, verbose = FALSE, BaseSample = NULL, yrange = NULL, MAplotOpt=FALSE,
              whichPlots_v = c("box", "scatter", "density", "spread")) {
    #' A wrapper for several plotting functions that help summarize abundance data #Feedback since rawmat in this case is on log2 scale, will not yet accommodate methylation data
    #' @description summary.plots currently wraps a boxplot (1), scatterplot (2), density plot (3), sd vs. intensity plot (4), and MA plot (5). The MA plot behavior plots all sample columns in a matrix vs. one BaseSample.
    #' @param rawmat a matrix of raw data. Should have minimal background addition and be scaled according to normalized data. In other words, rawmat should be the alograw slot of a normMatrix() call.
    #' @param normmat a matrix of experimental data in columns (with headers!)
    #' @param mynorm a label for the current normalization method #Feedback not sure what this data type is
    #' @param samp.labels a vector of brief, pithy display labels for each sample
    #' @param samp.classes a vector of tags for each sample, used to pick colors
    #'  a unique tag per sample colors by sample (e.g., c("01", "02", "03", "04"))
    #'  a shared tag by experimental groupings colors by experimental group (e.g. c("male", "female", "female", "male"))
    #' @param colorspec a vector of color specifiers for colorRampPalette  
    #' @param plotdata a list of info relevant to labeling and saving the plot. Consists of three elements:
    #'  plotdir:  plot destination directory
    #'  plotbase:  base filename for the plot
    #'  plottitle:  title for all plots
    #' @param histbins default number of bins for density plot
    #' @param expand.2D a numeric multiplier of histbins for 2D histogram
    #' @param SDrange maximum SD to plot in spread plots
    #' @param filesep a string to designate file separators #Feedback not sure
    #' @param plotIDOffset a number the specifies offset value for the plot ID. Default value is 0.
    #' @param BaseSample a string of the column/sample name that will be used as the sample to which other samples will be compared for MA plots.
    #' @param yrange a numerical vector of length 2 that sets the y-axis range. Default is NULL.
    #' @param MAplotOpt a boolean that specifies whether the user wants to generate MA plots in this call
    #' @return if verbose is TRUE, a list of ___
    #' @examples 
    #' Examples go here
    #' @export
    
    # imports
    # browser()
    require(data.table)
    
    # test plotdir for filesep; add if absent
    if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
      plotdata$plotdir = paste0(plotdata$plotdir,filesep)
    }
    
    # Trouble shooting messages if verbose == TRUE
    if(verbose) {cat("Create summary plots for", plotdata$plottitle, '\n')}
    if(verbose) {cat("Plot to file", ifelse(plot2file, "enabled", "disabled"), '\n')}
    
    # set up plotting colors
    u.samp.classes = unique(samp.classes)
    colmap = colorRampPalette(colorspec)
    colvec = colmap(length(u.samp.classes))
    # map colors back to samples
    colvec = colvec[as.numeric(as.factor(samp.classes))]
    # R doesn't resort on unique: saves order of occurrence
    u.col.classes = unique(colvec) 
    
    # Boxplot of all data
    if ("box" %in% whichPlots_v){
    plotID = 1 + plotIDOffset
    plotDesc = 'boxplot'
    if(verbose) {message(sprintf("Plotting %s, ID = %i", plotDesc, plotID))}
    if(plot2file) {
      png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
          width=5,height=5.4,units="in",res=144)
    }
    boxplot(normmat, main = plotdata$plottitle, border=colvec, las = 2)
    if(plot2file) dev.off()
    } # boxplot
    
    # Scatterplot of all data
    if ("scatter" %in% whichPlots_v){
    plotID = 2 + plotIDOffset
    plotDesc = 'scatter'
    if(verbose) {message(sprintf("Plotting %s, ID = %i", plotDesc, plotID))}
    if(plot2file) {
      png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
          width=5,height=5.4,units="in",res=144)
    }
    yl = range(normmat,na.rm=T); xl = range(rawmat,na.rm=T)
    plot(x=0, y=0, type="n", xlim=xl, ylim=yl, pch='.', xlab='', ylab = '')
    title(xlab='log2(raw + background)',
          ylab = paste(mynorm, 'normalization'),
          main = plotdata$plottitle)
    for(j in 1:ncol(rawmat) ){
      points(rawmat[,j], normmat[,j], pch='.', col=colvec[j])
    }
    legend(x="topleft", legend=u.samp.classes, col=u.col.classes, pch=16, cex=.5)
    if(plot2file) dev.off()
    } # scatterplot
  
  # Density plot
  if ("density" %in% whichPlots_v){
    plotID = 3 + plotIDOffset
    plotDesc = 'density'
    if(verbose) {message(sprintf("Plotting %s, ID = %i", plotDesc, plotID))}
    if(plot2file) {
      png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
          width=5,height=5.4,units="in",res=144)
    }

    dl = NULL; yl = c(0,0); xl = range(normmat, na.rm=T)
    for( j in 1:ncol(normmat) ) {
      dl = c(dl,list(density(normmat[,j],from=xl[1],to=xl[2],n=histbins,na.rm=T)))
      dl[[j]]$y = dl[[j]]$y/sum(dl[[j]]$y)
      yl = range(yl,dl[[j]]$y)
    }
    plot(x=0,y=0,type="n",xlim=xl,ylim=yl,
         xlab=paste(mynorm,' normalized intensity'),
         ylab='Fraction of features (smoothed)',
         main = plotdata$plottitle)
    # browser()
    for(j in 1:ncol(normmat) ){
      lines(dl[[j]]$x,dl[[j]]$y,lty=(j%%ncol(normmat)) + 1, col=colvec[j])
    }
    legend(x="topright", legend=u.samp.classes, col=u.col.classes, lty=1, lwd=2, cex=.6)
    
    if(plot2file) dev.off()
  } # density
    
    # SD plot
    # browser()
    if ("spread" %in% whichPlots_v){
    plotID = 4 + plotIDOffset
    plotDesc = 'spread'
    if(verbose) {message(sprintf("Plotting %s, ID = %i", plotDesc, plotID))}
    if(plot2file) {
      png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
          width=5,height=5.4,units="in",res=144)
    }
    
    # calculate group-wise SDs per gene
    sdmat = NULL; int_mat = NULL
    # find rows informative on SD for all classes
    row_mk = logical(nrow(normmat)); row_mk[] = TRUE
    # Create a mask indicating which rows have at least 1 measurement per samp.class
    for( myclass in u.samp.classes ){
      if(sum(samp.classes==myclass)==1){ #Catches those cases where only one sample is in a class
        warning(paste0(myclass, " is represented by only one sample. Be aware that this may be problematic for downstream analyses."))
        row_mk = row_mk & !(is.na(normmat[,samp.classes==myclass]))
      }else{
        row_mk = row_mk & rowSums(!is.na(normmat[,samp.classes==myclass]))>1
        if(verbose) {message(sprintf("Masked rows for %s = %i", myclass, sum(row_mk)))}
      }
    }
    # calculate mean and SDs for intensities by samp.classes
    for( myclass in u.samp.classes ){
      if(sum(samp.classes==myclass)==1){ #Catches those cases where only one sample is in a class. Can't take rowMeans of one column.
        int_mat = cbind(int_mat, normmat[row_mk,samp.classes==myclass])
      }else{
        int_mat = cbind(int_mat, rowMeans(normmat[row_mk,samp.classes==myclass],na.rm=T))
      }
      if(sum(samp.classes==myclass)>1){
        sdmat = cbind(sdmat, sapply(which(row_mk),
                                    function(x){ sd( normmat[x,samp.classes==myclass,drop=F] ,na.rm=T) }) )
      } else {
        sdmat = cbind(sdmat, rep(0,nrow(normmat)) )
      }
    }
    colnames(sdmat) = u.samp.classes; colnames(int_mat) = u.samp.classes #Feedback my u.samp.classes is length 4 while dim sdmat is 0,5.?
    rownames(sdmat) = rownames(normmat)[row_mk]
    
    # colors
    # bright version of each color: increase to cmult% of distance to 100%
    cmult = 0.8
    clower = rgb(t( cmult*(255 - col2rgb(colvec)) + col2rgb(colvec) )/255 )
    # build color list
    collist = NULL; ncols = 64
    for(j in 1:ncol(sdmat) ){
      if(verbose) {message(sprintf("Class %s, n(sd=0) %i", colnames(sdmat)[j], sum(sdmat[, j]==0)))}
      # color axis for each group, lighter at low end & mostly transparent 
      col1 = unique(clower[samp.classes==colnames(sdmat)[j] ])[1]
      col2 = unique(colvec[samp.classes==colnames(sdmat)[j] ])[1]
      collist = c(collist, list(colorRampPalette( colors=c( adjustcolor(col1, alpha.f=.85), adjustcolor( col2, alpha.f=0.85) ), alpha=T)(ncols) ) )
    }

    # 2D histogram
    nbins = histbins*expand.2D
    # Set the y (SD) and x (mean) bin sizes based on the range of the masked abundance data
    # bin[n:(n+1)] defines the range of values included in the bin
    # plt[n] defines the center point for plotting
    yl = range(sdmat, na.rm=T)
    ybin = seq(from=yl[1],to=yl[2],length.out=nbins+1)
    yplt = filter(ybin,filter=c(.5,.5),sides=1)[2:(nbins+1)]
    xl = range(normmat[row_mk,], na.rm=T)
    xbin = seq(from=xl[1],to=xl[2],length.out=nbins+1) 
    xplt = filter(xbin,filter=c(.5,.5),sides=1)[2:(nbins+1)]
    cmax = log2(nrow(sdmat)/nbins)
    # calculate x, y, color
    freq = NULL;
    for(i in 1:ncol( sdmat ) ){
      # frequency of abundance vs SD as binned
      freq[[i]] = as.data.table(table( findInterval(int_mat[,i],xbin,all.inside=T), findInterval(sdmat[,i],ybin,all.inside=T) ))
      freq[[i]] = freq[[i]][,V1:=as.numeric(V1)]
      freq[[i]] = freq[[i]][,V2:=as.numeric(V2)]
      freq[[i]] = freq[[i]][,N:=log2(N+1)]
      # transform to color numbers
      freq[[i]] = freq[[i]][,Nn := round((ncols-1)*N/cmax)+1]
      freq[[i]][Nn>ncols,Nn := ncols]
      # trim empty bins
      freq[[i]] = freq[[i]][N>0,]
    }
    
    # make single x, y, color vecs -- alpha seems to have limits
    myx = NULL; myy = NULL; myc = NULL
    for(i in 1:ncol( sdmat ) ){
      myx = c(myx,freq[[i]]$V1)
      myy = c(myy,freq[[i]]$V2)
      myc = c(myc,collist[[i]][freq[[i]]$Nn])
    }
    randx = sample(length(myx))
    # adjust ylim if indicated
    if( exists("SDrange") & !is.null(SDrange) & SDrange>yl[1] ){
      yl[2] = SDrange
    }
    # plot
    plot(x=0,y=0,type="n",xlim=xl,ylim=yl,
         xlab=paste(mynorm,' normalized intensity'),
         ylab='Standard deviation',
         main = plotdata$plottitle)
    if(verbose) {
      message(sprintf("lengths: randx %i, myx %i, myy %i, myc %i",
                      length(randx), length(myx), length(myy), length(myc)))
      message(sprintf("lengths:  xplt %i, yplt %i",
                      length(xplt), length(yplt)))
      message(sprintf("lengths:  my[rand] x %i, y %i, c %i",
                      length(myx[randx]), length(myy[randx]), length(myc[randx])))
      message(sprintf("lengths:  plt[my[rand]] x %i, y %i, c %i",
                      length(xplt[myx[randx]]), length(yplt[myy[randx]]), length(myc[randx])))
      returnVal = list(randx=randx, myx=myx, myx.randx=myx[randx], xplt=xplt, xplt.myx.randx=xplt[myx[randx]], myy=myy)
    }
    
    points(xplt[myx[randx]],yplt[myy[randx]],col=myc[randx],pch=15,cex=.6)
    
    # add loess fits as in Cope et al. Bioinformatics 20:323-331 (2004), Fig.2
    for(i in 1:ncol(sdmat) ){
      lfit = lowess( y = sdmat[,i], x = int_mat[,i], f=0.2, delta=(1/histbins)*diff(range(int_mat[,i])) )
      # plot heavy lines a bit darker than regular colors
      # lowess returns list elements x and y in plotting order
      lines( lfit$x, lfit$y, col=adjustcolor(u.col.classes[i],red.f=.75,green.f=.75,blue.f=.75), lwd=3 )
    }
    
    # add legend
    legend(x="topright", legend=u.samp.classes, col=u.col.classes, pch=16, cex=.6)
    
    if(plot2file) dev.off()
    } # spread
                  
    ##Add MA plots and handling potential input of BaseSample
    if(MAplotOpt){
      # if(verbose) {message(sprintf("Plotting %s, ID = %i", plotDesc, plotID))} #Feedback #Attn. buggy?
      
      if(is.null(BaseSample)){
        BaseSample = colnames(rawmat)[1] #Default is the first sample
        message("No base sample to compare to all others provided. Defaulting to first sample.")
      }
      colnamesNotBaseSample = setdiff(colnames(rawmat), c(BaseSample))
      for (i in colnamesNotBaseSample){
        v1Inst = normmat[,grep(BaseSample, colnames(normmat))]
        v2Inst = normmat[,grep(i, colnames(normmat))]
        v12namesInst = list(v1 = BaseSample, v2 = i)
        plotdataInst = plotdata
        yrangeInst = yrange
        plot2fileInst = plot2file
        plotIDOffsetInst = plotIDOffset + 5
        # browser()
        message(paste0("Plotting ",BaseSample," vs. ",i)) #plotdataInst,yrangeInst,plot2fileInst,plotIDOffsetInst))
        maplot(v1=v1Inst, v2=v2Inst, v12names = v12namesInst, plotdata = plotdataInst, yrange = yrangeInst, plot2file = plot2fileInst, plotIDOffset = plotIDOffsetInst)
      }
    }
    # if(plot2file) {
    #   dev.off()
    # }
    if(verbose) {message("Summary plots completed")}
    if(!exists("returnVal")){
      returnVal = NULL
    }
    return(returnVal)
  } #END summary.plots

scatterplot <-
  function (normmat, attribs, plotdata, plot2file = FALSE, plotIDOffset = 0) {
    #' Display indivual replicates vs. average of experimental group in a scatterplot
    #' @description 
    #' Uses matrix of experimental observations from multiple testing groups. 
    #' Compares individual observations to corresponding group average.
    #' Creates one plot per column, and an average vs average
    #' @param normmat matrix of exerimental data. Rows = observations, columns = samples/experiments (requires named rows and columns)
    #' @param attribs vector of experimental categories of same length as ncol(normmat)
    #' @param plotdata list of vectors relevant to labeling and saving plot
    #'    plotdir:  plot destination directory
    #'    plotbase:  base filename for the plot
    #'    plottitle:  title for all plots
    #' @export
    
    plotID = ifelse(plotIDOffset == 0, '2a', paste0(plotIDOffset + 2, "a"))
    plotlims = c(floor(min(normmat, na.rm=TRUE)), ceiling(max(normmat, na.rm=TRUE)))
    expts = unique(attribs)
    n.expts = length(expts)
    
    for (expt in expts)  {
      avg = rowMeans(normmat[, attribs == expt])
      for (i in 1:ncol(normmat))  {
        
        if (attribs[i] == expt) {  # Only plot the data in the experimental group
          if(plot2file) {
            plotDesc = sprintf('Scatter.average_%s.vs.replicate_%s', expt, colnames(normmat)[i]) 
            png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID,
                                   plotdata$plotbase, plotDesc), width=5,height=5.4,units="in",res=144)
          }  # end of plotfile setup
          
          plot(avg, normmat[,i], pch='.', xlab='', ylab = '', xlim = plotlims, ylim = plotlims)
          avgranges = fivenum(avg)
          repranges = fivenum(normmat[, i])
          title(xlab = sprintf('Average %s', expt),
                ylab = sprintf('Replicate %s', colnames(normmat)[i]),
                main = plotdata$plottitle)
          mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                        repranges[1], repranges[2], repranges[3], repranges[4], repranges[5]),
                side = 2, line = 2, cex = 0.75)
          mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                        avgranges[1], avgranges[2], avgranges[3], avgranges[4], avgranges[5]),
                side = 1, line = 2, cex = 0.75)
          
          if(plot2file) dev.off()
        }  # End of plot replicates test
      }  # end of i loop through all columns
    }  # end of expt loop through experimental groups
    
    # Now plot averages against each other
    for (i in 1:(n.expts-1)) {
      for (j in (i+1):n.expts) {
        avgx = rowMeans(normmat[, attribs == expts[i] ])
        avgy = rowMeans(normmat[, attribs == expts[j] ])
        xrange = fivenum(avgx)
        yrange = fivenum(avgy)
        
        if(plot2file) {
          plotDesc = sprintf('Scatter.average_%s.vs.average%s', expts[i], expts[j])
          png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID,
                                 plotdata$plotbase, plotDesc),
              width=5,height=5.4,units="in",res=144)
        }  # end of plotfile setup
        
        plot(avgx, avgy, pch='.', xlab='', ylab='', xlim = plotlims, ylim=plotlims)
        title(xlab = sprintf('Average %s', expts[i]),
              ylab = sprintf('Average %s', expts[j]),
              main = plotdata$plottitle)
        mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                      xrange[1], xrange[2], xrange[3], xrange[4], xrange[5]),
              side = 1, line = 2, cex = 0.75)
        mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                      yrange[1], yrange[2], yrange[3], yrange[4], yrange[5]),
              side = 2, line = 2, cex = 0.75)
        if(plot2file) dev.off()
      }  # end of j loop
    }  # end of i loop
  }  # end of scatterplot()


plot_raw_v_norm <-
  function(linear_scale_raw_mat_with_rownames, str_of_raw_name=NULL, log2_scale_norm_mat_with_rownames, str_of_norm_name=NULL, xrng_vec=range(range(linear_scale_raw_mat_with_rownames, na.rm=T), range(2^log2_scale_norm_mat_with_rownames, na.rm=T)), yrng_vec=log2(xrng_vec), save_as_png=F, string_to_lead_file_name_with=NULL, color_by_vec_of_str=F, color_vec_of_strings=NULL){
    if(save_as_png==TRUE){
      png(filename=paste0(string_to_lead_file_name_with,"_raw_v_norm.png"),width=5,height=5.4,units="in",res=600)
      print(plot(as.matrix(linear_scale_raw_mat_with_rownames), as.matrix(log2_scale_norm_mat_with_rownames),log='x',pch='.',xlim=xrng_vec,ylim=yrng_vec, main=paste0("Raw vs. Norm: ", str_of_raw_name," vs. ", str_of_norm_name), xlab="Non-normalized", ylab="Normalized"))
      if(color_by_vec_of_str==TRUE & !is.null(color_vec_of_strings)){
        cool_cols = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(color_vec_of_strings)))
        for(j in 1:length(unique(color_vec_of_strings))){
          current_str = unique(color_vec_of_strings)[j]
          cols_of_current_str = which(color_vec_of_strings==current_str)
          for (k in cols_of_current_str){
            print(points(as.matrix(linear_scale_raw_mat_with_rownames)[,k], as.matrix(log2_scale_norm_mat_with_rownames)[,k], col=cool_cols[j], pch="."))
          }
        }
      }
      dev.off()
    }
    print(plot(as.matrix(linear_scale_raw_mat_with_rownames), as.matrix(log2_scale_norm_mat_with_rownames),log='x',pch='.',xlim=xrng_vec,ylim=yrng_vec, main=paste0("Raw vs. Norm: ", str_of_raw_name," vs. ", str_of_norm_name), xlab="Non-normalized", ylab="Normalized"))
    if(color_by_vec_of_str==TRUE & !is.null(color_vec_of_strings)){
      cool_cols = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(color_vec_of_strings)))
      for(j in 1:length(unique(color_vec_of_strings))){
        current_str = unique(color_vec_of_strings)[j]
        cols_of_current_str = which(color_vec_of_strings==current_str)
        for (k in cols_of_current_str){
          print(points(as.matrix(linear_scale_raw_mat_with_rownames)[,k], as.matrix(log2_scale_norm_mat_with_rownames)[,k], col=cool_cols[j], pch="."))
        }
      }
    }
  }
maplot <-
  function (v1, v2, v12names=list(v1="Sample 1", v2="Sample 2"), plotdata, yrange = NULL, plot2file = FALSE, plotIDOffset = 0) {
    #' Generates Bland-Altman plot: log ratio vs. average of two samples
    #' @description
    #' Note:  fivenumber or quantile vector will be included in the y-axis label, so it will be clear if the plot has been limited. 
    #' ... #Attn. #Feedback flesh out description more
    #' @param v1 a vector of abundance data for a sample of interest
    #' @param v2 a second vector of same length as v1.
    #' @param v12names a list with two elements: v1 and v2, which are strings describing vectors v1 and v2. These strings will appear in the plot. Defaults to calling them "Sample 1" and "Sample 2".
    #' @param plotdata a list with three elements of info relevant to labeling and saving the plot: 1. plotdir-plot destination directory. 2. plotbase-base filename for the plot. 3. plottitle-title for all plots
    #  plotdata is a list of info relevant to labeling and saving the plot
    #    plotdir:  plot destination directory
    #    plotbase:  base filename for the plot
    #    plottitle:  title for all plots
    #' @param yrange a numerical vector of length 2 that sets the y-axis range. Default is NULL.
    #' @param plot2file a boolean specifying whether the plot is to be saved to file according to the file-naming specification laid out in plotdata. Default is FALSE.
    #' @param plotIDOffset a number the specifies offset value for the plot ID. Default value is 0.
    #' @return #Attn. #Feedback
    #' @example 
    #' norm_mat = summ_ls$loess
    #' v1Inst = norm_mat[,1]
    #' v2Inst = norm_mat[,2]
    #' v12namesInst = list(v1=gsub('(^.+?)_.*','\\1', colnames(norm_mat)[1]),v2=gsub('(^.+?)_.*','\\1', colnames(norm_mat)[2]))
    #' plotdataInst = list(plotdir = "~/Desktop/", plotbase ="MA", plottitle = "MATestPlot")
    #' yrangeInst = NULL
    #' plot2fileInst = TRUE
    #' plotIDOffsetInst = 0
    #' maplot(v1=v1Inst, v2=v2Inst, v12names = v12namesInst, plotdata = plotdataInst, yrange = yrangeInst, plot2file = plot2fileInst, plotIDOffset = plotIDOffsetInst)
    #' @export
    
    # browser()
    plotID = plotIDOffset
    plotDesc = sprintf('MA_%s.vs.%s', v12names$v1, v12names$v2) 
    if(plot2file) {
      png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
          width=5,height=5.4,units="in",res=144)
    }
    
    M = v1 - v2
    A = 0.5*(v1 + v2)
    ranges = fivenum(M)
    if (!is.null(yrange)) {
      ylimits = c(-max(abs(yrange), na.rm = TRUE), max(abs(yrange), na.rm = TRUE))
      plot(A, M, pch = '.', ylim = ylimits)
    } else plot(A, M, pch = '.')
    title(xlab = 'A', ylab = 'M', main = plotdata$plottitle)
    mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                  ranges[1], ranges[2], ranges[3], ranges[4], ranges[5]),
          side = 2, line = 2, cex = 0.75)
    mtext (sprintf('%s\n%s', v12names$v1, v12names$v2),
           side = 1, line = 3, adj = 0, padj = 0, cex = .75)
    # browser()
    lfit = loess(M ~ A)  # Calculate a loess fit of the MA plot to show variance
    MC = predict (lfit, A)
    abline(h = 0, v = 0, col = 'grey60')  # Add an x axis at 0
    lines(A[order(A)], MC[order(A)], col = 'red')
    
    if(plot2file) dev.off()
  }

makePairwisePlotMatrix = function(norm_mat, plotdata, plot2file = FALSE, reso=600, ...){ #, , plotIDOffset = 0, ...){
  # Attn. Feedback Needs to be completed: 1) arguments need to be set up so that they dovetail with summary.plots, 2) I didn't have time to look into how to pass pch="."-like arguments into ggplot2 plots., 3) Needs to be entered into summary.plots above
  #' Wrapper for ggpairs, which generates a matrix of plots with a given data set
  #' @description 
  #' This is intended to show individual replicates vs. average of experimental group
  #' It will create one plot per column, and an average vs average
  #' @param ... addition parameters that ggpairs() can take
  #' normmat is a matrix of experimental data in columns (with headers!)
  #' attribs is a vector of experimental categories (one entry per column in normmat)
  #' plotdata is a list of info relevant to labeling and saving the plot
  #'  plotdir:  plot destination directory
  #'  plotbase:  base filename for the plot
  #'  plottitle:  title for all plots
  #'  @return 
  #'  @example 
  #'  normInst_mat = summ_ls$loess
  #'  plotdataInst = list(plotdir = "~/Desktop/", plotbase="BaseTest", plottitle = "titleTest")
  #'  plot2fileInst = TRUE
  #'  makePairwisePlotMatrix(norm_mat = normInst_mat, plotdata = plotdataInst, plot2file=plot2fileInst, pch=".")
  library(GGally)
  library(ggplot2)
  
  ##Make sure column names are kosher##
  colnames(norm_mat) = make.names(colnames(norm_mat))
  # browser()
  
  if (plot2file){
    tiff(filename = sprintf('%s%i_%s_%s.tiff', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc), width=5,height=5.4,units="in",res=reso)
    #Attn. sprintf('%s%i_%s_%s.tiff', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc)
    dev.off()
  }
  print(GGally::ggpairs(data=as.data.frame(norm_mat),title=plotdata$plottitle,columnLabels = gsub('(^.+?)_.*','\\1', colnames(norm_mat)), mapping = aes(size=0.2))) #+ geom_point(colour = "red", size = 1)) #,pch=pch
  if (plot2file){
    dev.off()
  }
  
}
