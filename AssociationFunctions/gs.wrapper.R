#**********************************************************************************
# Revised version of gs.wrapper to calculate gene set enrichment for a list of signatures
# This should include the functionality of the V0 gs.wrapper, but is not backwards compatible
#
# ags.wrapper = wrapper for gs.pval to run calculations across multiple gene signatures 
#               and annotation sources
#               Changed the file write so that each element of setlist gets it's own file...
#                - this may not be quite right. 
#               Added "settings" information to returned value to indicate anno source details
# ags.settings = currently hard coded defaults for each of the annotation sources
#                Would like to see this developed further as we update our anno sources
#                Any specifications here are overridden by values in the gs.wrapper arguments
#                DEFAULT settings returned by gs.settings(infodump = TRUE)
#
# December 2016
#**********************************************************************************
#**********************************************************************************
ags.wrapper = function(setlist, backgroundset, annolist = NULL, anno_vals_ls = NULL,
                       ecut=0.05, ocut=5, anno.uni=TRUE,  
                       fileSettings=NULL, include.identifiers=FALSE,  return.OM=TRUE,
                       functiondir=NULL, resourcedir=NULL, verbose = FALSE)
{
  # Default annotation settings can be determined by typing gs.settings(infodump = TRUE)
  # 
  # Arguments:
  # setlist = list containing vectors of identifiers in backgroundset ID
  #    each list element is named with the display name for the gene set 
  #    Note that setlist MUST include a names() attribute.
  # backgroundset = data.table with fields ID, Entrez.ID and Anno.Symbol  
  #    Entrez.ID must be HUMAN IDs to use default annotation sources
  #    any other identifier types required by anno sources are in addnl columns
  #    required annotations by anno_sources --- gs.settings(infodump = TRUE)
  # annolist =  vector of annotation source variable names as strings
  #    NULL defaults to all default annotation sources
  # anno_val_ls = list of annotation settings by annolist element to override default values
  #    NULL defaults to all default values
  #    $anno_src = pertinent annotation source
  #      $gs.range = list of limits bounding the signatures that will be considered
  #          $min = minimum number of genes in an annotation set
  #          $max = maximum number of genes in an annotation set
  #      $maxe = maximum Evalue gene signature to return
  #      $idtype = anno_src gene identifier type; numeric preferred over symbol!
  #                note: currently "Entrez.ID" or "Anno.Symbol"
  #      $resourceset = file name for the default pathway annotation
  #                     Each file is an RData file. The filename should be the same
  #                     as that of the list varname inside the file.
  #      $annodir = overrides the "general" resource dir
  # ecut = general maximum E-value for annotations to return
  #    overridden for an anno_src by anno-val_ls[[anno_src]]$gs.range$maxe
  # ocut = general minimum number of genes in a pathway for reporting
  #    overridden for an anno_src by anno-val_ls[[anno_src]]$gs.range$mn
  # anno.uni = flag (default TRUE): 
  #    TRUE = use the annotation source's universe size; 
  #    FALSE = use the overlap of backgroundset & anno universes; more conservative
  #
  # fileSettings = structure of data used for saving data and settings
  #    NULL default - no output file written
  #    $directory = Directory to write files to. Defaults to the working directory
  #    $baseFilename = Base filename for any files; append file-specific info (including .csv)
  #                    Default = NULL, no files written
  # include.identifiers = TRUE returns the identifiers from the setlist element in filename
  #
  # return.OM = flag indicating whether occurrence matrix is to be returned
  #    occurrence matrix == logical matrix, genes x sets
  #    used for coldmap generation
  # functiondir = optional directory with functions to be loaded (default wd)
  #    currently loads PathwayAnalysis.R
  # resourcedir = general directory with annotation (.RData) files 
  #    Overridden for an anno_src by anno_val_ls[[anno_src]]$resourcedir
  #
  # verbose = FALSE toggles print statements that may be useful when there are problems...
  #
  # Returns a list by setlist element names:
  #   setlist element name 
  #     $GS_ls = List by setlist element
  #       $anno_src = list by annotation source
  #         data.table with gene signatures
  #     $OM_ls = List by setlist element
  #       $anno_src
  #         matrix of occurence matrix
  #     $AS_ls = list by anno_src indicating anno_vals used for processing and phyper release info
  #   
  
  if(is.null(functiondir)) functiondir=getwd()
  if(is.null(resourcedir)) resourcedir=getwd()
  if(is.null(fileSettings$directory)) fileSettings$directory=getwd()
  
  # imports
#  source(file.path(functiondir,'gs.pval.R'))
  source(file.path(functiondir,'PathwayAnalysis.R'))
  tmp=suppressPackageStartupMessages(require(data.table))

  GS_ls = list()
  OM_ls = list()
  AS_ls = list()
  
  # check arguments
  # setlist should be list of nonzero-length vectors
  if(!is.list(setlist) || sum(lengths(setlist)>0) != length(setlist) ){
    stop("setlist must be a list of nonzero-length vectors")
  }
  # get ready to save occurrence matrices if requested
  if( return.OM ){ OM_ls = NULL }

  # Set up the desired list of annotation sources
  if(is.null(annolist)) {
    annolist = names(gs.settings(infodump = TRUE))  
  }

  # Clear the data.table to be used for optional file output
  table2write = NULL

  # Process each set from setlist
  sigsets = names(setlist)
  for (mysig in sigsets) {
    # Process by anno
    for (myanno in annolist) {
      # Set up anno parameters
      params = gs.settings(anno_src = myanno, anno_vals = anno_vals_ls[[myanno]], resourcedir = resourcedir)
    
      message(sprintf("*** Anno Source %s:  any anno_vals specified? = %s", myanno, !is.null(anno_vals_ls[[myanno]])))
      # load annotation data
      # Identify the name of the loaded data base. Note that all include the anno_src name
      # Note: this should be cleaned up as the anno sources are cleaned up to a more deterministic method :-)
      pre_load_vars = NULL
      post_load_vars = NULL
      pre_load_vars = ls()
      load(file.path(params$annodir, params$resourceset))
      post_load_vars = ls()
      resource_varname = post_load_vars[which(!(post_load_vars %in% pre_load_vars))]
      if (verbose) { message("resource_varname = ", resource_varname)}
      myresource = get(resource_varname[grep(myanno, resource_varname)])
      
      # Convert Entrez IDs to numeric for efficient processing.
      if(params$idtype == "Entrez.ID") {
        myresource = lapply(myresource,function(x) as.numeric(x))
      }
      
      # Setup enrichment results data.table
      out.dt = data.table(rn="",E.Value="",P.Value="",Set.Size="",Overlap="",
                          Sig.Size="",Uni.Size="",Source="",Condx="",FoldEnrich="")
      if(include.identifiers){
        out.dt = cbind(out.dt,ID.list="",alt.ID.list="")
      }
      # Setup settings data.table
      settings.dt = data.table(AnnotationSource = "", AnnoSrc_IDs="",
                               min.n.pwayGenes="", max.n.pwayGenes="", MaxSignificance= "",
                               Species="", Version="",AnnoDirectory = "",
                               Stats.PackageVersion="", nOverlapGenes="")
    
      # Set up ID type(s) to use; include both the "native" ID for the anno_src
      #  and an alternative Symbol (human readable) ID; when the native ID for the anno_src
      #  is already Symbol, then use Entrez as the alternative
      alt.idtype = "Anno.Symbol" 
      if(params$idtype == "Anno.Symbol") {
        alt.idtype = "Entrez.ID"
      }

      # Set up universes - anno_src has a list of genes, the experiment has a list of genes
      # If anno.uni, the universe is the list of genes in the anno_src
      # If !anno.uni, universe is the intersection of genes in the backgroundset and anno_src
      #   !anno.uni is more conservative in the p-value department
      # Note that this may be different based on idtype, as there is not always a 1:1 mapping 
      # between Entrez and Symbol
      
      # measure the universes for the experimental set and annotation source
      anno_uni = unique(unlist(myresource, F, F))
      expt_uni = unique(backgroundset[, get(params$idtype)])
      uni_both = expt_uni[which(expt_uni %in% anno_uni)]
      message(sprintf("Annotation Source %s:\t %i pathways,\t%i unique genes, expt_uni = %i, uni_both = %i",
                      myanno, length(myresource), length(anno_uni), length(expt_uni), length(uni_both)))
      if( anno.uni ){
        uni.size=length(anno_uni)
      } else {
        uni.size = length(uni_both) 
      }
      message(sprintf("anno.uni = %s, Universe size for %s analysis: %i", anno.uni, myanno, uni.size))
    
      # annotate each set
      anno = gs.pval(gs = myresource, rnames=backgroundset[, get(params$idtype)], 
                     is.sig = backgroundset$ID %in% setlist[[mysig]], uni.size=uni.size,
                     min.size=params$gs.range['min'], max.size=params$gs.range['max'],
                     include.identifiers=include.identifiers, 
                     alt.rnames=backgroundset[,get(alt.idtype)],return.OM=return.OM)
      
      # Convert the GS to a data.table with signature data as appropriate
      myGS = anno$GS
      sigmask = (myGS[,"E.Value"] < params$maxe &
                 myGS[,"Overlap"] >= params$gs.range['min'])
      sig.n = sum(sigmask, na.rm = T)
      print(paste(mysig, myanno,"uni.size",attr(myGS,"Uni.Size"),"sig.n", sig.n))

      if(sig.n) {
        # convert slice of GS to data.frame to preserve structure if n==1
        out.dt = data.table(myGS, keep.rownames = T)[sigmask]
        names(out.dt)[names(out.dt) == "rn"] = "Set.Name"
        out.dt[, "Sig.Size" := attr(myGS, "Sig.Size")]
        out.dt[, "Uni.Size" := attr(myGS, "Uni.Size")]
        out.dt[, "Source" := myanno]
        out.dt[, "Condx" := mysig]
        out.dt[, "FoldEnrich" := (Overlap / Sig.Size)/(Set.Size/Uni.Size)]
        if(include.identifiers){
          out.dt[, "ID.list" := attr(myGS,"ID.list")[sigmask] ]
          out.dt[, "alt.ID.list" := attr(myGS,"alt.ID.list")[sigmask] ]
        }
      }
      
      # Filter OM to only include pathways included in the GS (out.dt) and genes in those pways
      outOM.dt = anno$OM[,colnames(anno$OM) %chin% out.dt$Set.Name]
      if (is.vector(outOM.dt)) {
        outOM.dt = outOM.dt[outOM.dt > 0]
      } else {
        outOM.dt = outOM.dt[rowSums(outOM.dt) > 0,]
      }

      GS_ls[[mysig]][[myanno]] = out.dt
      if (return.OM) {
        OM_ls[[mysig]][[myanno]] = outOM.dt
      }

      # Create a merged file for output if a file specified
      if(!is.null(fileSettings) & sig.n > 0) {
        if(is.null(table2write)) {
          table2write = out.dt
        } else {
          table2write = rbind(table2write, out.dt)
        }
      }
      
      # Need to clear out the resource references prior to load
      # This is essential for the GO set with three separate annos
      rm(list = resource_varname)
    
      # Update association source settings file
      AS_ls[[myanno]] = params
      AS_ls[[myanno]]["stats.pkgVersion"] = packageVersion("stats")
      AS_ls[[myanno]]["assoc_uni"] = ifelse(anno.uni, "Unique genes in annotation source",
                                            "Overlap unique genes in annotation source and background set")
    
    }  # end of myanno loop
  }  # end of mysig loop through each list of significant genes
  
  # Save enrichment results and settings
  if(!is.null(fileSettings)){
    # Write Enrichment 
    outfile = paste(fileSettings$baseFilename, mysig, "csv", sep = '.')
    outfilename = paste(fileSettings$directory, outfile, sep = '/')
    write.table(table2write, file=outfilename, row.names=F, quote=F, sep = ',')
    
    # Write Settings - TO BE IMPLEMENTED!!!
    # TAGL!!! Do this...
    # Create write_AS() and read_AS() functions to convert.
    # Need to convert to a data.table. Format specified in your notes. 
#    outfile = paste(fileSettings$baseFilename, assoc_settings, "csv", sep = '.')
#    outfilename = paste(fileSettings$directory, outfile, sep = '/')
#    write.table(settings2write, file=outfilename, row.names=F, quote=F, sep = ',')
  }
  
  if(return.OM) {
    return(list(GS = GS_ls, OM = OM_ls, Settings = AS_ls))
  } else {
    return(list(GS = GS_ls, Settings = AS_ls))
  }
}


#***********************************************************************************
# Assign default values as appropriate to included gs annotation sources
# Default values are defined in this function
# TBD means ToBeDeveloped
#
# Arguments:
#  anno_src = the name of the annotation source. Can be one of the defaults, or a new source
#  anno_vals = optional list of named objects; set values will be used; NULL revert to default
#              note: defaults differ by the annotation source
#    $gs.range = list of limits bounding the annotation signatures that will be considered
#        $min = minimum number of genes in an annotation set to consider
#        $max = maximum number of genes in an annotation set to consider
#    $maxe = maximum Evalue gene signature to return
#    $idtype = anno_src gene identifier type; numeric preferred over symbol!
#              note: currently "Entrez.ID" or "Anno.Symbol"
#    $resourceset = file name for the default pathway annotation
#                   Each file is an RData file. The filename should be the same
#                   as that of the list varname inside the file.
#    $resourcedir = overrides the "master" resource dir... Hook for future development
#    $species = TBD - species indicator; default is human
#    $version = TBD - anno_src release indicator; default is latest
#  resourcedir = optional general directory containing the annotations
#                defaults to the current working directory
#                overridden by anno_vals$resourcedir
#  ecut = general maximum adjusted evalue (adjusted pvalue) to report 
#         overridden by specific anno_vals$gs.range$maxe
#  ocut = general minimum number of genes in an anno pathway to merit reporting
#         overridden by specific anno_vals$gs.range$min
#         note that the default ocut value overrides this for some anno_src
#  infodump = if TRUE, returns annodefaults
#
# Returns a list of default/specified values
#  annoset = list of annotations to be used for the anno_src
#            all values listed in anno_vals above


gs.settings = function(anno_src, anno_vals = NULL, resourcedir=NULL, ecut=0.05, ocut=5,
                       infodump = FALSE)
{
  annoset = list()
  # Set the default general value for general resourcedir if not specified
  if (is.null(resourcedir)) {
    resourcedir = getwd()
  }

  # TBD:  Would like to see this replaced with some sort of dynamic table
  # Assign default values per anno source
  annodefaults = list(kegg = list(gs.range = list(min = ocut, max = Inf),
                                  maxe = ecut,
                                  idtype = 'Entrez.ID',
                                  resourceset = "hsa.kegg.gs.RData",
                                  species = "human",
                                  version = "bcorelatest",
                                  annodir = resourcedir),
                      msigdb5 = list(gs.range = list(min = 25, max = 2500),
                                     maxe = ecut,
                                     idtype = 'Entrez.ID',
                                     resourceset = "msigdb.v5.0.entrez.gmt.RData",
                                     species = "human",
                                     version = "bcorelatest",
                                     annodir = resourcedir),
                      GOBP = list(gs.range = list(min = 25, max = 2500),
                                  maxe = ecut,
                                  idtype = 'Entrez.ID',
                                  resourceset = "GO_processed.gmt.RData",
                                  species = "human",
                                  version = "bcorelatest",
                                  annodir = resourcedir),
                      GOMF = list(gs.range = list(min = 25, max = 2500),
                                  maxe = ecut,
                                  idtype = 'Entrez.ID',
                                  resourceset = "GO_processed.gmt.RData",
                                  species = "human",
                                  version = "bcorelatest",
                                  annodir = resourcedir),
                      GOCC = list(gs.range = list(min = 25, max = 2500),
                                  maxe = ecut,
                                  idtype = 'Entrez.ID',
                                  resourceset = "GO_processed.gmt.RData",
                                  species = "human",
                                  version = "bcorelatest",
                                  annodir = resourcedir),
                      pC = list(gs.range = list(min = 25, max = 2500),
                                maxe = ecut,
                                idtype = 'Anno.Symbol',
                                resourceset = "pathwayCommons.gmt.RData",
                                species = "human",
                                version = "bcorelatest",
                                annodir = resourcedir),
                      nci = list(gs.range = list(min = ocut, max = Inf),
                                 maxe = ecut,
                                 idtype = 'Anno.Symbol',
                                 resourceset = "NCI-Nature_Curated.gmt.RData",
                                 species = "human",
                                 version = "bcorelatest",
                                 annodir = resourcedir),
                      BioCarta = list(gs.range = list(min = ocut, max = Inf),
                                      maxe = ecut,
                                      idtype = 'Anno.Symbol',
                                      resourceset = "BioCarta.gmt.RData",
                                      species = "human",
                                      version = "bcorelatest",
                                      annodir = resourcedir),
                      Reactome = list(gs.range = list(min = ocut, max = Inf),
                                      maxe = ecut,
                                      idtype = 'Anno.Symbol',
                                      resourceset = "Reactome.gmt.RData",
                                      species = "human",
                                      version = "bcorelatest",
                                      annodir = resourcedir),
                      FF10 = list(gs.range = list(min = ocut, max = Inf),
                                     maxe = ecut,
                                     idtype = 'Entrez.ID',
                                     resourceset = "FF10_Entrez.RData",
                                     species = "human",
                                     version = "bcorelatest",
                                     annodir = resourcedir)
                     )

  # Return default settings for reference
  if (infodump) {
    return(annodefaults)
  }

  valuelist = names(annodefaults$kegg)
  for(valelem in valuelist) {
    if (valelem == "gs.range") {
      rangelist = names(annodefaults$kegg$gs.range)
      for(elem in rangelist) {
        annoset$gs.range[[elem]] = ifelse(is.null(anno_vals$gs.range[[elem]]),
                                      annodefaults[[anno_src]]$gs.range[[elem]],
                                      anno_vals$gs.range[[elem]])
      }
    } else {
        annoset[[valelem]] = ifelse(is.null(anno_vals[[valelem]]),
                                annodefaults[[anno_src]][[valelem]], anno_vals[[valelem]])
    }
  }
  return(annoset)
}
