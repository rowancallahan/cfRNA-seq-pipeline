#*******************************************************************************************
# Functions to summarize AssociationFunction results
#
# Implemented
# clusterGO = clusters GO hierarchy by direct parent/child relationships
# coldClusters    Ranks pathway clusters from coldMap
# makeColdMap = Creates a coldmap from an occurrence matrix (OM)
# mergeOM = merges occurrence matrices. Could also be used to trim OM
#
# To be implemented...
# gsHeatMap = Creates a heatmap from ColdMap clusters - may not be necessary
#
#*******************************************************************************************

#******* makeColdMap ***********************************************************************
# Makes a "cold map" from a true/false occurrence matrix
# Assumes that the OM has been trimmed to an appropriate size prior to the function call
#
# Arguments
#  OM: occurrence matrix returned by ags.wrapper. 
#    Rows = gene names
#    Columns = pathway names
#  plotdata: list of info relevant to labeling and saving the plot
#   plotdir:  plot destination directory
#   plotbase:  base filename for the plot
#   plottitle:  title for all plots
#   pngres:  resolution for png files; defaults here to 600 dpi
#   pngsize:  c(width, length); defaults to c(15, 10)

makeColdMap = function (OM, plotdata = list(plotdir = NULL, plotbase = NULL, plottitle = NULL, pngres = NULL, pngsize= NULL)) { 
  
  # Set default print values
  if (is.null(plotdata$pngres)) {plotdata$pngres = 600}
  if (is.null(plotdata$pngsize)) {plotdata$pngsize = c(15,10)}
  
  # Convert OM from logical to numeric
  mynumOM = matrix(as.numeric(OM), nrow = nrow(OM), ncol = ncol(OM))
  rownames(mynumOM) = rownames(OM)
  colnames(mynumOM) = colnames(OM)
  
  message(sprintf("OM is a matrix with %i rows and %i Columns", nrow(OM), ncol(OM)))
  
  # Check assignments
  if(sum(colSums(OM) == colSums(mynumOM)) != ncol(OM)){
    message("Column data does not match in OM conversion from logical to numeric")
  }
  if(sum(rowSums(OM) == rowSums(mynumOM)) != nrow(OM)) {
    message("Row data does not match in OM conversion from logical to numeric")
  }
  
  n = length(unique(as.vector(mynumOM)))
  colorspec = paste("Blues",n,sep=":")
  
  message(sprintf("pngres = %i, width = %i, height = %i", plotdata$pngres, plotdata$pngsize[1], plotdata$pngsize[2]))
  message(sprintf("plotting to %s/%s", plotdata$plotdir, plotdata$plotbase))
  
  png(filename = sprintf('%s/%s_coldmap.png', plotdata$plotdir, plotdata$plotbase),
      width=plotdata$pngsize[1],height=plotdata$pngsize[2],units="in",res=plotdata$pngres)
  
#  myColdMap = makeHeatmap(ratiomat = mynumOM, attribs = NA, plottitle = plotdata$plottitle, clim.pct = 0.99, cexRow = 1,
  
  myColdMap = makeHeatmap(ratiomat = mynumOM, attribs = NA, plottitle = plotdata$plottitle, clim.pct = 1, cexRow = 1,
                          colorbrew = colorspec, labcoltype = "colnames", labRowType = 'genelabels')
  closeDevices()
  
  return(myColdMap)
}

#******* clusterGO *******************************************************************
# Converts a GO gene signature (GS) and occurrence matrix (OM) into a clustered GS/OM (CGS/COM) 
# Two criteria must be met for each level of clustering
#   1.  There must be a direct Child "is_a" or "part_of" Parent relationship between terms
#   2.  All Child genes must be present in the Parent term (this may be true by default, but check)
# The Cluster name will be identified by the "eldest" term, with the addition of 
#     (C-n) indicating that n offspring terms have been Clustered in that term
# An additional CGOxx term will be returned detailing the clustered GS/OM
#
# Arguments
#  GS - gene signature data.table returned by ags.wrapper(); must be GO annotation!
#  OM - occurrence matrix returned by ags.wrapper
#  ontology - limited to CC, MF, and BP. Other selections will fail.
# plotdata: list of info relevant to saving the clustering information
#   Default = NULL - don't write a file
#   plotdir:  plot destination directory; unspecified will write to getwd()
#   plotbase:  base filename for the plot
#  AS = NULL; optional - specify a version of the annotation set to use.
#     Note: annotation version number selection is not yet supported...
#
# Returns
# allCGS - clustered gene signature data.table of all parent:offspring relationships
# CGS - clustered gene signature data.table of parent terms
# COM - clustered occurrence matrix
# CAS - clustering association settings

clusterGO = function (GS, OM, ontology, clusterannot = c("parent", "eval", "enrichment"), plotdata = NULL, AS = NULL) {
  # Check the arguments
  exitErrors = FALSE
  if (!is.data.table(GS)) {
    message("GS argument should be data.table returned by ags.wrapper()")
    if (is.list(GS)) {message("List argument should be broken up to individual data.table")}
    exitErrors = TRUE
  }
  if(!is.matrix(OM)) {
    message("OM argument should be a matrix, as returned by ags.wrapper(). Note that this is an option, not default.")
    if (is.list(OM)) {message("List argument should be broken up to individual data.table")}
  }
  if(!(ontology %chin% c("CC", "BP", "MF"))) {
    message("ontology argument must be CC, BP, or MF")
    exitErrors = TRUE
  }
  if (exitErrors) {return("ERROR: check arguments")}
  
  # Load the associated GO pathway information
  # Note for now this is hard coded. Eventually, there will be a single data set with all information in it.
  GOdirectory = "/Users/lusardi/Documents/association_resources/"
  GOfile = "GO.2017.02.20_all.RData"
  load(paste(GOdirectory, GOfile, sep = '/'))
  go.Offspring = GOrelations[[paste("GO", ontology, sep = '')]]$Offspring
  go.Parents = GOrelations[[paste("GO", ontology, sep = '')]]$Parents
  go.Children = GOrelations[[paste("GO", ontology, sep = '')]]$Children
  
#  if (!("GO.db" %chin% installed.packages())) {
#    source("https://bioconductor.org/biocLite.R")
#    biocLite("GO.db", ask = FALSE)
#  }
#    library(GO.db)
  
  # Unpack the GS name (there may be a more elegant way to do this... I don't like hard coding!)
  set.Names = as.data.table(GS[,Set.Name])
  colnames(set.Names) = "Set.Name"
  set.Names[, "Ontol" := substr(Set.Name, 1,1)]
  set.Names[, "GOID" := paste(substr(Set.Name, 3,4), substr(Set.Name, 6,12), sep = ":")]
  set.Names[, "Description" := substring(Set.Name, 14)]
  
  # Merge set.Names with GS for easy reference
  expandedGS = merge(x = GS, y = set.Names, by = "Set.Name")
  
  # check for ontology match
  if (sum(set.Names$Ontol == substr(ontology,2,2)) != nrow(set.Names)) {
    message("ERROR:  clusterGO only handles one ontology at a time currently (MF, CC, or BP)")
    message(sprintf("GS contains %i different ontologies %s, %s specified in arguments",
                    length(unique(set.Names$Ontol)), paste(unique(set.Names$Ontol), collapse = ","),
                    ontology))
    return("ERROR: Ontology mismatch - no clustering performed")
  }
  
  # Get a list of parents
#  go.Parents = NULL
#  go.Offspring = NULL
#  go.Children = NULL
#  my.Parents = NULL
#  my.Offspring = NULL
#  my.Children = NULL
#  go.Parents = switch(ontology, MF = as.list(GOMFPARENTS), CC = as.list(GOCCPARENTS), BP = as.list(GOBPPARENTS))
#  go.Children = switch(ontology, MF = as.list(GOMFCHILDREN), CC = as.list(GOCCCHILDREN), BP = as.list(GOBPCHILDREN))
#  go.Offspring = switch(ontology, MF = as.list(GOMFOFFSPRING), CC = as.list(GOCCOFFSPRING), BP = as.list(GOBPOFFSPRING))
  my.Parents = go.Parents[set.Names$GOID]
  my.Offspring = go.Offspring[set.Names$GOID]
  my.Children = go.Children[set.Names$GOID]
  
  # Identify "senior" terms (no parent terms in our data set)
  # First, the terms that have no parents defined in the database (will be empty if there are none)
  Seniors = set.Names$GOID[!(set.Names$GOID %chin% names(my.Parents))]
  # Next, the terms that have no parents named "is_a" or "part_of" AND in set.Names$GOID
  nsig.Parents = lapply(my.Parents, function(x) sum((x %chin% set.Names$GOID) & (names(x) %chin% c("is_a", "part_of"))))
  Seniors = c(Seniors, names(nsig.Parents[nsig.Parents == 0]))
  
  # Identify terms that cluster under each senior term
  # Each list includes the senior member
  # Note:  "OFFSPRING" does not include relationship names; some are indirect, so recursive identification of direct child relations used
  clusterList = list()
  parentOffspring_dt = data.table(ParentGO = character(), OffspringGO = character(),
                                  minEflag = logical(), maxEnrichFlag = logical(),
                                  clustern = integer())
  for (mySenior in Seniors) {
    # Pathways with no offspring end up as "logical" list elements with NA value, causing an error
    if(!is.logical(my.Offspring[[mySenior]])) { 
      clusterList[[mySenior]] = c(mySenior, my.Offspring[[mySenior]][(my.Offspring[[mySenior]] %chin% set.Names$GOID)]) 
      
      minEVal = expandedGS[GOID %chin% clusterList[[mySenior]], E.Value] ==
                            min(expandedGS[GOID %chin% clusterList[[mySenior]], E.Value])
      maxEnrich = expandedGS[GOID %chin% clusterList[[mySenior]], FoldEnrich] ==
                             max(expandedGS[GOID %chin% clusterList[[mySenior]], FoldEnrich])
    } else { 
      clusterList[[mySenior]] = mySenior
      minEVal = TRUE
      maxEnrich = TRUE
    }
    # Create a parent/offspring summary table 
    parentOffspring_dt = rbind(parentOffspring_dt, cbind(ParentGO = rep(mySenior, length(clusterList[[mySenior]])),
                                                         OffspringGO = clusterList[[mySenior]],
                                                         minEflag = minEVal, maxEnrichFlag = maxEnrich,
                                                         clustern = rep(length(minEVal), length(minEVal))))
  }
   
  # Fill in the rest of the summary table
  parentOffspring_dt[, ':=' (Ontology = rep(ontology, nrow(parentOffspring_dt)),
                             Level = ifelse(ParentGO == OffspringGO, "Parent", "Offspring"),
                             OffspringDescription = sapply(OffspringGO, function(x) set.Names[GOID == x, Description]),
                             Set.Name = sapply(OffspringGO, function(x) set.Names[GOID == x, Set.Name]))]
  allCGS = setorder(merge(parentOffspring_dt, GS, by = "Set.Name"), ParentGO)
  
  # Confirm each pathway is included in at least one cluster
  allpways = unique(as.character(unlist(clusterList)))
  if(sum(allpways %chin% set.Names$GOID) != length(set.Names$GOID)) {
    message(sprintf("Pathway clustering is off... %i pathways considered, %i pathways clustered",
                    nrow(set.Names), length(allpways)))
    message("returning all data for evaluation")
    return(allCGS = allCGS)
  } 
  
  # Confirm that the Offspring genes are all contained within the Parent gene list
  for (mySenior in Seniors) {
    seniorGenes = OM[, set.Names[GOID == mySenior, Set.Name]]
    for (mypway in clusterList[[mySenior]]) {
      offspringGenes = myOM[, set.Names[GOID == mypway, Set.Name]]
      if(sum(offspringGenes & (seniorGenes == offspringGenes)) != sum(offspringGenes)) {
        message(sprintf("Genes in Pathway %s are not all contained in Pathway %s", mypway, mySenior))
      }
    }
  }
  
  # Create an expanded name set, and rename the CGS and COM columns according to clusterAnnot
  cluster.Names = set.Names[GOID %chin% Seniors,]
  
  # Create the clustered OM
  COM = OM[, allCGS[Level == "Parent", Set.Name]]
  return(list(allCGS = allCGS, COM = COM, CGS = allCGS[Level == "Parent",]))
}
#*******************************************************************************************

#******* coldClusters *******************************************************************
# coldClusters returns a table consisting of cluster sizes and scores based on dendrogram heights
# used to identify clusters of pathways regulated by the same genes for heatmapping
# Height = 0 will return pathways with identical gene signatures, at the max value, all pathways clustered
#
# Arguments
#   coldmap = structure returned by makeheatmap run on an OM 
#   OM = Occurrence Matrix used to generate the coldmap
#   minCllustSize = minimum number of pathways in a cluster
#   maxPctGenes = excludes clusters representing > % of genes
#   quantCut = sets the max distance quantile for clustering of clusters
#
# Returns
#   coldCluster_dt = data table with cluster scoring information
#
coldClusters = function(coldmap, OM, minClustSize = 1, maxPctGenes = 1, exploreClust = FALSE) {
  pathway_hcd = coldmap$Colv
  pathway_hc = as.hclust(pathway_hcd)
  heights = unique(pathway_hc$height)
  
#  library(stats)
  
  # Create optional plots exploring cluster heights/k-cuts
  if (exploreClust) {
    # Create a cluster plot directory if needed
    h.dir = paste(getwd(), "h.clusters", sep = "/")
    if (!dir.exists(h.dir)) {
      dir.create(h.dir)
    }
    
    # Easily distinguished color set
    colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
    
    index = 1
    for (myheight in heights) {
      # Assign membership by cluster height
      clustAssign = cutree(pathway_hc, h = myheight)
      ncolors = length(unique(clustAssign))
      labelColors = rep(colors, ncolors %/% length(colors) + 1)[1:ncolors]
      
      # From https://rpubs.com/gaston/dendrograms
      # Apply colors to each node in the dedndrogram
      colLab <- function(n) {
        if (is.leaf(n)) {
          a <- attributes(n)
          labCol <- labelColors[clustAssign[which(names(clustAssign) == a$label)]]
          attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
        }
        n
      }
      clusDendro = dendrapply(pathway_hcd, colLab)
      png_res = 450
      png(filename = sprintf('%s/%i.%s_%s.png', h.dir, index, "h.cluster", myheight), width=15,height=10,units="in",res=png_res)
      par(mar = c(30, 4, 4, 2))
      plot(clusDendro, cex = 0.3, main = sprintf("h = %1.3f", myheight), horiz = F)
      abline(a = myheight, b = 0, col = "red")
      
      dev.off()
      index = index + 1
    }
  }
  
  k.dir = paste(getwd(), "k.clusters", sep = "/")
  if (!dir.exists(k.dir)) {
    dir.create(k.dir)
  }
  index = 1
  for (mykcut in 1:length(pathway_hc$height)) {
    # Assign membership by cluster height
    clustAssign = cutree(pathway_hc, k = mykcut)
    ncolors = length(unique(clustAssign))
    labelColors = rep(colors, ncolors %/% length(colors) + 1)[1:ncolors]
    
    # From https://rpubs.com/gaston/dendrograms
    # Apply colors to each node in the dedndrogram
    colLab <- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clustAssign[which(names(clustAssign) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
      }
      n
    }
    clusDendro = dendrapply(pathway_hcd, colLab)
    png_res = 144
    png(filename = sprintf('%s/%i.%s_%s.png', k.dir, index, "k.cluster", mykcut), width=15,height=10,units="in",res=png_res)
    par(mar = c(30, 4, 4, 2))
    plot(clusDendro, cex = 0.3, main = sprintf("k = %i", mykcut))
    
    dev.off()
    index = index + 1
  }
  
  
  return(pathway_dend)
  
  
  # Create a table of clusters based on different height cuts
  cluster_dt = data.table(heightCut = numeric(0), cluster = integer(0), nPways = integer(0),
                          nGenes = integer(), ClusterScore = numeric(0), CS_pway = numeric(0),
                          Pathways = vector(), Genes = vector())
  
  for (myheight in heights) {
    clusters = cutree(pathway_dend, h = myheight)
    nclusters = unique(clusters)
    message(sprintf("myheight = %1.3f, nclusters = %i", myheight, length(nclusters)))
    for (mycluster in nclusters) {
      nPways = sum(clusters == mycluster)
      if (nPways >= minClustSize) {
        pways = pathway_dend$labels[cutree(pathway_dend, h = myheight) == mycluster]
        pway_str = paste(pways, sep = '', collapse = ';')
        if (is.na(match(pway_str, cluster_dt$Pathways))) {
          ngenes = 0
          geneIDs = logical(length(nrow(OM)))
          for (mypw in pways) {
            ngenes = ngenes + sum(OM[,mypw])
            geneIDs = geneIDs | (OM[, mypw] == 1)
          }
          geneList = rownames(OM)[geneIDs]
          gene_str = paste(geneList, sep = '', collapse = ';')
          ClusterScore = ngenes/nrow(OM)
          row2add = data.table(heightCut = myheight, cluster = mycluster, nPways = nPways,
                               nGenes = sum(geneIDs), ClusterScore = ClusterScore, CS_pway = ClusterScore/nPways,
                               Pathways = pway_str, Genes = gene_str)
          cluster_dt = rbindlist(list(cluster_dt, row2add))
        }
      }
    }
  }
  cluster_dt[, clustID := .I]
  
  # Create a levelplot of distance matrices to show cluster hierarchy
  clusterInclude_dt = cluster_dt[nGenes <= nrow(OM)*maxPctGenes,]
  rows2eval = nrow(clusterInclude_dt)
  clusterOM = matrix(data = logical(), nrow = nrow(OM), ncol = rows2eval)
  rownames(clusterOM) = rownames(OM)
  message("clusterInclude_dt has ", nrow(clusterInclude_dt), " rows\n")
  for (cluster in 1:rows2eval) {
    clustGenes = unlist(strsplit(clusterInclude_dt[cluster, Genes], split = ';'))
    clustLogic = rownames(clusterOM) %chin% clustGenes 
    clusterOM[, cluster] = clustLogic
  }
  OMdist = as.matrix(dist(t(clusterOM),diag=T))
  levelplot(OMdist[1:ncol(OMdist),ncol(OMdist):1])
  
  # Calculate clustering based on quantile cutoff
  
  
  return(list(cluster_dt = cluster_dt, clusterInclude_dt = clusterInclude_dt, clusterOM = t(clusterOM), OMdist = OMdist))
  
  return(cluster_dt)
  
  
}

#******* mergeOM ***********************************************************************
# Merges occurrence matrices (OM). Genes/pathways are included according to geneRules, minPways, minGenes
#
# Arguments
#  OM_ls - list of OMs. List elements must be named
#  mergeOnt - vector OMs to merge. Vector elements must be named, name must match the OM_ls element names
#           - element value must be the idtype from agswrapper - in Settings[[Ontology]]$idtype
#           - for example:  mergeOnt = c("kegg" = "Entrez.ID", "GOBP" = "Uniprot")
#  geneRules - Defaults to "intersection", includes only genes in all OMs specified in mergeOnt
#            - Also accepts "union", which includes all genes in any OM specified in mergeOnt
#            - Also accepts a vector of gene names to inlcude in the merged OM
#  background - this should be the same list used for pway analysis correlating assay targets, Uniprot, Entrez, etc
#  geneAnno - the annotation source that should be used for the mergedOM; must correspond to one of the columns in background
#
# Not implemented.... 
#  minPways - minimum number of pathways that a gene must be present in for inclusion in merged OM
#           - defaults to 1
#  minGenes - minimum number of genes that a pathway must have to include in the merged OM
#           - defaults to 1
#
# Returns
#  mergedOM - logical matrix with colnames = pathway names, rownames = genenames

mergeOM = function(OM_ls, mergeOnt, geneRules = "intersection", background, geneAnno, minPways = 1, minGenes = 1) {
  # Check OM and merge arguments
  if (sum(names(mergeOnt) %chin% names(OM_ls)) != length(mergeOnt)) {
    message("ERROR (mergeOM): elements of mergeOnt MUST match names of OM_ls")
    return("mergeOM ERROR")
  }
  
  # Confirm that the geneAnno specification is present in background
  if (!(geneAnno %chin% colnames(background))) {
    message("ERROR (mergeOM): geneAnno must match a column in background")
    return("mergeOM ERROR")
  }
  
  # Confirm that the background set includes all of the OM idtypes
  idtypes = unique(mergeOnt)
  idok = TRUE
  for (myid in idtypes) {
    if (!(myid %chin% colnames(background))) {
      message("ERROR (mergeOM): mergeOnt element names must match a column in background")
      idok = FALSE
    }
    if (!idok) {
      message("ERROR (mergeOM): Exiting due to missing idtypes")
      return("mergeOM ERROR")
    }
  }
  
  # Create a list of genes to include in the mergeOM
  incGenes = unique(background[, get(geneAnno)])
  # Create a list of uniformly annotated OMs
  OM2Merge = list()
  
  if (geneRules %chin% c("intersection", "union")) {
    for (myOnt in names(mergeOnt)) {
      myOM = OM_ls[[myOnt]]
      ontAnno = mergeOnt[myOnt]
      if (!(ontAnno %chin% colnames(background))) {
        message(sprintf("ERROR (mergeOM): idtype %s for %s not present in background", ontAnno, myanno))
        return("mergeOM ERROR")
      }
      
      if (ontAnno != geneAnno) {
        # make sure that there is a single geneAnno for each ontAnno
        # limit the background set to the two annotations in question so that duplicates for other annotations are excluded
        trimBkgnd = unique(background[, mget(c(ontAnno, geneAnno))])
        if (length(trimBkgnd[get(ontAnno) %in% rownames(myOM), get(ontAnno)]) > nrow(myOM)) {
          # deal with extra rows. Take duplicates in the geneAnno column and concatenate them with a "|"
          extras = trimBkgnd[, .N, by = get(ontAnno)][N > 1, get]
          geneAnno_term = unlist(lapply(extras, function(X) paste(trimBkgnd[get(ontAnno) == X, get(geneAnno)], collapse = '|')))
          trimBkgnd = trimBkgnd[!(get(ontAnno) %in% extras), ]
          clumped = data.table(extras, geneAnno_term)
          colnames(clumped) = c(ontAnno, geneAnno)
          trimBkgnd = rbind(trimBkgnd, clumped)
        }
        rownames(myOM) = trimBkgnd[get(ontAnno) %in% rownames(myOM), get(geneAnno)]
      }
      OM2Merge[[myOnt]] = myOM
      
      # Get a list of genes in myOnt in geneAnno
      incGenes = switch (geneRules,
                         intersection = intersect(incGenes, rownames(myOM)),
                         union = union(incGenes, rownames(myOM)) )
      message(sprintf("%s has %i rows, incGenes now %i long", myOnt, nrow(myOM), length(incGenes)))
    }
  } else if (is.vector(geneRules) & length(geneRules) > 1) { 
    # Make sure the included genes match the geneAnno and are present in background
    if (sum(geneRules %chin% background[, get(geneAnno)] == length(geneRules))) { 
      incGenes = geneRules    # Use the list of genes specified in geneRules for the OM
    } else {
      message(sprintf("ERROR (mergeOM): gene list specfied in geneRules does not match %s in background", geneAnno))
      return("mergeOM ERROR")
    }
    
    # Create a list of ontologies to merge with a uniform gene ID row names
    for (myOnt in names(mergeOnt)) {
      myOM = OM_ls[[myOnt]]
      ontAnno = mergeOnt[myOnt]
      if (!(ontAnno %chin% colnames(background))) {
        message(sprintf("ERROR (mergeOM): idtype %s for %s not present in background", ontAnno, myanno))
        return("mergeOM ERROR")
      }
      rownames(myOM) = background[get(ontAnno) %chin% rownames(myOM), get(geneAnno)]
      OM2Merge[[myOnt]] = myOM
    }
  } else {
    message("ERROR (mergeOM):  geneRules should be intersection, union, or a vector of genes")
    return("mergeOM ERROR")
  }
  
  # Create the merged OM
  mergedOM = matrix(nrow = length(incGenes), ncol = 0)
  rownames(mergedOM) = incGenes
  for (myOnt in names(mergeOnt)) {
    myOM = OM2Merge[[myOnt]]
    myOM = myOM[rownames(myOM) %chin% incGenes,]
    mergedOM = cbind(mergedOM, myOM[match(rownames(mergedOM), rownames(myOM)),])
  }
  
  # Limit the merged OM as specified in minPways, minGenes
  # Not implemented
  
  return(mergedOM)
}

#******* End mergeOM ***********************************************************************