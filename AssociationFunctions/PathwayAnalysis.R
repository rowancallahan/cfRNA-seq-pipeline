#####################################################################################################
# Set of functions for pathway analysis
#
# add_enrichment_ratio_column   --- Adds a column to the output of gs.wrapper containing enrichment ratio
# add_rank_ratio_column         --- Adds a column to the output of gs.wrapper containing rank ratio
# Remove_gene_lists_of_length_one_from_setlists  --- Function to scrub elements with one gene from setlist
# gs.wrapper.deprecated (TL)                    --- Add to me
# GMT_from_BioPax               --- Add to me
# GMT_from_GOowl                --- Add to me
# gs.pval                       --- Add to me
#
# Author: Julja Burchard, Mark Fisher
# Started - 2016
#####################################################################################################


add_enrichment_ratio_column <-
function(pathway_annotation_dt){
  pathway_annotation_dt[,enrichment_ratio := (Overlap/Sig.Size)/(Set.Size/Uni.Size)]
  return(pathway_annotation_dt)
}
add_rank_ratio_column <-
function(pathway_annotation_dt){
  tmp_dt = pathway_annotation_dt[order(Condx, Source, E.Value)]
  final_dt = NULL
  for (i in unique(tmp_dt[["Condx"]])){
    for(j in unique(tmp_dt[["Source"]])){
      sub_dt = tmp_dt[tmp_dt[["Source"]]==j & tmp_dt[["Condx"]]==i,]
      if (nrow(sub_dt)>0){
        rank_order = sort(tmp_dt[tmp_dt[["Source"]]==j & tmp_dt[["Condx"]]==i,][["E.Value"]])
        ranks = c(1:length(rank_order))
        contributing_dt = as.data.table(cbind(sub_dt,ranks))
        final_dt = rbind(final_dt, contributing_dt)
      }
    }
  }
  colnames(final_dt) = c(colnames(pathway_annotation_dt), "Rank")
  return (final_dt)
}
Remove_gene_lists_of_length_one_from_setlists <-
function(setlist){
  #If one of the signature sets from among the list that you're passing to gs.wrapper has length <= 1, gs.wrapper will error out. This function filters those sets out.
  #setlist is the list of gene signature lists that you pass to gw.wrapper
  final_list = list()
  for (i in 1:length(setlist)){
    if(length(setlist[[i]])>1){
      final_list[[length(final_list)+1]] = setlist[[i]]
      names(final_list)[length(final_list)] = names(setlist)[i]
    }
  }
  return(final_list)
}
gs.wrapper.deprecated <-
function(setlist, backgroundset, include.identifiers=FALSE, 
                      filename=NULL, resourceset=NULL, 
                      functiondir=NULL, resourcedir=NULL, return.OM=FALSE,
                      gs.range=c(25,2500), ecut=0.05, ocut=5, anno.uni=TRUE,  
                      annolist = c('kegg','msigdb5','GOBP','GOMF','GOCC',
                                   'pC','nci','BioCarta','Reactome'), 
                      idtypes = c(rep('Entrez.ID',5),rep('Anno.Symbol',4)),
                      use.gsr = c(F,T,T,T,T,T,F,F,F) ){
  # function to run basic pathway annotation on an input list of gene sets, using predefined annotation sets.  
  # gs.range limits the sizes of gene sets used for annotation
  # ecut sets the maximum E-value for annotations to return
  # ocut sets the minimum identifier number for annotations to return
  # return.OM is a flag indicating whether occurrence matrix is to be returned
  #   occurrence matrix == logical matrix, genes x sets
  # anno.uni is a flag: if TRUE, use the annotation source's universe size; 
  #   or if FALSE, use the overlap of backgroundset & anno universes
  # functiondir = optional directory with functions to be loaded (default wd)
  # resourcedir = optional directory with genesets to be loaded (default wd)
  # please adjust locations of functions and annotation files for your system
  # if resourceset is NULL (default) a default set of annotation sources
  #   is loaded, with all arguments ready as defaults
  # otherwise:
  # resourceset is a vector of filenames containing annotation
  #   sources as lists of named elements, each of which contains a vector
  #   of identifiers. Each file is an RData file. The filename should be
  #   the same as that of the list varname inside the file.
  # idtypes is a vector of the same length as resourceset that contains
  #   backgroundset colnames corresponding to each resourceset's id type
  # annolist is a vector of annotation source variable names as strings
  # use.gsr is a logical vector setting whether to use gs.range to restrict
  #   annotation source by gene set size
  # backgroundset is a data.table with fields ID, Entrez.ID and Anno.Symbol  
  # Entrez.ID must be HUMAN IDs to use default annotation sources
  # any other identifier types required by anno sources are in addnl columns
  # setlist is a list containing vectors of identifiers in backgroundset ID
  #  each element is named with the display name for the gene set 
  # Note that setlist MUST include a names() attribute.
  
  if(is.null(functiondir)) functiondir=getwd()
  if(is.null(resourcedir)) resourcedir=getwd()
  
  # imports
  source(file.path(functiondir,'gs.pval.R'))
  tmp=suppressPackageStartupMessages(require(data.table))
  
  # check arguments
  # setlist should be list of nonzero-length vectors
  if(!is.list(setlist) || sum(lengths(setlist)>0) != length(setlist) ){
    stop("setlist must be a list of nonzero-length vectors")
  }
  # get ready to save occurrence matrices if requested
  if( return.OM ){ OM_ls = NULL }
  
  # load annotations
  if( is.null(resourceset) ){
    message("Using default annotation sets")
    # annotations with EntrezGeneIDs
    load(file.path(resourcedir,"hsa.kegg.gs.RData")); kegg=hsa.kegg.gs
    kegg = lapply(hsa.kegg.gs,function(x) as.numeric(x))
    load(file.path(resourcedir,"msigdb.v5.0.entrez.gmt.RData"))
    msigdb5 = lapply(msigdb5,function(x) as.numeric(x))
    load(file.path(resourcedir,"GO_processed.gmt.RData"))
    # annotations with geneSymbols; not a robust identifier type
    load(file.path(resourcedir,"pathwayCommons.gmt.RData"));
    load(file.path(resourcedir,"NCI-Nature_Curated.gmt.RData"));nci=nci_gmt; rm(nci_gmt)
    load(file.path(resourcedir,"BioCarta.gmt.RData"));BioCarta=BioCarta_gmt;rm(BioCarta_gmt)
    load(file.path(resourcedir,"Reactome.gmt.RData"));Reactome=Reactome_gmt;rm(Reactome_gmt)
    # define universes for each gene set source
    uni = list(uni.kegg = unique(unlist(kegg,F,F)),
      uni.pC = unique(unlist(pC,F,F)),
      uni.msigdb5 = unique(unlist(msigdb5,F,F)),
      uni.GOBP = unique(unlist(GOBP,F,F)),
      uni.GOMF = unique(unlist(GOMF,F,F)),
      uni.GOCC = unique(unlist(GOCC,F,F)),
      uni.nci = toupper(unique(unlist(nci,F,F))),
      uni.BioCarta = toupper(unique(unlist(BioCarta,F,F))),
      uni.Reactome = toupper(unique(unlist(Reactome,F,F))) 
    )
  } else { #non-default load
    uni = NULL
    for( i in 1:length(resourceset) ){
      resourcename = sub('.RData','',resourceset[i])
      message("loading ",file.path(resourcedir,resourceset[i]))
      stopifnot( file.exists(file.path(resourcedir,resourceset[i])) )
      load(file.path(resourcedir,resourceset[i]))
      uni[[ paste("uni",resourcename,sep='.') ]] = 
             unique(unlist(get(resourcename)))
    }
    for( src in annolist ){
      message(src," universe size: ",length(uni[[paste('uni',src,sep='.')]]))
    }
  }  
  # annotate up- and down-regulated selected genes in files
  out.dt = data.table(rn="",E.Value="",P.Value="",Set.Size="",Overlap="",Sig.Size="",Uni.Size="",Source="",Condx="",FoldEnrich="")
  if(include.identifiers){
    out.dt = cbind(out.dt,ID.list="",alt.ID.list="")
  }
  out.dt = out.dt[0,]
  # loop over anno sources 
  for( src in annolist ){
    a = which(annolist==src)
    # set anno parameters
    # range of gene set sizes
    if( use.gsr[a] ){ min.size=gs.range[1]; max.size=gs.range[2]
    }else{ min.size=ocut; max.size=Inf }
    # ID type(s) to use
    tmp = backgroundset[,get(idtypes[a])]
    alt.idtype = "Anno.Symbol" # default == Symbol for human readability
    if( idtypes[a]=="Anno.Symbol" ){ alt.idtype = "Entrez.ID" }
    # universe size to use
    if( anno.uni ){
      uni.size=length(uni[[paste('uni',src,sep='.')]])
    } else { uni.size = NULL }
    # annotate each set
    for( i in 1:length(setlist) ){
      signame = names(setlist)[i]; sig = setlist[[i]]
      if(include.identifiers){
        anno = gs.pval(gs=get(src), rnames=tmp, 
            is.sig = backgroundset$ID %in% sig, uni.size=uni.size,
            min.size=min.size, max.size=max.size, include.identifiers=TRUE, 
            alt.rnames=backgroundset$Anno.Symbol,return.OM=return.OM)
      } else {
        anno = gs.pval(gs=get(src), rnames=tmp, 
            is.sig = backgroundset$ID %in% sig, uni.size=uni.size, 
            min.size=min.size, max.size=max.size, include.identifiers=FALSE,
            return.OM=return.OM)
      }
      GS = anno$GS
      if( return.OM ){ 
        OM_ls = c(OM_ls, anno['OM'])
        names(OM_ls)[length(OM_ls)] = signame 
      }
      sig.n = sum(GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut,na.rm=T)
      print(paste(signame, src,"uni.size",attr(GS,"Uni.Size"),"sig.n", sig.n))
      if(sig.n){
        # convert slice of GS to data.frame to preserve structure if n==1
        out.dt = rbindlist(list(out.dt,data.table(GS,keep.rownames=T)[GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut,]),fill=T)
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Sig.Size"] = attr(GS,"Sig.Size")
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Uni.Size"] = attr(GS,"Uni.Size")
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Source"] = src
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Condx"] = signame
        out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"FoldEnrich"] = 
      ( as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Overlap"]) /
        as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Sig.Size"]) ) /
      ( as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Set.Size"]) /
        as.numeric(out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"Uni.Size"]) ) 
        if(include.identifiers){
          out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"ID.list"] = attr(GS,"ID.list")[GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut]
          out.dt[(nrow(out.dt)-sig.n+1):(nrow(out.dt)),"alt.ID.list"] = attr(GS,"alt.ID.list")[GS[,"E.Value"]<ecut & GS[,"Overlap"]>=ocut]
        }
      }
    }
  }
  names(out.dt)[names(out.dt)=="rn"] = "Set.Name"
  if(!is.null(filename)){
    write.csv(out.dt,file=filename,row.names=F,quote=F)
  }
  return(list(out.dt=out.dt, OM_ls=OM_ls))
  
}
GMT_from_BioPax <-
function(biopaxstring){
  # extracts gene sets from OWL file as gene lists in GMT format
  # currently set up for BioPax2 only
  # returns list in GMT format
  require('rBiopaxParser')
  # process argument
  if( is.character(biopaxstring) ){ # assume file name
    biopax = readBiopax(biopaxstring)
  } else { biopax = biopaxstring }
  # check structure
  if( is.null(names(biopax)) | sum(names(biopax)=="dt")==0 | 
      !"biopax" %in% class(biopax) ){
    stop(paste(biopaxstring,"appears not to be a BioPax object"))
  }
  pw_component_list = listPathwayComponents(biopax, listInstances(biopax,class="pathway")$id, returnIDonly = T)
  if (length(pw_component_list) == 0) {
      stop("Pathways seem to have no pathway components")
  }
  if( sum(biopax$dt$class %chin% c("control", "catalysis", "modulation", "templatereactionregulation")) == 0 ){
      stop("Pathways have no regulatory components" )
  }
  # find pathways and pathway components
  p2pc = subset(biopax$dt,class=='pathway' & property=='PATHWAY-COMPONENTS',c('id','property_attr_value'))
  p2pc$property_attr_value = sub('^#','',p2pc$property_attr_value)
  # find next ID downstream of relevant pathway components
  pc2x = subset(biopax$dt,property %chin% toupper(c('controller','controlled','left','right','product')) & id %chin% p2pc$property_attr_value,c('class','id','property','property_attr_value'))
  setnames(pc2x,'property_attr_value','next_id')
  setnames(pc2x,'id','property_attr_value')
  pc2x$next_id = sub('^#','',pc2x$next_id)
  # track next IDs back to pathways 
  p2next = merge(p2pc,pc2x,by='property_attr_value',all.x=T,allow.cartesian=T)
  setnames(p2next,'property_attr_value','pathway_component')
  # at this point, some IDs are one step from a physical entity, others farther
  # retrieve and merge next layer..
  pc2pe = subset(biopax$dt, property=='PHYSICAL-ENTITY' & id %chin% pc2x$next_id, c('class','id','property','property_attr_value'))
  pc2n = subset(biopax$dt, property %chin% toupper(c('controller','controlled','left','right','product')) & id %chin% pc2x$next_id, c('class','id','property','property_attr_value'))
  setnames(pc2n,old=c('id','property_attr_value','class','property'),new=c('next_id','id3','next_class','next_property'))
  setnames(pc2pe,'id','next_id')
  # second round of mapping for next-layer IDs mapping to non-physical entities
  pc2n$id3 = sub('^#','',pc2n$id3)
  pc2pe2= subset(biopax$dt, property=='PHYSICAL-ENTITY' & id %chin% pc2n$id3, c('class','id','property','property_attr_value'))
  # bring in first-layer IDs for mapping back to pathways
  setnames(pc2pe2,'id','id3')
  pc2pe2n = merge(pc2n,subset(pc2pe2,select=c('id3','class','property','property_attr_value')),by='id3',all.y=TRUE)
  # check retrieval on remaining non-physical mappings
  pc2n2= subset(biopax$dt, property %chin% toupper(c('controller','controlled','left','right','product')) & id %chin% pc2n$property_attr_value, c('class','id','property','property_attr_value')) # expect 0 rows
  if( length(pc2n2$id)>0 ){
    warning(paste(length(unique(pc2n2$id)),"recursive mappings not pursued"))}
  # combine harvested mappings; adding mappings doesn't depend on new entities
  pc2pe.2x = rbindlist(list(pc2pe,pc2pe2n),use.names=T,fill=T)
  pc2pe.2x$property_attr_value = sub('^#','',pc2pe.2x$property_attr_value)
  pc2pe.2x = pc2pe.2x[!data.table:::duplicated.data.table(pc2pe.2x,by=NULL)]
  # pull mapped entities
  mymk.pe = biopax$dt$property_attr == "rdf:datatype" & biopax$dt$property == "NAME" & biopax$dt$class !='complex' & ! sub('^#','',biopax$dt$property_attr_value) %chin% pc2pe.2x$property_attr_value & biopax$dt$id %chin% pc2pe.2x$property_attr_value
  # pull mapped complexes
  pc2c = subset(biopax$dt,property_attr == "rdf:resource" & property %chin% c("COMPONENTS","PHYSICAL-ENTITY") & ! sub('^#','',property_attr_value) %chin% pc2pe.2x$property_attr_value & id %chin% pc2pe.2x$property_attr_value,c('id','property_attr_value'))
  pc2c = pc2c[!duplicated(pc2c)]
  setnames(pc2c,old=names(pc2c),new=c('property_attr_value','id4'))
  pc2c$id4 = sub('^#','',pc2c$id4)
  pc2c2pe = subset(biopax$dt, property=='PHYSICAL-ENTITY' & id %chin% pc2c$id4, c('class','id','property','property_attr_value'))
  setnames(pc2c2pe,old=c('id','property_attr_value'),new=c('id4','component_prop_attr_val'))
  pc2c2pe$component_prop_attr_val = sub('^#','',pc2c2pe$component_prop_attr_val)
  # bring in previous IDs for mapping back to pathways
  pc2pe2c = merge(pc2c,pc2c2pe,by='id4',all=TRUE)
  # pull entities mapped within complexes
  mymk.cpe= biopax$dt$property_attr == "rdf:datatype" & biopax$dt$property == "NAME" & biopax$dt$class !='complex' & ! sub('^#','',biopax$dt$property_attr_value) %chin% pc2pe2c$component_prop_attr_val & biopax$dt$id %chin% pc2pe2c$component_prop_attr_val
  # entity to name translation
  pe2name = subset(biopax$dt,mymk.pe|mymk.cpe,c('id','property_value'))
  setnames(pe2name,old=names(pe2name),new=c("property_attr_value","name"))
  # bring final mappings and entity names back to pathways
  pc2pe2name1 = merge(subset(pc2pe.2x,!pc2pe.2x$property_attr_value %chin% pc2pe2c$property_attr_value,c('property_attr_value','next_id','id3')),pe2name,all.x=T,by='property_attr_value')
  setnames(pe2name,'property_attr_value','component_prop_attr_val')
  pc2pe2name2 = merge(merge(subset(pc2pe.2x,pc2pe.2x$property_attr_value %chin% pc2pe2c$property_attr_value,c('property_attr_value','next_id','id3')),subset(pc2pe2c,select=c('property_attr_value','id4','component_prop_attr_val')),all.x=T,by='property_attr_value',allow.cartesian=T),pe2name,all.x=T,by='component_prop_attr_val')
  pc2pe2name = rbindlist(list(pc2pe2name1,subset(pc2pe2name2,select=names(pc2pe2name1))),use.names=T,fill=T)
  pc2pe2name = pc2pe2name[!data.table:::duplicated.data.table(pc2pe2name,by=NULL)]
  # bring mappings and names back to pathway IDs
  p2name = merge(p2next,pc2pe2name,all.x=T,by='next_id',allow.cartesian=T)
  p2name = p2name[!is.na(p2name$name)]
  # convert to GMT format
  pw_list = listInstances(biopax,class="pathway")
  mymk = pw_list$id %chin% unique(p2name$id)
  if( sum(mymk) < length(pw_list$id) ){
    warning(paste(sum(mymk),"of",length(pw_list$id),"pathways gave mapped named physical entities"))
    pw_list = pw_list[mymk,]
  }
  pw_gmt = vector(mode="list",length=length(pw_list$id))
  names(pw_gmt) = pw_list$name
  for(i in 1:length(pw_list$id) ){
    pw_gmt[[i]] = unique(p2name$name[p2name$id==pw_list$id[i]])
  }

  return(pw_gmt)
}
GMT_from_GOowl <-
function(GOowlstring,GOannostring,trim_regex=NULL){
  # TODO: add alternate reading, and argument protection
  # reads in GO ontology in GOowl format and as parsed in GO.db
  #  and reads in species-specific gene mappings from GOanno file
  # optionally trims gene identifiers using trim_regex
  # returns list of GMT-format lists @ for GO BP, MF and CC
  require("XML")
  require("data.table")
  # import processed GO tree
  require("GO.db")
  GOtags = list(P = 'BP', F = 'MF', C = 'CC')
  # read go OWL
  tmp = xmlRoot(xmlTreeParse(GOowlstring))
  tmp = tmp$children[xmlSApply(tmp,xmlName)=='Class']
  # fill GO ontology data table
  GOtab = data.table(rowcount=1:length(tmp),id="",label="",hasOBONamespace="", subClassOf="",key="rowcount")
  for(i in 1:length(tmp)){
    tmp2 = xmlApply(tmp[[i]],xmlValue)
    GOtab[i,'deprecated'] = length(xmlElementsByTagName(tmp[[i]],'deprecated'))>0
    # parse rest if not deprecated (data may not be available)
    if( !GOtab[i,get('deprecated')] ){
      GOtab[i,c("id","label","hasOBONamespace")] = tmp2[c("id","label","hasOBONamespace")]
      GOtab[i,'subClassOf'] = paste(sub('^.*\\/(GO_[0-9]+).*$','\\1',xmlElementsByTagName(tmp[[i]],'subClassOf')),sep='',collapse=';')
    }
  }
  GOtab = subset(GOtab,subset=deprecated==F,select=c("id","label","hasOBONamespace","subClassOf"))
  rm(tmp,tmp2)
  # read in tab-delimited GO-gene annotation in GAF format
  #  rows with NOT in column 4 refer to FALSE mappings to exclude
  tmp = fread(GOannostring,skip='\t')
  tmp = tmp[!grepl('NOT',tmp$V4)]; setkeyv(tmp,c('V5','V9','V2','V1'))
  tmp = subset(tmp,select=paste('V',c(1,2,3,5,9),sep=''))
  tmp = tmp[!data.table:::duplicated.data.table(tmp,by=NULL)]
  # join data tables to form GMT
  GOtype = unique(GOtab$hasOBONamespace)
  names(GOtype) = toupper(sub('^[^_]+_(.).*$','\\1',GOtype))#matches GOtab
  GO_gmt = vector(mode="list",length=length(GOtype))
  names(GO_gmt) = paste("GO",GOtags,sep='')
  for(j in 1:length(GOtype)){
    GOdown = as.list(get(paste("GO",GOtags[names(GOtype)[j]],"OFFSPRING",sep='')))
    tmp2 = GOtab[hasOBONamespace==GOtype[j]]
    setnames(tmp2,old=c('id','hasOBONamespace'),new=c('V5','V9'))
    tmp2$V9 = names(GOtype)[j]
    tmp2 = merge(tmp2,tmp,by=c('V5','V9'))
    # group ID columns
    tmp3 = as.data.table(cbind(paste(tmp2$V9,tmp2$V5,tmp2$label,sep=":"),paste(tmp2$V1,tmp2$V2,tmp2$V3,sep="|"),tmp2$V5))
    lfun = function(x){if(length(x)==0){I(list(''))}else{I(list(x))}}
    GOj = unique(tmp3$V1); GOjID = sub('^.:(GO:[0-9]+):.*$','\\1',GOj)
    GO_gmt[[j]] = vector(mode='list',length=length(GOj))
    # stuff list with genes assoc w each GO term and its offspring
    for(i in 1:length(GOj)){
      GOset = c(GOjID[i],GOdown[[ GOjID[i] ]])
      GO_gmt[[j]][i] = lfun( unique(tmp3$V2[tmp3$V3 %in% GOset]) )
    }
    names(GO_gmt[[j]]) = GOj
  }
  if( !is.null(trim_regex) ){
    for(j in 1:length(GOtype)){
      GO_gmt[[j]] = lapply(GO_gmt[[j]],function(x){sub(trim_regex,'\\1',x)})
    }
  }
  return(GO_gmt)
}
gs.pval <-
function(gs, rnames, is.sig, uni.size=NULL, min.size=0, max.size=Inf, include.identifiers=FALSE, alt.rnames=NULL, return.OM=FALSE){
  # Calculate hypergeometric p- and E-values of enrichment for gene sets
  # Code adapted from EnrichmentBrowser
  # Numeric gene identifiers are recommended for analysis speed
  # gs = list read from GMT file: elems contain sets, names are set names
  # rnames = unique identifiers of genomic features==rows of input data
  # is.sig = logical marking identifiers to retain
  # uni.size = size of potential gene universe
  #   if NULL, uni.size is counted directly in enrichment calculations
  #   if size of annotation source is supplied, it's used instead of counting
  #   for near-genome size sets, true anno source size may be best background
  #   if measurement platform is limited or biased, direct calc from rnames
  #    applied to gs may be best as that will give the intersection size
  # include.identifiers causes overlap gene identifiers to be returned
  # alt.rnames are identifers alternative to rnames, as vector in same order
  #   provided for returning as overlap gene alternate identifiers
  # return.OM is a flag indicating whether occurrence matrix is to be returned
  #   occurrence matrix == logical matrix, genes x sets
  
  # function for collapsing identifier lists
  semistring = function(myvec){paste(names(myvec)[as.logical(myvec)],sep='',collapse=';')}

  # check arguments
  if( !is.list(gs) | length(gs)==0 | length(names(gs)) != length(gs) ){
    stop('"gs" does not appear to be a gene set') }
  if( length(rnames) != length(is.sig) ){
    stop('feature identifier "rnames" and mask "is.sig" must have same length')}
  if( !is.null(alt.rnames) & length(rnames) != length(alt.rnames) ){
    stop('feature identifier "rnames" and alternative identifier "alt.rnames" must have same length')
  }
  if( length(unique(rnames)) != length(rnames) ){
    # convert mask to ids
    sig = rnames[is.sig]
    # shrink rnames and is.sig to size of unique elements
    mymk = !duplicated(rnames)
    rnames = rnames[mymk]; is.sig = rnames %in% sig
    if(!is.null(alt.rnames) ){ alt.rnames = alt.rnames[mymk] }
  }
  if( !is.logical(include.identifiers) ){
    include.identifiers = as.logical(include.identifiers)
    warn('flag include.identifiers converted to logical') }
  if( mode(gs[[1]]) != mode(rnames) ){
    # convert to character for safety
    rnames = as.character(rnames)
    gs = lapply(X=gs, FUN=function(x){as.character(x)})
  }
  # use faster %chin% function in package 'data.table' for character gene IDs
  if( is.character(rnames) ){
    tmp=suppressPackageStartupMessages(require(data.table)) }

  # trim & transform gene set to matrix, and align features with input data
  gs.sizes = lengths(gs,use.names=F)
  gs = gs[gs.sizes>=min.size & gs.sizes<=max.size]
  gs.sizes = gs.sizes[gs.sizes>=min.size & gs.sizes<=max.size]
  if( length(gs)==0 ){ stop("No gene set with valid size!") }
  if( is.character(rnames) ){
    cmat = as.matrix(as.data.frame(x=lapply(X=gs, FUN=function(x){ rnames %chin% x })))
  } else {
    cmat = as.matrix(as.data.frame(x=lapply(X=gs, FUN=function(x){ rnames %in% x })))
  }
  rownames(cmat) = rnames # rows correspond to input identifier universe
  has.set = which(rowSums(cmat) > 0) # some input IDs may not be in anno source
  cmat = cmat[has.set,]   # limit calculations to intersect of input & anno
  if(length(has.set)<length(rnames)){
    rnames = rnames[has.set]
    is.sig = is.sig[has.set]
    if( !is.null(alt.rnames) ){ alt.rnames = alt.rnames[has.set] }
  }
  if( is.null(uni.size) ){# universe size not given: calculate from intersect
    uni.size = nrow(cmat)
  }
  # calculate hypergeometric p-values
  nr.sigs = sum(is.sig)
  sig.cmat = cmat & is.sig
  ovlp.sizes = colSums(sig.cmat)
  gs.sizes = colSums(cmat)
  # R uses: successes in sample, successes in population, _non_successes in population, sample size
  gs.ps = phyper(ovlp.sizes,gs.sizes,uni.size-gs.sizes,nr.sigs,lower.tail=FALSE)
  gs.ps[ovlp.sizes==0] = 1 #phyper says finding 0 of 0-3 is p<0.05?!
  # generate return structure
  gs.idx = sort.list(gs.ps)
  # Holm controls family-wise error rate, is valid under arbitrary assumptions
  #  and has more power than simple Bonferroni correction
  gs.es = p.adjust(gs.ps,method="holm")
  res.tbl = cbind(gs.es,gs.ps,gs.sizes,ovlp.sizes)
  colnames(res.tbl) = c("E.Value","P.Value","Set.Size","Overlap")
  rownames(res.tbl) = colnames(cmat)
  res.tbl = res.tbl[gs.idx,]
  attr(res.tbl,"Sig.Size") = nr.sigs
  attr(res.tbl,"Uni.Size") = uni.size
  if(include.identifiers){
    attr(res.tbl,"ID.list") = apply(X=sig.cmat,MARGIN=2,FUN=semistring)[gs.idx]
    if( !is.null(alt.rnames) ){
      tmp = sig.cmat; rownames(tmp) = alt.rnames
      attr(res.tbl,"alt.ID.list") = apply(X=tmp,MARGIN=2,FUN=semistring)[gs.idx]
    }
  }
  return( list(GS=res.tbl, OM=sig.cmat) )
}
