#######################################################################################################################
# Gene Set Parsing Files
# These are version2 parsers, adding to the versions in PathwayAnalysis.R
#
# GMT_from_GOowl2 = enhanced parser. For each ontology, add the following tables:
#         _hierarchy = kludge to save a contemporary list of hierarchical relationships between GO terms
#         _

############ GMT_from_GO_owl ########################################################################
GMT_from_GOowl <-
function(GOdirectory, GOowlstring, GOannostring, GOsavebase, trim_regex=NULL, verbose = F){
# TODO: add alternate reading, and argument protection
# reads in GO ontology in GOowl format and as parsed in GO.db
#  and reads in species-specific gene mappings from GOanno file
# TODO: add Evidence Code parsing
# TODO: Consider adding a switch to allow direct, indirect, or all hierarchies
#       Currrently - only allows direct relataionships
#
# Arguments
#  GOdirectory - annotation source directory
#  GOowlstring - OWL format xml file with ontology relationships.
#    -- need to decide whether this is duplicating effort from GO.db
#  GOannostring - GAF (gene annotation file), version 2.0 assumed
#    -- Note: any line beginning with ! should be ignored
#    -- Documented in http://geneontology.org/page/go-annotation-file-format-20
#  GOsavebase - base name of the file to save annotation data
#  trim_regex - optionally trims gene identifiers
#  verbose - if TRUE, writes info messages to terminal
  
# Saves and Returns
#  RData file with:
#    GO_gmt - list of GMT-format lists for GO, BP, MF, and CC
#    GOrelations - list by Aspect of Offspring, Parents, and Children terms
#    GO_ontol_dt - data table of ontology details by GOID
#    GOsettings - List of details used to create the above data
#      GOdir - source data directory)
#      GOowl - Name of the file for GO term details
#      GOgoa - Name of the file for GO gene/term associations
#      GOdb - information returned from GO.db

  require("XML")
  require("data.table")
  # import processed GO tree
  require("GO.db")
#  GOtags = list(P = 'BP', F = 'MF', C = 'CC')
  GOtags_dt = data.table(Aspect = c('P', 'F', 'C'),
                         Tag = c('BP', 'MF', 'CC'),
                         OBO = c('biological_process', 'molecular_function', 'cellular_component'),
                         Ont = c('GOBP', 'GOMF', 'GOCC'))
  
  # OWL file has the relationships among GO terms - the ontology. May duplicate GO.db
  # read go OWL
  tmp = xmlRoot(xmlTreeParse(GOowlstring))
  tmp = tmp$children[xmlSApply(tmp,xmlName)=='Class']
    
  # fill GO ontology data table
  GOontol_dt = data.table(rowcount=1:length(tmp),id="",label="",hasOBONamespace="", subClassOf="",key="rowcount")
  for(i in 1:length(tmp)){
    tmp2 = xmlApply(tmp[[i]],xmlValue)
    GOontol_dt[i,'deprecated'] = length(xmlElementsByTagName(tmp[[i]],'deprecated'))>0
    # parse rest if not deprecated (data may not be available)
    if( !GOontol_dt[i,get('deprecated')] ){
      GOontol_dt[i,c("id","label","hasOBONamespace")] = tmp2[c("id","label","hasOBONamespace")]
      GOontol_dt[i,'subClassOf'] = paste(sub('^.*\\/(GO_[0-9]+).*$','\\1',xmlElementsByTagName(tmp[[i]],'subClassOf'))    ,sep='',collapse=';')
    }
  }
  GOontol_dt = subset(GOontol_dt,subset=deprecated==F,select=c("id","label","hasOBONamespace","subClassOf"))
  # Adjust id column name
  setnames(GOontol_dt,old='id',new='GOID')
  if (verbose) {
    message(sprintf("Summary of ontology %s by column", GOowlstring))
    cols = colnames(GOontol_dt)
    for (i in cols) {
      cat("\n", i, "\n")
      print(GOontol_dt[, .N, by = get(i)][order(-N)])
    }
  }
  rm(tmp,tmp2)
  
  # Annotation File - associates Genes with an Ontology Term
  # read in tab-delimited GO-gene annotation in GAF format
  #  rows with NOT in column 4 refer to FALSE mappings to exclude
  # gaf-version: 2.0 column names
  gaf.2.0_colnames = c(v1 = "DB", v2 = "DB_ID", v3 = "DB_Symbol", v4 = "Qualifier", v5 = "GOID",
                       v6 = "DB_Refs", v7 = "EvidenceCode", v8 = "WithFrom", v9 = "Aspect", v10 = "DB_Name",
                       v11 = "DB_Synonym", v12 = "DB_ObjType", v13 = "Taxon", v14 = "Date", v15 = "AssignedBy",
                       v16 = "AnnotExt", v17 = "GeneProductID")
  
  annotations = fread(GOannostring,skip='\t')
  names(annotations) = gaf.2.0_colnames
  
  # Exclude annotations with qualifier "NOT"
  annotations = annotations[!grepl('NOT',annotations$Qualifier)]
  
  if (verbose) {
    message(sprintf("Summary of %s by column", GOannostring))
    for (i in gaf.2.0_colnames) {
      cat("\n", i, "\n")
      print(annotations[, .N, by = get(i)][order(-N)])
    }
  }
  setkeyv(annotations,c('GOID','Aspect','DB_ID','DB'))
  cols2keep = gaf.2.0_colnames[c(1,2,3,5,7,9)]
  annotations = subset(annotations,select=cols2keep)
  # Ask JB - wh duplicated.data.table rather than unique?
  # Seem to get identical results either way.
  annotations = annotations[!data.table:::duplicated.data.table(annotations,by=NULL)]
  
  # join data tables to form GMT
  # Separate by ontology
  GOontols = unique(GOontol_dt$hasOBONamespace)
  GOontol_dt[, Aspect := unlist(lapply(GOontol_dt$hasOBONamespace, function(X) GOtags_dt[OBO %chin% X, Aspect]))]
  GO_gmt = vector(mode="list",length=length(GOontols))
  # Gather obsolete terms for verification...
  obsoleteTerms = vector(mode="list",length=length(GOontols))
  names(GO_gmt) = GOtags_dt[OBO %chin% GOontols, Ont]
  names(obsoleteTerms) = GOtags_dt[OBO %chin% GOontols, Ont]
  
  # For debugging - just to GOCC - it's smallest
  for(myontol in names(GO_gmt)) {
    message("Calculating ", myontol," summary terms")
    message("This will take a while...")
    # Get a list of all offspring from GO.db
    GOoffspring = as.list(get(paste(myontol, "OFFSPRING", sep = '')))
    GOparents = as.list(get(paste(myontol, "PARENTS", sep = '')))
    GOchildren = as.list(get(paste(myontol, "CHILDREN", sep = '')))
    
    # Get the relevant ontology terms and name columns consistently with the annotations
    ontol = GOontol_dt[Aspect == GOtags_dt[Ont == myontol, Aspect]]
    mergedOA = merge(ontol,annotations,by=c('GOID','Aspect'))
    
    # Add combined terms for backwards compatability
    mergedOA[, DB3 := paste(DB, DB_ID, DB_Symbol, sep = '|')]
    mergedOA[, GO3 := paste(Aspect, GOID, label, sep = ':')]
    
    lfun = function(x){if(length(x)==0){I(list(''))}else{I(list(x))}}
    
    # Use ontol, not mergedOA, which only has terms with a gene directly assigned to it.
    # -- If a GO term is a summary of multiple child terms, it wont have any "solo" genes
    uniqueGOID = unique(ontol$GOID)
    
    # stuff list with genes assoc w each GO term and its direct offspring 
    directs = c("is_a", "part_of")
    
    for (myid in uniqueGOID) {
      GOoff = GOoffspring[[myid]]
      # Each term includes itself...
      GOset = myid
      # Create a vector of excluded offspring terms
      GOexclude = NULL
      
      # If a term is obsolete, GOoffspring will be NULL; report and move on; hopefully as we resolve our versions, this will not happen.
      if (is.null(GOoff)) {
        message(sprintf("Term %s, %s is obsolete or has been replaced and excluded from the dataset", myid, ontol[GOID == myid, label])) 
        obsoleteTerms[[myontol]] = unique(c(obsoleteTerms[[myontol]], myid))
        next   # Skip to next myid
      }
      
      # If a term has no offspring, it's the lowest level of the ontology, GOoffspring = NA, & there are no terms to merge.
      if (length(GOoff) > 0 & sum(is.na(GOoff)) > 0) {
        GO_gmt[[myontol]][myid] = lfun(unique(mergedOA[GOID %in% GOset, DB3]))
        next   # Skip to next myid
      }
        
      # Confirm that all GOoffspring have an 'is_a' or 'part_of' relationship all the way to myid
      left2test = copy(GOoff)
      while (length(left2test) > 0) {
        for (i in 1:length(left2test)) {
          myterm = left2test[i]
          myparents = GOparents[[myterm]][names(GOparents[[myterm]]) %chin% directs]
          if (sum(myparents %chin% GOset) > 0) {  # direct children get in. Take their kids too.
            mychildren = GOchildren[[myterm]][names(GOchildren[[myterm]]) %chin% directs]
            GOset = unique(c(GOset, myterm, mychildren))
            left2test = left2test[!(left2test %chin% c(myterm, mychildren))]
            break
          } else if (sum(myparents %chin% left2test) == 0) { # if no direct parents in remaining offspring - exclude
            GOexclude = c(GOexclude, myterm)
            left2test = left2test[!(left2test %chin% myterm)]
            break
          } else if (sum((myparents %chin% GOexclude) == (myparents %chin% left2test)) == length(myparents)) { # if all parents already excluded, then exclude
            GOexclude = c(GOexclude, myterm)
            left2test = left2test[!(left2test %chin% myterm)]
            break
          } 
        }  # End of for loop through left2test
      }  # end of while loop testing GOoffspring for inclusion 
        
      # confirm all terms still accounted for... 
      if (sum((GOoff %chin% GOset) | (GOoff %chin% GOexclude)) != length(GOoff)) {
        message(sprintf("BAD NEWS!!! You lost GOoff terms - fix your code; myid %s, myterm %s", myid, myterm))
        return()
      }
        
      if (sum(myid %chin% obsoleteTerms[[myontol]]) == 0) {
        GO_gmt[[myontol]][myid] = lfun(unique(mergedOA[GOID %in% GOset, DB3]))
      } 
    }  # End of for myid loop through uniqueGOID
  }
  
  if( !is.null(trim_regex) ){
    for(j in 1:length(GOtype)){
      GO_gmt[[j]] = lapply(GO_gmt[[j]],function(x){sub(trim_regex,'\\1',x)})
    }
  }
  
  # Create a complete GOrelations file
  GOrelations = list()
  for(myontol in c("GOBP", "GOMF", "GOCC")) {
    GOrelations[[myontol]][['Parents']] = as.list(get(paste(myontol, "PARENTS", sep = '')))
    GOrelations[[myontol]][['Offspring']] = as.list(get(paste(myontol, "OFFSPRING", sep = '')))
    GOrelations[[myontol]][['Children']] = as.list(get(paste(myontol, "CHILDREN", sep = '')))
  }
  
  GOsettings = list(GOdirectory = GOdirectory, GOowl = GOowlstring, GOgoa = GOannostring, GOdb = GO_dbInfo())
  GOReturn = list(GO_gmt = GO_gmt, GOrelations = GOrelations, GOontol_dt = GOontol_dt, GOsettings = GOsettings)
  alldatafile = paste(GOsavebase, "all.RData", sep = "_")
  save(GO_gmt, GOrelations, GOontol_dt, GOsettings, file = paste(GOdirectory, alldatafile, sep = '/'))
  # Keep a version with just the gmt data for backwards compatibility
  gmtdatafile = paste(GOsavebase, "gmt.RData", sep = "_")
  save(GO_gmt, file = paste(GOdirectory, gmtdatafile, sep = '/'))
  # Save versions by Association
  GOCC = GO_gmt$GOCC
    gmtdatafile = paste(GOsavebase, "GOCC", "gmt.RData", sep = "_")
    save(GOCC, file = paste(GOdirectory, gmtdatafile, sep = '/'))
  GOMF = GO_gmt$GOMF
    gmtdatafile = paste(GOsavebase, "GOMF", "gmt.RData", sep = "_")
    save(GOMF, file = paste(GOdirectory, gmtdatafile, sep = '/'))
  GOBP = GO_gmt$GOBP
    gmtdatafile = paste(GOsavebase, "GOBP", "gmt.RData", sep = "_")
    save(GOBP, file = paste(GOdirectory, gmtdatafile, sep = '/'))
  
  return(GOReturn)
}

############ GMT_from_GO_owl ########################################################################

