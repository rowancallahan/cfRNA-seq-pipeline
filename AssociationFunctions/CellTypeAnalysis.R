##################################################################################################################################
# Set of functions useful for cell specificity analysis using the Fantom5 phase1 database.
#
# getTopNCellTypes                   --- Generates list of data.tables with top N cell types
# getSortedDTwithSplitOfFirstCol     --- Generate and return a data table of the fantom5 data SORTED chromosome number and THEN by Start position
# ConvertFFlsOflsOfDtToDt            --- Converts list of list of data tables to one data.table of summarized output
# FFTopNCellTypesWrapper             --- Wraps all of the other functions in this script besides fetchTssRangesWithXFoldChangeInFFOntologyID
# fetchTssRangesWithXFoldChangeInFFOntologyID    --- Generates a list of matrices of entries with log2-fold ratios greater than x (or abs(log2(ratio)>x)), one element for each query/fantom5 ontology ID. Wraps generateLg2RatioToAveAllMat.
# generateLg2RatioToAveAllMat        --- Generate and return a matrix containing log2-scale ratios to average all for the numerical portion of the fantom5 data
#
# Author: Mark Fisher
# Started - October, 2016
#####################################################################################################

getTopNCellTypes <-
function(ff_dt=NULL,chromoNum, range_v, N, colStrWithChromLocations_str="00Annotation", ffDtPath=NULL, ffOntIdInstead=F, optRank_Mat=NULL, optURLdecode=F, ffSorted_dt=NULL) {
  #' Locates top N tpm counts from fantom5 database for each TSS within or closest to query range.
  #' @description Takes as input a genome position (chr,Start,Stop), looks for entries in ff_dt with an intersect with the query, and returns the top N tpm columns FOR EACH hit as a separate data.table. In other words, the output is a LIST of data.tables. If no matching entries are found, the function will search for the nearest entry and return the top N tpm columns of that. If there are two nearest neighbors, both will be returned.
  #' 
  #' @param ff_dt a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user.
  #' @param chromoNum a numerical or string of the chromosome (nb NOT "chr1", but "1") of the query position.
  #' @param range_v range_v a numeric() vector containing the Start and Stop bp position of the query (e.g., range_v = c(100,123))
  #' @param N a numeric() argument that retreives the top N largest tpm-containing columns in a given row of ff_dt
  #' @param colStrWithChromLocations_str a string containing the column name of the annotation column of ff_dt (i.e., the column that contains TSS positional information). In ff_dt, note that the position information needs to be parsed out, which this function('s dependency) does automatically. Default is "00Annotation".
  #' @param ffDtPath is a string containing the path to ff_dt if 
  #' @param ffOntIdInstead a boolean designating whether the users wants to substitute the full URL names of the cell types/samples with just the fantom5 ontology IDs instead.
  #' @param optRank_Mat an optional matrix containing the ranked values of all of the fantom5 tpm entries. Most useful and time-saving if lots of calls to the sort() function expected (i.e., if lots of hits are expected or if this function is being called in a for loop as with FFTopNCellTypesWrapper).
  #' @param optURLdecode a boolean designating whether the user wants to apply a URL-decoder function on the colnames of topCells_dt (i.e., on the column names of the sorted/ranked cell types/samples)
  #' @param ffSorted_dt an optional data.table containing ff_dt sorted by chromosome and then by start position. Most useful and time-saving if lots of calls to the sort() function expected (i.e., if lots of hits are expected or if this function is being called in a for loop as with FFTopNCellTypesWrapper).
  #' @return A list of data.tables with each element of the list named with the intersection of query and ff_dt as well as original ff_dt position information. Each data.table in the list contains the cell type/sample info. in rank order from 1:N in decreasing order, as well as the tpm count corresponding to it.
  #' @examples
  #' require(data.table)
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' names(ffDummy_dt) = c("Annotation", "tpmOne","tpmTwo.CNhs10619.10014-101C5", "tpmThree")
  #' cNum = 10
  #' rng_vec = c(1,1000)
  #' topN = 2
  #' ChromLocCol_str = "Annotation"
  #' test1 = getTopNCellTypes(ff_dt=ffDummy_dt,chromoNum=cNum, range_v=rng_vec, N=topN, colStrWithChromLocations_str=ChromLocCol_str, ffDtPath=NULL,ffOntIdInstead=F)
  #' test1.2 = getTopNCellTypes(ff_dt=ffDummy_dt,chromoNum=cNum, range_v=rng_vec, N=topN, colStrWithChromLocations_str=ChromLocCol_str, ffDtPath=NULL,ffOntIdInstead=T)
  #' test1.3 = getTopNCellTypes(ff_dt=ffDummy_dt,chromoNum=cNum, range_v=rng_vec, N=topN, colStrWithChromLocations_str=ChromLocCol_str, ffDtPath=NULL,ffOntIdInstead=T, optRankMatPath="/home/exacloud/lustre1/CompBio/genomic_resources/fantom5/human/ProcessedAndUsefulFiles/AllFFEntriesRanked.txt")
  #' ##Test with chromoNum as "X"##
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chrX:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' cNum = "X"
  #' test2 = getTopNCellTypes(ff_dt=ffDummy_dt,chromoNum=cNum, range_v=rng_vec, N=topN, colStrWithChromLocations_str=ChromLocCol_str, ffDtPath=NULL)
  #' ##Test with no ff_dt provided## #Attn
  #' test3 = getTopNCellTypes(ff_dt=NULL,chromoNum=cNum, range_v=rng_vec, N=topN, colStrWithChromLocations_str="00Annotation", ffDtPath="grep -v '^#' /Users/fishema/Git_repos/AbundanceFunctions/IgnoreThese/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt", optURLdecode=T)
  #' ff_dt = fread("grep -v '^#' /Users/fishema/Git_repos/AbundanceFunctions/IgnoreThese/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt")
  #' test4 = getTopNCellTypes(ff_dt=ff_dt,chromoNum=cNum, range_v=rng_vec, N=topN, colStrWithChromLocations_str="00Annotation", ffDtPath=NULL, optURLdecode=T)
  #' #Attn test on the battery of calls at the end of this script
  #' @export
  
  require(data.table)
  #browser()
  
  #If no ff_dt is provided, use the path provided instead
  if(is.null(ff_dt)){
    ff_dt = fread(ffDtPath)
  }
  
  #Check that N isn't greater than the number of cell types/samples available
  if(N > length(grep('tpm', names(ff_dt)))){
    return(NULL)
  }else{
    #Sort and parse the identifier column of ff_dt to enable looking for nearby TSS neighbors
    if(is.null(ffSorted_dt)){
      ffSorted_dt = getSortedDTwithSplitOfFirstCol(ff_dt=copy(ff_dt), nameOfCol1=colStrWithChromLocations_str)
    }
    
    #Check that range_v is valid
    if (length(range_v)!=2){
      stop("Too few or two many items in range_v. Should be a vector of length 2") #Attn
    }else{
      #Ceorce the range to go in correct order
      range_v = range(range_v)
      
      #Find which rows in the chromosome column of ffSorted_dt match our query chromosome
      chromHits = which(ffSorted_dt[,chr]==as.character(chromoNum))
      if(length(chromHits)<1){
        print(paste0("Invalid or non-matching chromosome number: ",chromoNum))
        return(NULL)
      }else{
        startRowChrom = chromHits[1]
        endRowChrom = chromHits[length(chromHits)]
        
        #Find the first row where Start of our query is >= Start in dt
        directHit = which(ffSorted_dt[startRowChrom:endRowChrom,Stop]>=range_v[1] & range_v[2]>= ffSorted_dt[startRowChrom:endRowChrom,Start])
        
        #Sort each row (or optionally retrieve ranked entries passed from wrapper function), and pull out the top N hits into a data.table called topCells_dt. Add topCells_dt to a list called matches_ls.
        #I think this has to be a for loop, because I have to sort each row independently, and there's no way to retain colnames under those circumstances
        if(length(directHit)>0){
          matches_ls = list()
          for (dh in 1:length(directHit)){
            
            #If rank_mat has been passed, use it instead #Attn needs to be adjusted for tie break
            if(!is.null(optRank_Mat)){
              rankColNums = as.numeric(optRank_Mat[startRowChrom:endRowChrom,,drop=F][directHit[dh],grep('tpm', colnames(optRank_Mat)),drop=F])<=N
              totalSuccessfullyRankedCols = sum(rankColNums)
              #If there was a tie, you may be a few short of a full N's worth of entries. If this happens, populate with the next N-totalSuccessfullyRankedCols hits in order of appearance #Feedback not sure whether this was the best thing to do, but I like it better than populating a Top20 matrix with more than 20 entries.
              if(totalSuccessfullyRankedCols<N){
                RankNums = as.numeric(optRank_Mat[startRowChrom:endRowChrom,,drop=F][directHit[dh],grep('tpm', colnames(optRank_Mat)),drop=F])
                SupplementaryRankColNumsToUse = which(RankNums==min(RankNums[!rankColNums]))[1:(N-totalSuccessfullyRankedCols)]
                #countAsTrue_v = c(which(rankColNums==TRUE),SupplementaryRankColNumsToUse)
                rankColNums[SupplementaryRankColNumsToUse] = TRUE
              }
              rankOrder = as.numeric(optRank_Mat[startRowChrom:endRowChrom,,drop=F][directHit[dh],grep('tpm', colnames(optRank_Mat)),drop=F][,rankColNums])
              #replace any greater than N with 0 and then with max +1
              rankOrder[rankOrder > N] = 0
              rankOrder[rankOrder ==0] = (max(rankOrder)+1):(max(rankOrder)+1+sum(rankOrder ==0))
              #there could still be ties in the <20 range, and these should be dealt with as well.
              rankOrder = rank(rankOrder,ties.method = "random")
              tmp = ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],grep('tpm', names(ffSorted_dt)),with=F]
              topCells_dt = tmp[,rankColNums,with=F]
              NamesOrdered = names(topCells_dt)[rankOrder]
              setcolorder(topCells_dt,NamesOrdered)
            }else{
              topCells_dt = sort(ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],grep('tpm', names(ffSorted_dt)),with=F], decreasing = TRUE)[,1:N,with=F]
            }
            
            #If ontology ID desired instead, return just those values instead
            if(ffOntIdInstead ==T){
              OntIDs = gsub('^.*CNhs.*\\.(.*$)','\\1',names(topCells_dt))
              names(topCells_dt) = OntIDs
            }
            
            #If URL decode desired, call the URL decode function
            if(optURLdecode==T){
              names(topCells_dt) = sapply(names(topCells_dt),URLdecode) #unlist(lapply(names(topCells_dt),URLdecode)) #Attn URLdecode hasn't been written yet
            }
            
            matches_ls[[length(matches_ls)+1]] = topCells_dt
            
            #Name the match in matches_ls with latest start, earlier stop (i.e., the intersect range) and the original TSS identifer position in ff_dt.
            LatestStart= max(range_v[1],ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],Start])
            EarliestStop = min(range_v[2],ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],Stop])
            OriginalFFdtIdentifier= ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],][[colStrWithChromLocations_str]]
            names(matches_ls)[length(matches_ls)] = paste(chromoNum,LatestStart,EarliestStop,"matched_region_OriginalColName",OriginalFFdtIdentifier,sep="_")
          }
          return(matches_ls)
        }
        if(length(directHit)==0){
          ##Look for whether the query start is beyond the very last Stop in ffSorted_dt for this chromosome
          if (range_v[1] > ffSorted_dt[startRowChrom:endRowChrom,][nrow(ffSorted_dt[startRowChrom:endRowChrom,]),Stop]){
            nearestUpstreamNeighborIdx = nrow(ffSorted_dt[startRowChrom:endRowChrom,])
            nearestDownstreamNeighborIdx = nearestUpstreamNeighborIdx
          }else{
            ##Look for first occassion where Start in the ffSort_dt is greater than Stop in the query
            nearestDownstreamNeighborIdx = which(ffSorted_dt[startRowChrom:endRowChrom,Start]>range_v[2])[1]
            if (nearestDownstreamNeighborIdx>1){ #if it's not the first index (which is the start of the chromsome or even of ffSorted_dt, let's look one back from there)
              nearestUpstreamNeighborIdx = nearestDownstreamNeighborIdx-1
            }else{
              nearestUpstreamNeighborIdx=nearestDownstreamNeighborIdx #otherwise, we'll keep it as is
            }
          }
          #is the Stop of the nearestUpstreamNeighbor closer to start of the query, or is the Start of nearestDownstreamNeighbor closer to the stop of the query?
          QstopDiffDownStart = ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Start] - range_v[2]
          QstartDiffUpStop= range_v[1] - ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Stop]
          if (QstopDiffDownStart < QstartDiffUpStop){ #start of nearest downstream neighbor is closer to the stop of the query
            newPositions = range(ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Start], ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Stop])
            return(getTopNCellTypes(ff_dt=ff_dt,chromoNum=chromoNum, range_v=newPositions, N=N, colStrWithChromLocations_str=colStrWithChromLocations_str,ffDtPath,ffOntIdInstead,optRank_Mat,optURLdecode,ffSorted_dt))
          }
          if(QstartDiffUpStop < QstopDiffDownStart){ #stop of the nearest upstream is closer to the start of the query
            newPositions = range(ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Start], ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Stop])
            return(getTopNCellTypes(ff_dt=ff_dt,chromoNum=chromoNum, range_v=newPositions, N=N, colStrWithChromLocations_str=colStrWithChromLocations_str,ffDtPath,ffOntIdInstead,optRank_Mat,optURLdecode,ffSorted_dt))
          }
          if(QstopDiffDownStart == QstartDiffUpStop){ #they are equidistant
            newPositions = range(ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Stop],ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Start])
            return(getTopNCellTypes(ff_dt=ff_dt,chromoNum=chromoNum, range_v=newPositions, N=N, colStrWithChromLocations_str=colStrWithChromLocations_str,ffDtPath,ffOntIdInstead,optRank_Mat,optURLdecode,ffSorted_dt))
          }
        }
      }
    }
  }
}
getSortedDTwithSplitOfFirstCol <-
function(ff_dt,nameOfCol1="00Annotation"){
  #' Generate and return a data table of the fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #' 
  #' @param ff_dt a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user.
  #' @param nameOfCol1 a string matching the column in ff_dt that contains the TSS ID information (e.g., "chr10:10..20,-"). This function will parse out the chromosome, start, and stop.
  #' @return A data.table of fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' sorted_dt = getSortedDTwithSplitOfFirstCol(ff_dt = ffDummy_dt,nameOfCol1="Annotation")
  #' @export
  
  
  #Parse the fantom5 annotation (TSS identifier) column
  splits_ls = lapply(ff_dt[,get(nameOfCol1)], function(x) unlist(strsplit(x,':|\\.+|,'))[1:3])
  splits_mat = matrix(unlist(splits_ls), nrow=length(splits_ls), byrow = T)
  
  #Get a placeholder vector of column names of ff_dt to later arrange the column order. You'll have to get a deep copy of this because data.table does everything by reference
  cNamesTmp = copy(names(ff_dt))
  
  #Remove chr string and substitute X and Y chromosome for 23 and 24, respectively, temporarily while the rows get sorted.
  ff_dt[,chr := gsub('chr','',splits_mat[,1])]
  ff_dt[,chr := gsub('X','23',ff_dt[,chr])]
  ff_dt[,chr := gsub('Y','24',ff_dt[,chr])]
  
  #Populate start and stop columns with the parsed start and stop locations. Re-arrange the column order.
  ff_dt[,Start := splits_mat[,2]]
  ff_dt[,Stop := splits_mat[,3]]
  setcolorder(ff_dt,c("chr", "Start", "Stop",cNamesTmp))
  
  #coerce chr, Start, and Stop into numeric with the idea being that things that turn into NAs *should* be NAs (e.g., see the first few rows of ff_dt)
  ff_dt$chr = as.numeric(ff_dt$chr)
  ff_dt$Start = as.numeric(ff_dt$Start)
  ff_dt$Stop = as.numeric(ff_dt$Stop)
  
  #Sort ff_dt by chromosome and then by start position
  ff_dt_return = ff_dt[order(chr,Start)]
  
  #Restore X and Y chromosomes to their former glory (i.e., their previous chromosome identifiers)
  ff_dt_return[,chr := gsub('23','X',ff_dt_return[,chr])]
  ff_dt_return[,chr := gsub('24','Y',ff_dt_return[,chr])]
  return(ff_dt_return)
}
ConvertFFlsOflsOfDtToDt <-
function(listOfListOfDataTables_input,N){
  #' Generate a data.table with N+2 columns. The first column ("query_info") contains the chromosome, start, and stop positions of the query, separated by underscore delimiters. The second column ("match_info") contains the chromosome, start, and stop position of the INTERSECT region of the query position and any match within ff_dt. It also includes the original column name of its match in ff_dt following an "OriginalColName" string and an underscore delimiter. Failing a match, the match_info column will report the entire position of the nearest ff_dt TSS. If there are two nearest ff_dt TSSs, two separate matches will be reported, each reporting the entire position of its ff_dt TSS. Note that the query position may intersect more than one ff_dt TSS. In this case, each intersect will be reported as a separate row in the output. The output data.table will therefore have nrow() >= nrow(query_dt).
  #' 
  #' @param listOfListOfDataTables_input a list of list of data.tables generated by the repeated calls to the getTopNCellTypes function in the FFTopNCellTypesWrapper function (i.e., the FFtopNcells_lsOfdt_ls object from the wrapper function).
  #' @param N a numeric() argument the retreives the top N largest tpm-containing columns in a given row of ff_dt
  #' @return A list containing two elements: 1) a data table with ncol() N+2. The first column ("query_info") contains the chromosome, start, and stop positions of the query, separated by underscore delimiters. The second column ("match_info") contains the chromosome, start, and stop position of the INTERSECT region of the query position and any match within ff_dt. It also includes the original column name of its match in ff_dt following an "OriginalColName" string and an underscore delimiter. Failing a match, the match_info column will report the entire position of the nearest ff_dt TSS. If there are two nearest ff_dt TSSs, two separate matches will be reported, each reporting the entire position of its ff_dt TSS. Note that the query position may intersect more than one ff_dt TSS. In this case, each intersect will be reported as a separate row in the output. The output data.table will therefore have nrow() >= nrow(query_dt). 2) A list of list of data tables. The length() of the list will be nrow(query_dt). Each element of the list will contain a list of data tables, where each element of this list will represent one match (as previously mentioned, there will always be at least one "match"--direct or nearest--and sometimes more than one match). Each match will correspond to a data.table with nrow()=1 and ncol()=N. The names() of the data.table will correspond to the cell types/samples in ff_dt, and the entries in the table will be tpm counts from ff_dt, sorted from highest to lowest. If a chromosome number in the query data table doesn't match a chromosome number in ff_dt, it simply gets skipped.
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' q_dt = data.table(chrP = c(10,10,10,10,11), startP = c(1,31,1,10000,1), stopP=c(10000,33,3,10030,100))
  #' N = 2
  #' output_ls = FFTopNCellTypesWrapper(ff_dt = ffDummy_dt, query_dt = q_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 2, savePath = "~/Desktop/", colStrWithChromLocations_str="Annotation")
  #' listOfListOfDataTables_ls = output_ls[[2]]
  #' Converted_dt = ConvertFFlsOflsOfDtToDt(listOfListOfDataTables_ls,N)
  #' @export
  
  #Start with empty data.table
  master_dt = NULL
  
  #Fresh copy of listOfListOfDataTables_input to avoid manipulating it by reference
  listOfListOfDataTables = copy(listOfListOfDataTables_input)
  
  #Loop through each data table in the list of list of data tables
  for(m in 1:length(listOfListOfDataTables)){
    for(s in 1:length(listOfListOfDataTables[[m]])){
      
      #And create a "row" which has to be a data.table for rbind to work containing a query_info column
      contributing_dt = data.table(query_info = names(listOfListOfDataTables)[m] )
      
      #Add a match_info column to that same "row"
      contributing_dt[,match_info := names(listOfListOfDataTables[[m]][s])]
      
      for (tpmcol in 1:N){
        #And then add a new column for each of the N most-higlhy-ranked tpm counts, and populate with the CAGE peak IDs corresponding to these highly-ranked counts
        contributing_dt [ ,paste0("CellType_",tpmcol) := names(listOfListOfDataTables[[m]][[s]])[tpmcol]]
      }
      
      #Tack the "row" onto the master data.table
      master_dt = rbind(master_dt,contributing_dt)
      
      #And reset the "row"
      rm(contributing_dt)
    }
  }
  
  #When you're done looping through the list of lists of data.tables, return the resulting master data.table that's been generated
  return(master_dt)
}
FFTopNCellTypesWrapper <-
function(ff_dt, query_dt, ChrColName, StartColName, StopColName, N, savePath, colStrWithChromLocations_str="00Annotation", ffDtPath=NULL,optRankMatPath=NULL,ffOntIdInstead=F,optURLdecode=F,ffSorted_dt=NULL, optQueryIDcolName=NULL, optRankMat = NULL){ #, queryDtENSRidColName="ENSRid", queryDtRegDescriptionColName="gwasEnsemblRegDescription", queryDtVarLenAcrossENSRidColName="VarLenAcrossENSRid_x"){
  #' Wrap the output of getTopNCellTypes over a data.table of queries instead of one at a time and convert the output of that function into a data.table of cell type/sample names instead of the list of list of data.tables using the ConvertFFlsOflsOfDtToDt function.
  #' 
  #' @param ff_dt a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user. If no data.table is provided, the function will look for the hard-coded path in exacloud.
  #' @param query_dt a data.table with at least three columns: one corresponding to a chromosome number, one corresponding to the start position, and one correpsonding to the stop position. Rows populated with positions (e.g., regulatory element positions) of interest to the user.
  #' @param ChrColName a string containing the name() of the column in query_dt containing the chromosome position information.
  #' @param StartColName a string containing the name() of the column in query_dt containing the start position information.
  #' @param StopColName a string containing the name() of the column in query_dt containing the stop position information.
  #' @param N a numeric() argument the retreives the top N largest tpm-containing columns in a given row of ff_dt
  #' @param savePath a string that specifies the destination directory for the two return items (a data table and list of list of data tables)
  #' @param ffDtPath an optional path provided by the user to supply the fantom5 tpm data.
  #' @param optRankMatPath an optional path to a pre-existing ranked matrix. Importantly, the RankMat MUST be sorted by chromosome and then by start position!
  #' @param ffOntIdInstead a boolean designating whether the users wants to substitute the full URL names of the cell types/samples with just the fantom5 ontology IDs instead.
  #' @param optURLdecode a boolean designating whether the user wants to apply a URL-decoder function on the colnames of topCells_dt (i.e., on the column names of the sorted/ranked cell types/samples)
  #' @param ffSorted_dt an optional data.table containing ff_dt sorted by chromosome and then by start position. Most useful and time-saving if lots of calls to the sort() function expected (i.e., if lots of hits are expected or if this function is being called in a for loop as with FFTopNCellTypesWrapper).
  #' @param optQueryIDcolName an optional string containing the column name of an identifier column in query_dt. If the string is provided, it will be concatenated onto the query_info string
  #' @param optRank_Mat an optional matrix containing the ranked values of all of the fantom5 tpm entries. Importantly, the RankMat MUST be sorted by chromosome and then by start position! Most useful and time-saving if lots of calls to the sort() function expected (i.e., if lots of hits are expected or if this function is being called in a for loop as with FFTopNCellTypesWrapper).
  #' @return A list containing two elements: 1) a data table with ncol() N+2. The first column ("query_info") contains the chromosome, start, and stop positions of the query, separated by underscore delimiters. The second column ("match_info") contains the chromosome, start, and stop position of the INTERSECT region of the query position and any match within ff_dt. It also includes the original column name of its match in ff_dt following an "OriginalColName" string and an underscore delimiter. Failing a match, the match_info column will report the entire position of the nearest ff_dt TSS. If there are two nearest ff_dt TSSs, two separate matches will be reported, each reporting the entire position of its ff_dt TSS. Note that the query position may intersect more than one ff_dt TSS. In this case, each intersect will be reported as a separate row in the output. The output data.table will therefore have nrow() >= nrow(query_dt). 2) A list of list of data tables. The length() of the list will be nrow(query_dt). Each element of the list will contain a list of data tables, where each element of this list will represent one match (as previously mentioned, there will always be at least one "match"--direct or nearest--and sometimes more than one match). Each match will correspond to a data.table with nrow()=1 and ncol()=N. The names() of the data.table will correspond to the cell types/samples in ff_dt, and the entries in the table will be tpm counts from ff_dt, sorted from highest to lowest. If a chromosome number in the query data table doesn't match a chromosome number in ff_dt, it simply gets skipped.
  #' @examples
  #' require(data.table)
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0))
  #' q_dt = data.table(ensrID = c("ENSR001342", "ENSR001453", "ENSR00923412", "ENSR12399912", "ENSR00125632"), chrP = c(10,10,10,10,11), startP = c(1,31,1,10000,1), stopP=c(10000,33,3,10030,100))
  #' optQueryIDcolNameInst = "ensrID"
  #' output_ls = FFTopNCellTypesWrapper(ff_dt = ffDummy_dt, query_dt = q_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 2, savePath = "~/Desktop/", colStrWithChromLocations_str="Annotation", optQueryIDcolName=optQueryIDcolNameInst)
  #' output2_ls = FFTopNCellTypesWrapper(ff_dt = ffDummy_dt, query_dt = q_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 2, savePath = "~/Desktop/", colStrWithChromLocations_str="Annotation", optRankMatPath="/Users/fishema/Desktop/Projects/Bcore_clients/Klein_AMD/data/AllFFEntriesRankedSortedWithParsedCols.txt") #Attn test this
  #' optRank_dt = getSortedDTwithSplitOfFirstCol(copy(ffDummy_dt),"Annotation")
  #' numericPartOfFf_mat = data.matrix(optRank_dt[,grep('tpm', names(optRank_dt)),with=F])
  #' ffRanks_mat = matrix(nrow=(nrow(numericPartOfFf_mat)), ncol=length(grep('tpm', colnames(ffDummy_dt))))
  #' colnames(ffRanks_mat) = colnames(numericPartOfFf_mat)[grep('tpm', colnames(numericPartOfFf_mat))]
  #' for (r in 1:nrow(numericPartOfFf_mat)){
  #'   ffRanks_mat[r,] = rank(-numericPartOfFf_mat[r,], na.last = T)
  #' }
  #' optRank_Mat = ffRanks_mat
  #' optRank_Mat = rbind(colnames(optRank_Mat), optRank_Mat)
  #' optRank_Mat = cbind(c("Fantom5ID",ffDummy_dt[["Annotation"]]),optRank_Mat)
  #' write.table(optRank_Mat,file="~/Git_repos/AbundanceFunctions/testRankMat.txt", sep="\t",row.names = F, col.names = F)
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0))
  #' sortedFF_dt = getSortedDTwithSplitOfFirstCol(copy(ffDummy_dt),"Annotation")
  #' output3_ls = FFTopNCellTypesWrapper(ff_dt = ffDummy_dt, query_dt = q_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 2, savePath = "~/Desktop/", colStrWithChromLocations_str="Annotation", optRankMatPath="~/Git_repos/AbundanceFunctions/testRankMat.txt")
  #' output4_ls = FFTopNCellTypesWrapper(ff_dt = NULL, query_dt = q_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 20, ffDtPath="grep -v '^#'/Users/fishema/Downloads/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt",savePath = "~/Desktop/", colStrWithChromLocations_str="00Annotation", optRankMatPath="/Users/fishema/Desktop/Projects/Bcore_clients/Klein_AMD/data/AllFFEntriesRankedSortedWithParsedCols.txt")
  #' @export
  
  require(data.table)
  
  #If ff_dt was not provided, to load it from user-provided path
  if(is.null(ff_dt)){ 
    if(!is.null(ffDtPath)){
      ff_dt = fread(ffDtPath)
    }else{
      stop("No Fantom5 data path or matrix provided.")
    }
  }
  
  #If ranked matrix was provided for speed, load it
  if(!is.null(optRankMatPath)){
    Rank_mat = as.matrix(fread(optRankMatPath))
    optRank_MatInst = Rank_mat
  }else{
    if(!is.null(optRankMat)){
      optRank_MatInst = as.matrix(optRankMat) #Attn. what kind of object is this? Make sure this doesn't get coerced to a data.frame.
    }else{
      optRank_MatInst = NULL
    }
  }
  
  #Check that ff_dt and query_dt are data.tables. If not, try to coerce them to data.tables
  if(! (is.data.table(ff_dt) & is.data.table(query_dt))){
    ff_dt = as.data.table(ff_dt)
    query_dt = as.data.table(query_dt)
  }
  
  #Sort the ff_dt once before passing through getTopNCellTypes loop
  ffSortedInst_dt = getSortedDTwithSplitOfFirstCol(ff_dt=copy(ff_dt), nameOfCol1=colStrWithChromLocations_str)
  
  #Initiate and populate the list of list of data.tables
  FFtopNcells_lsOfdt_ls = list()
  ffDtPathInst = ffDtPath
  optRankMatPathInst = optRankMatPath
  ffOntIdInsteadInst = ffOntIdInstead
  optURLdecodeInst = optURLdecode
  for(r in 1:nrow(query_dt)){
    if(r%%50==0) {print(r)} #print r every 50th iteration of the loop to reassure an eager analyst
    FFtopNcells_lsOfdt_ls[[r]] = getTopNCellTypes(ff_dt=ff_dt,chromoNum=as.numeric(query_dt[r,get(ChrColName)]), range_v=c(as.numeric(query_dt[r,get(StartColName)]),as.numeric(query_dt[r,get(StopColName)])), N=N, colStrWithChromLocations_str=colStrWithChromLocations_str, ffDtPath=ffDtPathInst, ffOntIdInstead=ffOntIdInsteadInst, optRank_Mat=optRank_MatInst, optURLdecode=optURLdecodeInst, ffSorted_dt=ffSortedInst_dt)
    #............................getTopNCellTypes(ff_dt=ff_dt,chromoNum=chromoNum, ..............................range_v=newPositions, .........................................................................N=N, colStrWithChromLocations_str=colStrWithChromLocations_str,ffDtPath,...............ffOntIdInstead,....................optRank_Mat,.................optURLdecode,..................ffSorted_dt)
    #if the r-th element of the list was successfully populated (it wouldn't be if, say, there was no chromosome number in ff_dt that matched the query), let's go ahead and give that element a name
    if(length(FFtopNcells_lsOfdt_ls) ==r){
      if(!is.null(optQueryIDcolName)){
        names(FFtopNcells_lsOfdt_ls)[r] = paste(query_dt[r,get(ChrColName)],query_dt[r,get(StartColName)],query_dt[r,get(StopColName)],"query",query_dt[r,get(optQueryIDcolName)],sep="_")
      }else{
        names(FFtopNcells_lsOfdt_ls)[r] = paste(query_dt[r,get(ChrColName)],query_dt[r,get(StartColName)],query_dt[r,get(StopColName)],"query",sep="_") #Attn
      }
      
    }
  }
  
  #Save the list of list of data.tables to a text file. Not for client viewing, but useful for analyst to check.
  sink(file=paste0(savePath,"FFtop",N,"cells_lsOfdt_ls.txt"))
  print(FFtopNcells_lsOfdt_ls)
  sink()
  
  #Convert the list of list of data.tables to one data.table more amenable to downstream analyses
  Converted_dt=ConvertFFlsOflsOfDtToDt(FFtopNcells_lsOfdt_ls,N)
  write.table(Converted_dt, file=paste0(savePath, "Fantom5Top",N,"CellTypes.txt"), sep="\t", row.names = F)
  
  #Function will return both the data.table Converted_dt and the list of list of data.tables.
  output_ls = list(Converted_dt, FFtopNcells_lsOfdt_ls)
  names(output_ls) = c(paste0("Data table of query and match output for top ",N," cell types" ), "List of list of data tables output containing actual tpms")
  return(output_ls)
}
fetchTssRangesWithXFoldChangeInFFOntologyID <-
function(ffID_v, ff_mat, x, onlyOverExprs=T, TSSidColStr ="00Annotation" , EntrezIDcolStr="entrezgene_id"){
  #' For each fantom5 ontology ID or query string of interest in ffID_v, generates a matrix of just entries with log2-fold ratios greater than x (or abs(log2(ratio)>x)), decorated with chromosome, start position, stop position, entrezID, and description. If no entries match at the >x (or abs()>x) threshold, populate that matrix with NULL. Assembles these matrices into a list and returns the list.
  #' 
  #' @param ffID_v a vector containing keywords or -even better- fantom5 ontology IDs by which to narrow the search for 
  #' @param ff_mat a matrix version of the fantom5 tpm count data (hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt). Importantly, the 1st column must be the fantom5 TSS identifier containing positional information and there must be a column containing entrezgene_id.
  #' @param onlyOverExprs a boolean specifying whether the user wants to inspect only over-expressed
  #' @param x a number representing the desired log2-fold over- and possibly under-expression desired
  #' @param TSSidColStr a string containing the column name of the TSS ID column of ff_mat (default is "00Annotation")
  #' @param EntrezIDcolStr a string containing the column name of the entrez ID column of ff_mat (default is "entrezgene_id")
  #' @return a list of matrices, one matrix per list element, of just entries with log2-fold ratios greater than x (or abs(log2(ratio)>x)), decorated with chromosome, start position, stop position, entrezID, and description.
  #' 
  #' @examples
  #' ffDummy_mat = as.matrix(data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3), entrezID=c("ID:1","ID:2","ID:3","ID:4","ID:5","ID:6","ID7","ID:8","ID:9","ID:10","ID:11"),description=c("Description1","Description2","Description3","Description4","Description5","Description6","Description7","Description8","Description9","Description10","Description11")))
  #' ffID_v1 = c("Two", "Three")
  #' EntrezIDcolStr1= "entrezID"
  #' TSSidColStr1 = "Annotation"
  #' x_var = 3
  #' hits_ls = fetchTssRangesWithXFoldChangeInFFOntologyID(ffID_v=ffID_v1, ff_mat=ffDummy_mat, x=x_var, onlyOverExprs=T,TSSidColStr=TSSidColStr1,EntrezIDcolStr=EntrezIDcolStr1)
  #' @export
  
  lg2Ratio_mat = generateLg2RatioToAveAllMat(ff_mat)
  
  #search for those rows with >abs(x) fc for one of any columns containing ids in ffID_v and return those log ratios, along with the chr, start, and stop
  hit_ls = list()
  
  #loop through the fantom5 ontology IDs or query strings and if you find entries with log2-scale ratios >x, assemble a matrix of just these entries, decorated with chromosome, start position, stop position, entrezID, and description
  for (ffID in ffID_v){
    if(onlyOverExprs==TRUE){
      HighExprsInColOfInterestLogic_v = lg2Ratio_mat[,grep(ffID, colnames(lg2Ratio_mat))] > x
    }else{
      HighExprsInColOfInterestLogic_v = abs(lg2Ratio_mat[,grep(ffID, colnames(lg2Ratio_mat))]) > abs(x)
    }
    
    #If there were TSSs in our cell type/sample of interest >x
    if(sum(HighExprsInColOfInterestLogic_v, na.rm=T)>0){
      #create a matrix of that cell type sample column
      logRatio_mat = lg2Ratio_mat[HighExprsInColOfInterestLogic_v,grep(ffID, colnames(lg2Ratio_mat)),drop=F]
      
      #parse the identifier column of ff_mat
      splits_ls = lapply(ff_mat[HighExprsInColOfInterestLogic_v,TSSidColStr,drop=F], function(x) unlist(strsplit(x,':|\\.+|,'))[1:3])
      splits_mat = matrix(unlist(splits_ls), nrow=length(splits_ls), byrow = T)
      
      entrezIDsplit_ls = lapply(ff_mat[HighExprsInColOfInterestLogic_v, EntrezIDcolStr, drop=F], function(x) unlist(strsplit(x,':'))[2])
      entrezID_mat = matrix(unlist(entrezIDsplit_ls), nrow=length(entrezIDsplit_ls), byrow = T)
      description_mat = ff_mat[HighExprsInColOfInterestLogic_v, "description", drop=F]
      logRatioJoin_mat = cbind(splits_mat,entrezID_mat,description_mat, logRatio_mat)
      colnames(logRatioJoin_mat) = c("Chr","Start","Stop","EntrezID","Description", colnames(logRatio_mat))
      hit_ls[[length(hit_ls)+1]] = logRatioJoin_mat=logRatioJoin_mat
      names(hit_ls)[length(hit_ls)] = colnames(lg2Ratio_mat)[grep(ffID, colnames(lg2Ratio_mat))]
    }
  }
  return(hit_ls)
}
generateLg2RatioToAveAllMat <-
function (ff_DtOrMat){
  #' Generates and returns a matrix containing log2-scale ratios to average all for the numerical portion of the fantom5 data
  #' 
  #' @param ff_DtOrMat a data.table, matrix, or data.frame of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user.
  #' @return A matrix of log2-scale ratios to average all for the numerical portion of the fantom5 data. Zeroes in the original tpm count matrix are replaced with the smallest non-zero value within the same CAGE peak ID to avoid log2(0)-related errors/issues.
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' lgRatio_mat = generateLg2RatioToAveAllMat(ffDummy_dt)
  #' @export
  
  #Coerce ff_mat to matrix if not already and possible
  ff_mat = as.matrix(ff_DtOrMat)
  
  #Subset ff_mat by the portion that is numeric (the tpm columns), coerce to numeric, check for 0s,NAs???? and replace with col-wise non-zero min.
  numericPartOfFf_mat = data.matrix(ff_mat[,grep('tpm', colnames(ff_mat))])
  numericPartOfFf_mat = apply(numericPartOfFf_mat, 2, function(x) as.numeric(x))
  numericPartOfFf_mat_test = apply(numericPartOfFf_mat, 2, function(x) replace(x,x %in% c(0), min(x[x>0]))) #include NA?
  #min(numericPartOfFf_mat[numericPartOfFf_mat > 0],na.rm=T)
  
  #log2 transform the tpm count containing columns
  lg2FF_mat = log2(numericPartOfFf_mat_test)
  
  #take rowmeans of the same columns
  rMeans = rowMeans(lg2FF_mat)
  
  #find ratio to average for each of the columns (a difference in log2 scale)
  lg2Ratio_mat= lg2FF_mat - rMeans
  
  return(lg2Ratio_mat)
}
