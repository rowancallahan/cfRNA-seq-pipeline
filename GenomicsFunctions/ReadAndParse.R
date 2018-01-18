#####################################################################################################
# Set of functions useful for reading and parsing file formats
#
# readENSgtf             --- Read ensembl gtf format
#
# Author: Julja Burchard
# Started - 2016
#####################################################################################################

readENSgtf <-
  function (filename, gtf.colnames=c('seqname','source','feature','start','end','score','strand','frame','attribute'), feature.col=9, comma.sub="|", curly.open.sub='<',curly.end.sub='>',comment.char='#'){
    # arguments
    # filename = required pathname to file to be read
    # gtf.colnames = ENSEMBL GTF format field names. please check for accuracy
    
    ## imports
    require(data.table)
    require(jsonlite)
    
    ## constants
    # int for indexing
    a = as.integer(feature.col)
    
    ## read in data
    # default fread settings skip tracklines, identify typeof each column
    gtf = fread(filename)
    names(gtf) = gtf.colnames
    
    ## parse attribute column: transform to JSON and use JSON parsers
    
    # first, protect commas or curly braces if any within fields
    # commas
    mymk = grepl(',',gtf[,get(gtf.colnames[9])])
    if( sum(mymk)>0 ){
      set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub(',',comma.sub,gtf[mymk,get(gtf.colnames[a])]))
    }
    # opening curlies
    mymk = grepl('\\{',gtf[,get(gtf.colnames[9])])
    if( sum(mymk)>0 ){
      set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub('\\{',curly.open.sub,gtf[mymk,get(gtf.colnames[a])]))
    }
    # closing curlies
    mymk = grepl('\\}',gtf[,get(gtf.colnames[9])])
    if( sum(mymk)>0 ){
      set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub('\\}',curly.end.sub,gtf[mymk,get(gtf.colnames[a])]))
    }
    
    # next, clear non-JSON formatting
    # remove any comments
    mymk = grepl('[^\'"]#[^"\'{},]+',gtf[,get(gtf.colnames[9])])
    if( sum(mymk)>0 ){
      set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub(',?#[^"\'{},]+','',gtf[mymk,get(gtf.colnames[a])]))
    }
    # quote unquoted strings as found in Ensembl GFF3
    mymk = grepl('=[^\'"]\\w+[^"\']',gtf[,get(gtf.colnames[9])])
    if( sum(mymk)>0 ){
      set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub('(=[^\'"]\\w+[^"\'])','"\\1"',gtf[mymk,get(gtf.colnames[a])]))
    }
    
    # lastly, adapt GTF/GFF collapsed fields to JSON format
    set(x=gtf,j=a,value=paste("{",gtf[,get(gtf.colnames[a])],"}",sep=''))
    set(x=gtf,j=a,value=gsub('; ?',',',gtf[,get(gtf.colnames[a])]))
    set(x=gtf,j=a,value=gsub('([,{}])([A-Za-z0-9_.-]+) ','\\1"\\2" : ',gtf[,get(gtf.colnames[a])]))
    set(x=gtf,j=a,value=gsub(',}','},',gtf[,get(gtf.colnames[a])]))
    # begin and end properly
    gtf[1,(gtf.colnames[a]) := gsub('^','[',get(gtf.colnames[a]))]
    gtf[nrow(gtf),(gtf.colnames[a]) := gsub(',$',']',get(gtf.colnames[a]))]
    
    ## read JSON
    # JSON
    gtf.attributes = as.data.table(fromJSON(gtf[,get(gtf.colnames[a])]))
    # combine tables
    gtf = cbind(gtf,gtf.attributes)
    
    return( gtf ) 
  }
