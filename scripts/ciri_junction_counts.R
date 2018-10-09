files=snakemake@input

for(i in 1:length(files)){
  samp_counts = read.table(files[[i]], header = F, stringsAsFactors  = F,sep="\t",skip=1)
  samp_counts=samp_counts[,c(1,5)]
  colnames(samp_counts)=c("circRNA_ID","NumReads")
  samp_counts$NumReads=as.integer(samp_counts$NumReads)
  print(files[[i]])
  print(dim(samp_counts))
  samp=gsub("_ciriout.txt","",files[[i]])
  colnames(samp_counts)=c("circRNA_ID",samp)
  if(i==1){
    all_counts=samp_counts
  } 
  else
  {all_counts=merge(all_counts, samp_counts,by="circRNA_ID",all=T) }

}

all_counts$counts_circ=apply(all_counts,1,function(x) {length(which(!is.na(x)))-1})
circ_counts=all_counts[,c("circRNA_ID","counts_circ")]
all_counts=all_counts[,!colnames(all_counts) %in% c("counts_circ")]
write.table(all_counts,file=snakemake@output[[1]],sep="\t",quote=F,row.names=F)
write.table(circ_counts,file=snakemake@output[[2]],sep="\t",quote=F,row.names=F)
