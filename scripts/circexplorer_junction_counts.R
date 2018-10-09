samples = snakemake@config[['samples']]

files= snakemake@input

for(i in 1:length(files)){
  samp_counts = read.table(files[[i]], header = F, stringsAsFactors  = F,sep="\t")
  samp_counts=samp_counts[,c(1,2,3,13)]
  colnames(samp_counts)=c("chr","start","end","NumReads")
  samp_counts$start=samp_counts$start+1
  samp_counts$circRNA_ID=paste(samp_counts$chr,samp_counts$start,sep=":")
  samp_counts$circRNA_ID=paste(samp_counts$circRNA_ID,samp_counts$end,sep="|")
  samp_counts$NumReads=as.integer(samp_counts$NumReads)
  samp_counts=  samp_counts[,c(-1,-2,-3)]
  print(files[[i]])
  samp=gsub("circularRNA_known.txt","",files[[i]])
  colnames(samp_counts)=c(samp,"circRNA_ID")
             
  if(i==1) {
    all_counts=samp_counts
   }
    
  else
  {all_counts=merge(all_counts, samp_counts,by="circRNA_ID",all=T) }
    
}   
write.table(all_counts,file=snakemake@output[[1]],sep="\t",quote=F,row.names=F)
