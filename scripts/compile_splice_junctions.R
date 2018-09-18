
project_id = snakemake@params[['project_id']]

path="/results/"
files=dir(path,pattern="_circexp")

for(i in 1:length(files)){
  samp_counts = read.table(paste(path,files[i],"/circularRNA_known.txt",sep=""), header = F, stringsAsFactors  = F,sep="\t")
  samp_counts=samp_counts[,c(1,2,3,13)]
  colnames(samp_counts)=c("chr","start","end","NumReads")
  samp_counts$start=samp_counts$start+1
  samp_counts$circRNA_ID=paste(samp_counts$chr,samp_counts$start,sep=":")
  samp_counts$circRNA_ID=paste(samp_counts$circRNA_ID,samp_counts$end,sep="|")
  samp_counts$NumReads=as.integer(samp_counts$NumReads)
  samp_counts=  samp_counts[,c(-1,-2,-3)]
  print(files[i])
  samp=gsub("circularRNA_known.txt","",files[i])
  colnames(samp_counts)=c(samp,"circRNA_ID")
             
  if(i==1) {
    all_counts=samp_counts
   }
    
  else
  {all_counts=merge(all_counts, samp_counts,by="circRNA_ID",all=T) }
    
}   
write.table(all_counts,file=paste(project_id,'circexplorer_junctioncounts.txt',sep='_'),sep="\t",quote=F,row.names=F)
_round_sequencing/PP_cohort_sample_info_log_JE.2018.06.06.txt
