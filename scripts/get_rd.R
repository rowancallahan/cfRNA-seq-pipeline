files= snakemake@input

for(i in 1:length(files)){
  samp_counts = read.table(files[[i]], stringsAsFactors  = F,sep="\t")
  samp=gsub(".read_distribution.txt","",basename(files[[i]]))
  total_assigned_tags <- as.numeric(unlist(strsplit(samp_counts$V1[3], " "))[14])
  CDS_tags <- as.numeric(unlist(strsplit(samp_counts$V1[6], " "))[24])
  UTR5_tags <- as.numeric(unlist(strsplit(samp_counts$V1[7], " "))[23])
  UTR3_tags <- as.numeric(unlist(strsplit(samp_counts$V1[7], " "))[71])
  Intron_tags <- as.numeric(unlist(strsplit(samp_counts$V1[8], " "))[24])

  Exon_Fraction <- (CDS_tags + UTR5_tags + UTR3_tags)/total_assigned_tags
  Intron_Fraction <- Intron_tags/total_assigned_tags
  samp_data <- as.data.frame(rbind(Exon_Fraction, Intron_Fraction))
  names(samp_data) <- samp
             
  if(i==1) {
    all_data <- samp_data
   } else {
    all_data <- cbind(all_data, samp_data) 
  }  
}   
write.table(all_data,file=snakemake@output[[1]],sep="\t",quote=F)
