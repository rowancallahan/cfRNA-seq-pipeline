#!/usr/bin/env Rscript
############################################################
# HOW TO RUN
# Rscript --vanilla run_fusionannotation.R starfusion_output path_to_AnnotationFunctions_dir output_name
# Rscript --vanilla run_fusionannotation.R star-fusion.fusion_candidates.final.abridged.FFPM /home/exacloud/lustre1/users/patterja/Functions/git fusions_annotated.txt
# Rscript --vanilla run_fusionannotation.R /Users/patterja/Workspace/FusionAnnotation/star-fusion.fusion_candidates.final.abridged.FFPM /Users/patterja/Workspace/FusionAnnotation fusions_annotated.txt
############################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ParseArguments #STAR-Fusion
args = commandArgs(trailingOnly=TRUE)
#R file that takes in oncofuse.txt and


starfusion = args[1]
AnnotationFunction_repo = args[2]
#starfusion="/Users/patterja/Workspace/FusionAnnotation/star-fusion.fusion_candidates.final.abridged.FFPM"
#AnnotationFunction_repo="/Users/patterja/Workspace/FusionAnnotation/"
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "fusions_annotated.txt"
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load Libraries

library(data.table)
source(paste0(AnnotationFunction_repo,"/fusionannotation.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#location of DataResources
resource_vardir="/Users/patterja/Workspace/FusionAnnotation/"
#resource_vardir="/home/users/patterja/BioCoders/DataResources/Variation/"

fusion_bedpe_dt=starfusion_to_bedpe(starfusion)

actionvarsPath=paste0(resource_vardir,"OncoKB/030117/allActionableVariants.txt")
annovarsPath=paste0(resource_vardir,"OncoKB/030117/allAnnotatedVariants.txt")

#see combineOncoKBrefs function
allvars_oncokb_dt= combineOncoKBrefs(annovarsPath = annovarsPath, actionvarsPath = actionvarsPath)

#annotate with OncoKB, see annoOncoKB_fusions function
fusion_bedpe_annoOncokb_dt=annoOncoKB_fusions(allvars_oncokb_dt = allvars_oncokb_dt, fusion_bedpe_dt = fusion_bedpe_dt)

fusion_bedpe_annoCivic_dt=annoCivic_fusions(fusion_bedpe_dt=fusion_bedpe_annoOncokb_dt, resource_vardir = resource_vardir)

#write.table(, file=args[2], row.names=FALSE)
write.table(fusion_bedpe_annoCivic_dt, file=args[3], row.names=FALSE, quote=FALSE, sep="\t")
  