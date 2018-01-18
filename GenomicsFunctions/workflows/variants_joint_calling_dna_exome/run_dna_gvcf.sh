#!/bin/bash

# USAGE: bash run_dna_gvcf.sh [1|2] [SAMPLE|VCF_NAME]
# BRIEF INSTRUCTIONS:
# Make sure your PROJ folder is set below.
# Make sure the proper region definiton file is being pointed to by PROBES (dna_gvcf_source.sh).
# Use option #1, along with sample name (no extension) to run each sample through the gVCF creation pipeline.
# After these are done, you can use option #2, along with the name of your output combined VCF.
# These are meant to be sent out via HTCondor, so the submit script will have to be filled as well.
# Other cleanup options available below.

# Requires: parse_read_groups.sh (GenomicsFunctions/fastq_utilities)
# Requires: dna_gvcf_source.sh, dna_gvcf.sh, dna_gvcf_genotype.sh (GenomicsFunctions/workflows/variants_joint_calling_dna_exome)
# CODED BY: John H. Letaw (letaw@ohsu.edu)
# TODO: getopts

VERSION="0.3.0"

# This needs to be set!                                                                                                        
PROJ="/home/exacloud/lustre1/BioCoders/ProjectCollaborations/sjkim_rop/analysis"

if [ $1 == "1" ]; then
    if [ $# == 2 ]; then
	source $PROJ/scripts/dna_gvcf_source.sh $2 "VCF_NAME"
	source $PROJ/scripts/parse_read_groups.sh
	source $PROJ/scripts/dna_gvcf.sh
    else
 	echo "Incorrect number of arguments, must have SAMPLE."
	exit
    fi
elif [ $1 == "2" ]; then
    if [ $# == 2 ]; then
	source $PROJ/scripts/dna_gvcf_source.sh "SAMPLE" $2
	source $PROJ/scripts/dna_gvcf_genotype.sh
    else
 	echo "Incorrect number of arguments, must have VCF_NAME."
	exit
    fi
# Remove directory structure, not including temp directory.
elif [ $1 == "99999" ]; then
    rm -rf $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/br $PROJ/pr $PROJ/hc $PROJ/genotype_gvcfs $PROJ/vqsr
elif [ $1 == "rmtemp" ]; then
    rm -rf $PROJ/temp
elif [ $1 == "rmlogs" ]; then
    rm -rf $PROJ/logs
else
    echo "Incorrect first argument, must be either 1 or 2."
    exit
fi
