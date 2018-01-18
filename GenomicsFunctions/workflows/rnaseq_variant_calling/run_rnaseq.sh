#!/bin/bash

# Start script.

WORKFLOW="RSMSVC"
VERSION="0.3.0"

# Create your project directory with one directory called 'raw' that contains
# input FASTQ files to be processed.  Also, a 'logs' directory should be included
# before invoking this script.

## Option --outSAMattrRGline deals with adding read groups.
## Check duplicate marking, doesn't look quite right
## Try to name filename something other than Aligned.
## Probably don't care about gene counts.
## --generate_md5 if we want to create more md5's and 
## check after each step.

usage() { echo "Usage: $0 [-p <project directory>] [-s <input fastq>]" 1>&2; exit 1; }

while getopts p:s:h FLAG; do
    case "${FLAG}" in
	h)
	    usage;;
        p)
	    PROJ=${OPTARG};;
	s)
	    SAMPLE=${OPTARG}
    esac
done

# FILENAME=$PROJ/raw/$SAMPLE;;

if [ ! -n "$PROJ" ]; then
    echo "Project directory must be specified!"
    exit 0
fi

# Load environmental variables from config file.
source $PROJ/scripts/rvc_config.sh
source $PROJ/scripts/parse_read_groups.sh

# check for GATK indexing of genome FASTA file
if [ ! -e $GENOME_FASTA\.fai ]; then
    DICTFILE=`echo $GENOME_FASTA|sed -e 's/fa$/dict/'`
    $JAVA8 $JAVAVM -jar $PICARD CreateSequenceDictionary R=$GENOME_FASTA O=$DICTFILE
    samtools faidx $GENOME_FASTA
fi

# # check for indexing of dbSNP VCF file
# if [ ! -e $DBSNP\.tbi ]; then
#     tabix -p vcf $DBSNP
# fi

# Create STAR's genome index files if not present                                                                               
# only saw 2 threads in use when 15 were requested, so limiting ask to 2 threads                                                
# however there was a file issue so a retest is in order                                                                        
if [ ! -d $INDEX_DIR ]
then
  mkdir $INDEX_DIR
  cd $INDEX_DIR # STAR writes to wd
  STAR --runMode genomeGenerate --genomeDir $INDEX_DIR --genomeFastaFiles $GENOME_FASTA --sjdbGTFfile $GTF --sjdbOverhang $READLEN --runThreadN 2
fi

# Create PROJ subdirectories.                                                                                                   
mkdir -p $PROJ/output
mkdir -p $PROJ/output/$SAMPLE

source $PROJ/scripts/rnaseq_multi_sample_variant_calling.sh
