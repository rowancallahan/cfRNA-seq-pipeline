#!/bin/bash
### RNA-seq DE pipeline script. (GSNAP based, adjusted to run STAR)
# modified to run either single-end or paired-end RNAseq
###
# Set directories at top of script.  If you are unsure about using this by pushing the enter button, perform  each step individually.
# Program and genome directories should be adjusted for current installations 
# Raw FASTQ files should be one per sample and placed within the raw subdirectory of your project folder.
###
# John Letaw 12/31/14 and Julja Burchard
# letaw@ohsu.edu
#
###
if [ "$#" -ne 1 ] || [ $1 == "-h" ] 
then
  echo "Usage: $0 nth-job:user:projectdir:listfile:paired:stranded"
  exit 1
fi

### !!set me!! ###
# base of paths
ROOT_DIR="/home/exacloud/lustre1/CompBio"
### The reference genome and gene annotation.
####!!!>>>CHECK THIS<<<!!!####
genome="hg38"
REF_FA="$ROOT_DIR/genomic_resources/genomes/$genome/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
GTF="/home/exacloud/lustre1/CompBio/genomic_resources/gtf/Homo_sapiens.GRCh38.84.gtf"
# file parameters  #
fastq="fastq.gz"   # file ending for FASTQ files
raw="data/fastq"   # location for FASTQ files
out="results/STAR" # location for output directories

# constants
BIN_DIR="$ROOT_DIR/bin"
### Where the genome FASTAs are.
INDEX_DIR="$ROOT_DIR/genomic_resources/genomes/$genome/star"
### Location of samtools
SAMTOOLS="/home/exacloud/lustre1/CompBio/bin/samtools"
### Location of PICARD
PICARD="/opt/installed/picard/picard-tools-1.110"
### Prefix used for index files, and splicesites files.
GENOME="$genome"
### How many processing threads?
THREADS=15
### for Picard, if run
RGPL="illumina"

# arguments
var=$1
eval `echo $var| sed -e 's/^\([^:]\{1,\}\)\:\([^:]\{1,\}\)\:\([^:]\{1,\}\)\:\([^:]\{1,\}\):\([^:]\{1,\}\):\([^:]\{1,\}\)$/n=\1 u=\2 mydir=\3 fil=\4 pr=\5 strand=\6/'`
myfile="$mydir/$fil" #listfile of samples, one per line
mynum=`expr $n + 1`

PROJECT_DIR=$mydir
SAMP_NAME=`head -$mynum $myfile|tail -1`
PROJ_NAME=`echo $mydir| sed -e 's/^.\{1,\}\/\([^\/]\{1,\}\)$/\1/'`
READLEN=`expr $(/bin/gzip -cd $PROJECT_DIR/$raw/$SAMP_NAME.$fastq| head -2|tail -1| wc -c) - 1` # minus one takes off 1 char for wc -c counting the newline
###  is this a paired-end or single-end read library? It could be yes or no
PAIRED=$pr
###  is this a strand-specific assay?  It could be yes, no, or reverse.  See HTSeq-count help for more.
STRANDED=$strand
### for Picard, if run
RGDS=$PROJ_NAME

# set up environment
if [ -e /home/users/$u/.bash_profile ]
then
    . /home/users/$u/.bash_profile
else
  if [ -e /home/users/$u/.profile ]
  then
    . /home/users/$u/.profile
  fi
fi

### Unset variables FWD and REV in case they are set in environment.
unset REV
unset FWD

### Set nullglob to avoid globs with zero elements appearing to be one element.
shopt -s nullglob

### go home
cd $PROJECT_DIR



### Clean old directory structure if 666 is specified on command line.
if [ $SAMP_NAME == 666 ]
then
    rm -r $PROJECT_DIR/$out $PROJECT_DIR/$out\_rg $PROJECT_DIR/$out\_dups $PROJECT_DIR/$out\_htseq
    exit 1
fi

### Create project directory structure if not present.

if [ ! -d $PROJECT_DIR/$out ]
then
  mkdir -p $PROJECT_DIR/$out 
fi
# being careful not to incite warning with shared directories, tho mkdir won't overwrite
if [ ! -d $PROJECT_DIR/logs ]
then
  mkdir -p $PROJECT_DIR/logs
fi

### Create STAR's genome index files if not present
# only saw 2 threads in use when 15 were requested, so limiting ask to 2 threads
# however there was a file issue so a retest is in order
if [ ! -d $INDEX_DIR ]
then
    mkdir $INDEX_DIR
    cd $INDEX_DIR # STAR writes to wd
    STAR --runMode genomeGenerate --genomeDir $INDEX_DIR --genomeFastaFiles $REF_FA --sjdbGTFfile $GTF --sjdbOverhang `expr $READLEN - 1` --runThreadN 2
fi

### set up filenames for paired or unpaired reads
# variable will contain space-delimited list if more than one file matches
#  exactly 2 files per sample are recommended for paired-end reads
# STAR takes space-separated pairs of fastq files for paired-end reads; split files are comma-delimited
# fastqc takes a space-delimited list of files to process
SAMP_FASTQ=`ls $PROJECT_DIR/$raw/$SAMP_NAME*$fastq`
	
### Run STAR
if [ ! -d $PROJECT_DIR/$out/$SAMP_NAME ]
then
    mkdir $PROJECT_DIR/$out/$SAMP_NAME
    cd $PROJECT_DIR/$out/$SAMP_NAME #STAR writes to wd

  if [ $STRANDED == "no" ]
  then
      oSsF="intronMotif"; oFIM="RemoveNoncanonicalUnannotated"
  else
      oSsF="None"; oFIM="None"
  fi
  # consider adding wiggle-track output? raw or norm, read 5' or all, stranded.. STAR can also generate wiggle or bedgraph from position-sorted BAM output
    STAR --genomeDir $INDEX_DIR --readFilesCommand zcat --outFileNamePrefix $SAMP_NAME\_ --readFilesIn $SAMP_FASTQ --runThreadN $THREADS --twopassMode Basic --quantMode GeneCounts --outReadsUnmapped Fastq --outSAMstrandField $oSsF --outFilterIntronMotifs $oFIM --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 56000000000 --bamRemoveDuplicatesType UniqueIdentical

fi

### Prepare for QC by modifying bam file

### go home
cd $PROJECT_DIR

### Add Read group.
# STAR sorted and handled duplicates, but RNASeQC expects readgroup marking 
# samtools and GATK tools also expect indexing
# the lines commented out are not set up for sorted BAM input or choices between paired and unpaired reads

RGID=$(/bin/gzip -cd $PROJECT_DIR/$raw/$SAMP_NAME.$fastq | head -1 | cut -d ':' -f 1-4);RGID=`echo "'"$RGID"'"` # protecting SRR whitespace from shell
RGPU=$(/bin/gzip -cd $PROJECT_DIR/$raw/$SAMP_NAME.$fastq | head -1 | cut -d ':' -f 1-3);RGPU=`echo "'"$RGPU"'"` # protecting SRR whitespace from shell
RGLB=$(/bin/gzip -cd $PROJECT_DIR/$raw/$SAMP_NAME.$fastq | head -1 | cut -d ':' -f 3-4)
RGSM=$RGPU # using run information instead of sample name

java -jar -Xms76G -Xmx76G $PICARD/AddOrReplaceReadGroups.jar INPUT=$PROJECT_DIR/$out/$SAMP_NAME/$SAMP_NAME\_Aligned.out.bam OUTPUT=$PROJECT_DIR/$out\_rg/$SAMP_NAME\_rg.bam CREATE_INDEX="true" RGID="$RGID" RGLB="$RGLB" RGPL=$RGPL RGPU="$RGPU" RGSM="$RGSM" RGDS=$RGDS MAX_RECORDS_IN_RAM=10000000


### Mark Duplicates

#java -jar -Xms76G -Xmx76G $PICARD/MarkDuplicates.jar INPUT=$PROJECT_DIR/$out\_rg/$SAMP_NAME\_rg.sam OUTPUT=$PROJECT_DIR/$out\_dups/$SAMP_NAME\_rg_dups.sam M=$PROJECT_DIR/$out\_dups/$SAMP_NAME\_rg_dups.metrics CREATE_INDEX=TRUE MAX_FILE_HANDLES=25000 MAX_RECORDS_IN_RAM=10000000


### Count gene set/alignment overlap.

#htseq-count -s $STRANDED -f bam -r pos -i gene_id $PROJECT_DIR/$out/$SAMP_NAME/$SAMP_NAME\_Aligned.sortedByCoord.out.bam $GTF > $PROJECT_DIR/$out/$SAMP_NAME\_htseq.counts


