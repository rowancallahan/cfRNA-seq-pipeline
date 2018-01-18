#!/bin/bash
# These values need to be set in order for the workflow to function.

# Initialize BioCoders
BIOCODERS="/home/exacloud/lustre1/BioCoders"
. $BIOCODERS/USE_APPS --python2

# Java and resource allocation
JAVA8="/opt/installed/jdk1.8.0_111/bin/java"
JAVAVM="-Xms16g -Xmx54g"
CTHREADS=1
DTHREADS=24

# Resource files
GENOME_FASTA="/opt/installed/galaxy_genomes/hg19/Homo_sapiens_assembly19.fasta"
INDEX_DIR="/opt/installed/galaxy_genomes/hg19/star"
GTF="$BIOCODERS/DataResources/Transcriptomes/Human/grch37_13/ref_GRCh37.p13_top_level.gff3"
DBSNP="$BIOCODERS/DataResources/Variation/PopulationDatabases/dbSNP/146/grch37_p13_146.vcf"

# Jars
GATK="$BIOCODERS/Applications/GATK.jar"
PICARD="$BIOCODERS/Applications/picard.jar"

# Set regex for file searches based on sample.
READ1_REGEX='_[rR]1[._]'
READ2_REGEX='_[rR]2[._]'
FILES=($PROJ/raw/*)

# Temp
READLEN="75"
