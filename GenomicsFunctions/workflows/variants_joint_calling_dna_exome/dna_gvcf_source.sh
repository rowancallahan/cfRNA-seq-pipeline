#!/bin/bash

# NOTE: Most of these do not need to change if you are planning on working within BioCoders.
# PROBES does need to point to the correct region definiton file!
# Feel free to install this in a desriptive directory within BioCoders/DataResources.

# This script is meant to function within the BioCoders environment.
BIOCODERS="/home/exacloud/lustre1/BioCoders"
. $BIOCODERS/USE_APPS --python2

# Inputs
SAMPLE=$1
VCF_NAME=$2

# Make the $RAW and $LOGS directory before invoking this script.  Since you will most likely want to send
# this out via HTCondor, you must already have a logs directory in place.

# Specify number of worker threads.
DTHREADS=24
CTHREADS=4
MAX_READS=14000000

# Program directories
BIN_DIR="$BIOCODERS/Applications"
PICARD="$BIN_DIR/picard-tools-1.110"
GATK="$BIN_DIR/GATK.jar"
JAVA8="/opt/installed/jdk1.8.0_111/bin/java"

# Project directories
# Fill this with you project directory.  Subdirectories will be created here!
RESOURCE="$BIOCODERS/DataResources/Variation/PopulationDatabases"
HG="/opt/installed/galaxy_genomes/hg19/Homo_sapiens_assembly19.fasta"

# THIS NEEDS TO BE POINTING TO YOUR EXOME REGION FILE!!!
PROBES="$BIOCODERS/DataResources/IntervalDefinitions/Agilent_CRE/agilent_cre.interval_list"

# Create these before running the script.
RAW="$PROJ/raw"
LOGS="$PROJ/logs"

# Resource Files
DBSNP="$RESOURCE/dbSNP/146/grch37_p13_146.vcf"
MILLS="$RESOURCE/Mills/mills.vcf"
TG="$RESOURCE/1000G/1000g_phase1_broad_variant.vcf"
HAPMAP="$RESOURCE/Hapmap/hapmap_3.3.vcf"
OMNI="$RESOURCE/Omni/1000G_omni2.5.hg19.sites_nochr.vcf"

# Create directory structure.
mkdir -p $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/br $PROJ/pr $PROJ/hc $PROJ/temp $PROJ/genotype_gvcfs $PROJ/vqsr

# Copy resource files to temporary sample directories.  This helps prevent multiple processes from accessing the same files and causing problems.
# Probably just need to provide these files on each node, especially if pipelines will be used by general public.
mkdir -p $PROJ/temp/$SAMPLE
TEMP="$PROJ/temp/$SAMPLE"

# Rework this to not have to include filenames.
rsync -iv -a $HAPMAP $TEMP
rsync -iv -a $HAPMAP.idx $TEMP
HAPMAP="$TEMP/hapmap_3.3.vcf"

rsync -iv -a $OMNI $TEMP
rsync -iv -a $OMNI.idx $TEMP
OMNI="$TEMP/1000G_omni2.5.hg19.sites_nochr.vcf"

rsync -iv -a $MILLS $TEMP
rsync -iv -a $MILLS.idx $TEMP
MILLS="$TEMP/mills.vcf"

rsync -iv -a $DBSNP $TEMP
rsync -iv -a $DBSNP.idx $TEMP
DBSNP="$TEMP/grch37_p13_146.vcf"

rsync -iv -a $TG $TEMP
rsync -iv -a $TG.idx $TEMP
TG="$TEMP/1000g_phase1_broad_variant.vcf"
