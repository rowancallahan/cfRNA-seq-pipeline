#!/bin/bash

# RSEM v1.2.21

export PATH="/opt/installed:$PATH"
PROJ="/home/users/letaw/lustre1/projs/beataml2"
RAW="/home/users/letaw/lustre1/projs/beataml/raw"
THREADS=14

echo "Concatenating forward and reverse FASTQ's."

cat $RAW/$1/*R1* >> $RAW/$1\_R1.fastq.gz
cat $RAW/$1/*R2* >> $RAW/$1\_R2.fastq.gz

echo "Running RSEM."

/mnt/lustre1/CompBio/bin/rsem-calculate-expression -p $THREADS --paired-end <(zcat $RAW/$1\_R1.fastq.gz) <(zcat $RAW/$1\_R2.fastq.gz) $PROJ/hg19 $PROJ/output/$1
#<(zcat mate1.fastq.gz) <(zcat mate2.fastq.gz)
