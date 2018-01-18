#!/bin/bash

### Specify number of worker threads.
DTHREADS=20
CTHREADS=8
MAX_READS=14000000

### Program directories
BIN_DIR="/opt/installed"
PICARD="$BIN_DIR/picard/picard-tools-1.110"
#GATK="$BIN_DIR/GATK/GenomeAnalysisTK-3.1.jar"
GATK="/mnt/lustre1/CompBio/bin/GATK-3.5.jar"

### Project directories
PROJ="/mnt/lustre1/KleinAMD/letaw"
RESOURCE="/mnt/lustre1/CompBio/genomic_resources"
HG_OLD="$RESOURCE/genomes/hg19/Homo_sapiens_assembly19.fasta"
HG="$BIN_DIR/galaxy_genomes/hg19/Homo_sapiens_assembly19.fasta"
RAW="$PROJ/raw_bam"
LOGS="$PROJ/logs"

### Resource Files
DBSNP="$RESOURCE/variation/dbsnp/142/grch37_p13_142.vcf"
MILLS="$RESOURCE/variation/mills/mills.vcf"
TG="$RESOURCE/variation/1000g/1000g_phase1_broad_variant.vcf"
COSMIC="/mnt/lustre1/CompBio/genomic_resources/variation/cosmic/cosmic_75_coding_hg19_sorted.vcf"

SAMPLE=$1

if [ $SAMPLE == "99999" ]; then
    rm -rf $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/ir $PROJ/br $PROJ/pr $PROJ/hc $PROJ/temp
    rm $LOGS/*
    exit
fi

### Create directory structure.
mkdir -p $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/ir $PROJ/br $PROJ/pr $PROJ/hc $PROJ/temp $PROJ/data

### Set read groups
READGROUPS="@RG\tID:$SAMPLE\tPL:ILLUMINA\tLB:$SAMPLE\tSM:$SAMPLE"

### Copy resource files to temporary sample directories.  This helps prevent multiple processes from accessing the same files and causing problems.
### Probably just need to provide these files on each node, especially if pipelines will be used by general public.
mkdir -p $PROJ/temp/$SAMPLE
cp $MILLS $PROJ/temp/$SAMPLE
cp $MILLS.idx $PROJ/temp/$SAMPLE
cp $DBSNP $PROJ/temp/$SAMPLE
cp $DBSNP.idx $PROJ/temp/$SAMPLE
cp $TG $PROJ/temp/$SAMPLE
cp $TG.idx $PROJ/temp/$SAMPLE
cp $COSMIC $PROJ/temp/$SAMPLE
cp $COSMIC.idx $PROJ/temp/$SAMPLE

### Create directory to house sample output files.
mkdir -p $PROJ/data/$SAMPLE

### Align with BWA
if [ ! -e $PROJ/data/$SAMPLE/$SAMPLE\_sortsam.bam.md5 ]
then
    bwa mem -t $DTHREADS -k 19 -w 100 -d 100 -r 1.5 -c 10000 -A 1 -B 4 -O 6 -E 1 -L 5 -U 17 -T 30 -M -R $READGROUPS $HG_OLD $RAW/$SAMPLE\_R1.fastq $RAW/$SAMPLE\_R2.fastq > $PROJ/data/$SAMPLE/$SAMPLE.bam

### Sort alignments
    java -jar -Xms36g -Xmx54g $PICARD/SortSam.jar SO=coordinate I=$PROJ/data/$SAMPLE/$SAMPLE.bam O=$PROJ/data/$SAMPLE/$SAMPLE\_sortsam.bam MAX_RECORDS_IN_RAM=$MAX_READS CREATE_MD5_FILE=true TMP_DIR=$PROJ/temp/$SAMPLE
fi

### Mark Duplicates
if [ ! -e $PROJ/data/$SAMPLE/$SAMPLE\_markdups.bam.md5 ]
then
    java -jar -Xms36g -Xmx54g $PICARD/MarkDuplicates.jar MAX_FILE_HANDLES=200 MAX_RECORDS_IN_RAM=$MAX_READS I=$PROJ/data/$SAMPLE/$SAMPLE\_sortsam.bam O=$PROJ/data/$SAMPLE/$SAMPLE\_markdups.bam M=$PROJ/data/$SAMPLE/$SAMPLE.metrics CREATE_INDEX=true CREATE_MD5_FILE=true TMP_DIR=$PROJ/temp/$SAMPLE
fi

if [ ! -e $PROJ/data/$SAMPLE/$SAMPLE\_ir.bam.md5 ]
then

### Realign Target Creator
    java -jar $GATK -R $HG -T RealignerTargetCreator -I $PROJ/data/$SAMPLE/$SAMPLE\_markdups.bam -nt $DTHREADS -o $PROJ/data/$SAMPLE/$SAMPLE.interval_list -known:indels,vcf $PROJ/temp/$SAMPLE/mills.vcf --disable_auto_index_creation_and_locking_when_reading_rods

### Indel Realigner
    java -jar $GATK -R $HG -T IndelRealigner -I $PROJ/data/$SAMPLE/$SAMPLE\_markdups.bam -targetIntervals $PROJ/data/$SAMPLE/$SAMPLE.interval_list -maxInMemory $MAX_READS -known:indels,vcf $PROJ/temp/$SAMPLE/mills.vcf -o $PROJ/data/$SAMPLE/$SAMPLE\_ir.bam --generate_md5 --disable_auto_index_creation_and_locking_when_reading_rods
fi

if [ ! -s $PROJ/data/$SAMPLE/$SAMPLE\_br1.gatk_report ]
then
### Base Recalibrator Rd. 1
    java -jar $GATK -R $HG -T BaseRecalibrator -I $PROJ/data/$SAMPLE/$SAMPLE\_ir.bam --knownSites:indels,vcf $PROJ/temp/$SAMPLE/mills.vcf --knownSites:dbsnp,vcf  $PROJ/temp/$SAMPLE/grch37_p13_142.vcf --knownSites:snps,vcf $PROJ/temp/$SAMPLE/1000g_phase1_broad_variant.vcf -o $PROJ/data/$SAMPLE/$SAMPLE\_br1.gatk_report --disable_auto_index_creation_and_locking_when_reading_rods -nct $CTHREADS
fi

if [ ! -s $PROJ/data/$SAMPLE/$SAMPLE\_br2.gatk_report ]
then
### Base Recalibrator Rd. 2
    java -jar $GATK -R $HG -T BaseRecalibrator -I $PROJ/data/$SAMPLE/$SAMPLE\_ir.bam -knownSites:indels,vcf $PROJ/temp/$SAMPLE/mills.vcf -knownSites:dbsnp,vcf $PROJ/temp/$SAMPLE/grch37_p13_142.vcf --knownSites:snps,vcf $PROJ/temp/$SAMPLE/1000g_phase1_broad_variant.vcf -BQSR $PROJ/data/$SAMPLE/$SAMPLE\_br1.gatk_report -o $PROJ/data/$SAMPLE/$SAMPLE\_br2.gatk_report --disable_auto_index_creation_and_locking_when_reading_rods -nct $CTHREADS
fi

if [ ! -e $PROJ/data/$SAMPLE/$SAMPLE\_pr.bam.md5 ]
then
### Print Reads
    java -jar $GATK -R $HG -T PrintReads -I $PROJ/data/$SAMPLE/$SAMPLE\_ir.bam -BQSR $PROJ/data/$SAMPLE/$SAMPLE\_br2.gatk_report -o $PROJ/data/$SAMPLE/$SAMPLE\_pr.bam --generate_md5 -nct $CTHREADS
fi

### Calling variants.
### Add PON (Use normals from our tumor/normal pairs.

java -jar $GATK -R $HG -T MuTect2 -I $PROJ/data/$SAMPLE/$SAMPLE\_pr.bam -o $PROJ/data/$SAMPLE/$SAMPLE.vcf -cosmic $PROJ/temp/$SAMPLE/cosmic_75_coding_hg19_sorted.vcf -dbsnp $PROJ/temp/$SAMPLE/grch37_p13_142.vcf --disable_auto_index_creation_and_locking_when_reading_rods

###
if [ "$?" = "0" ]
then
    echo "Removing temporary directory."
    rm -rf $PROJ/temp/$SAMPLE
fi
