#!/bin/bash

### Specify number of worker threads.
DTHREADS=24
CTHREADS=4
MAX_READS=14000000

### Program directories
BIN_DIR="/opt/installed"
PICARD="$BIN_DIR/picard/picard-tools-1.110"
GATK="$BIN_DIR/GATK/GenomeAnalysisTK-3.1.jar"

### Project directories
PROJ="/home/users/letaw/lustre1/projs/trio_related"
RESOURCE="/mnt/lustre1/CompBio/genomic_resources"
HG="$RESOURCE/genomes/hg19/Homo_sapiens_assembly19.fasta"
RAW="$PROJ/raw"
LOGS="$PROJ/logs"

### Resource Files
DBSNP="$RESOURCE/variation/dbsnp/142/grch37_p13_142.vcf"
MILLS="$RESOURCE/variation/mills/mills.vcf"
TG="$RESOURCE/variation/1000g/1000g_phase1_broad_variant.vcf"

SAMPLE=$1

if [ $SAMPLE == "99999" ]; then
    rm -rf $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/ir $PROJ/br $PROJ/pr $PROJ/hc
    rm $LOGS/*
    exit
fi

### Create directory structure.
mkdir -p $PROJ/bwa $PROJ/sortsam $PROJ/markdups $PROJ/rtc $PROJ/ir $PROJ/br $PROJ/pr $PROJ/hc

### Set read groups
READGROUPS="@RG\tID:$SAMPLE\tPL:ILLUMINA\tLB:$SAMPLE\tSM:$SAMPLE"

### Align with BWA
bwa mem -t $DTHREADS -k 19 -w 100 -d 100 -r 1.5 -c 10000 -A 1 -B 4 -O 6 -E 1 -L 5 -U 17 -T 30 -M -R $READGROUPS $HG $RAW/$SAMPLE\_R1.fastq.gz $RAW/$SAMPLE\_R2.fastq.gz > $PROJ/bwa/$SAMPLE.bam

### Sort alignments
java -jar -Xms36g -Xmx54g $PICARD/SortSam.jar SO=coordinate I=$PROJ/bwa/$SAMPLE.bam O=$PROJ/sortsam/$SAMPLE\_sortsam.bam MAX_RECORDS_IN_RAM=$MAX_READS

### Mark Duplicates
java -jar -Xms36g -Xmx54g $PICARD/MarkDuplicates.jar MAX_FILE_HANDLES=200 MAX_RECORDS_IN_RAM=$MAX_READS I=$PROJ/sortsam/$SAMPLE\_sortsam.bam O=$PROJ/markdups/$SAMPLE\_markdups.bam M=$PROJ/markdups/$SAMPLE.metrics CREATE_INDEX=true

### Realign Target Creator
java -jar $GATK -R $HG -T RealignerTargetCreator -I $PROJ/markdups/$SAMPLE\_markdups.bam -nt $DTHREADS -o $PROJ/rtc/$SAMPLE.interval_list -known $MILLS

### Indel Realigner
java -jar $GATK -R $HG -T IndelRealigner -I $PROJ/markdups/$SAMPLE\_markdups.bam -targetIntervals $PROJ/rtc/$SAMPLE.interval_list -maxInMemory $MAX_READS -known $MILLS -o $PROJ/ir/$SAMPLE\_ir.bam

### Base Recalibrator Rd. 1
java -jar $GATK -R $HG -T BaseRecalibrator -I $PROJ/ir/$SAMPLE\_ir.bam -knownSites $MILLS -knownSites $DBSNP -knownSites $TG -o $PROJ/br/$SAMPLE\_br1.gatk_report -nct $CTHREADS
 
### Base Recalibrator Rd. 2
java -jar $GATK -R $HG -T BaseRecalibrator -I $PROJ/ir/$SAMPLE\_ir.bam -knownSites $MILLS -knownSites $DBSNP -knownSites $TG -BQSR $PROJ/br/$SAMPLE\_br1.gatk_report -o $PROJ/br/$SAMPLE\_br2.gatk_report -nct $CTHREADS

### Print Reads
java -jar $GATK -R $HG -T PrintReads -I $PROJ/ir/$SAMPLE\_ir.bam -BQSR $PROJ/br/$SAMPLE\_br2.gatk_report -o $PROJ/pr/$SAMPLE\_pr.bam -nct $CTHREADS

### Calling variants.
java -jar $GATK -R $HG -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -I $PROJ/pr/$SAMPLE\_pr.bam -o $PROJ/hc/$SAMPLE.gvcf -nct $CTHREADS --disable_auto_index_creation_and_locking_when_reading_rods


# ### Merge and genotype gVCFs.

# if [ "$SAMPLE" == "proband" ]
# then
#     echo "Merging files."

#     for log in $(ls -1 $LOGS/*.{1,2}.log)
#     do
# 	echo "Waiting for $log."
# 	until grep -q "return value 100" $log
# 	do
# 	    sleep 5
# 	done
#     done

#     java -jar $GATK -R $HG -T GenotypeGVCFs `for line in $(ls -1 $PROJ/*.gvcf); do echo -n "-V $line "; done` -D $DBSNP -o $PROJ/trio.vcf -nt $DTHREADS
# else
#     echo "DONE"
#     exit 100
# fi

# ### Hard filtering
# echo "Applying filters."
# java -jar $GATK -T VariantFiltration -R $HG -V $PROJ/trio.vcf -o $PROJ/trio_filtered.vcf -filterName "QDFilter" -filter "QD < 5" -filterName "QUALFilter" -filter "QUAL < 100" -filterName "DPFilter" -filter "DP < 10" -cluster 3 -window 15
 
# ### Apply phasing based on pedigree info.
# echo "Applying phasing."
# java -jar $GATK -T PhaseByTransmission -R $HG -V $PROJ/trio_filtered.vcf -ped $PROJ/pedigree -o $PROJ/final_trio.vcf -mvf $PROJ/trio_mend_viol
