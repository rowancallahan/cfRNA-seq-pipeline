#!/bin/bash

### Shell scripted version of exome analysis pipeline for matched tumor/normals.
### John Letaw 03/15/15

### Specify number of worker threads.
DTHREADS=20  ### Specify this.
CTHREADS=4

### Command line arguments
# MATE1=$1
# MATE2=$2
# SAM_BASE=$(basename $MATE1 .fastq)
# OUT_SAM=$SAM_BASE.sam
# SAMPLE=$3
REL=$4

### Chromosomes (to split by)
CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"

### Program directories
BIN_DIR="/opt/installed"
PICARD="$BIN_DIR/picard/picard-tools-1.110"
GATK="$BIN_DIR/GATK/GenomeAnalysisTK-3.1.jar"
MUTECT_BIN="/mnt/lustre1/CompBio/bin/mutect-1.1.7.jar"

### Project directories
PROJ="/mnt/lustre1/users/letaw/projs/test_exomes"  ### Specify this.
RESOURCE="/mnt/lustre1/CompBio/genomic_resources/"
RAW="$PROJ/raw"
LOGS="$PROJ/logs"
BWA="$PROJ/bwa"
SORTED="$PROJ/sorted"
DUPS="$PROJ/mark_dups"
INDEL="$PROJ/indel_realigned"
BQSR="$PROJ/bqsr"
MUTECT="$PROJ/mutect"
HC="$PROJ/hc"
# GG="$PROJ/genotype_gvcfs"
# VR="$PROJ/var_recal"
# PHASED="$PROJ/phased"
# DOC="$PROJ/doc"

### Reference files
HG="$RESOURCE/genomes/hg19/Homo_sapiens_assembly19.fasta"
#MILLS="$RESOURCE/variation/mills/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
#DBSNP="$RESOURCE/variation/dbsnp/138/dbsnp_138.hg19.vcf"
DBSNP="$RESOURCE/variation/dbsnp/137/dbSNP137.vcf"
THG="$RESOURCE/variation/1000g/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
HAPMAP="$RESOURCE/variation/hapmap/hapmap_3.3.hg19.sites.vcf"
OMNI="$RESOURCE/variation/omni/1000G_omni2.5.hg19.sites.vcf"
INTERVALS="/home/users/letaw/clinical/RichardsLab/resources/intervals/WES.interval_list"
#COSMIC="$RESOURCE/variation/cosmic/cosmic_hg19.vcf"
COSMIC="$RESOURCE/variation/cosmic/b37_cosmic_v54_120711.vcf"

### NA12878
#@NS500390:8:H2CWJAFXX:1:11101:2206:985 1:N:0:ACAGTG
### NA12891
#@NS500390:8:H2CWJAFXX:1:11101:17980:985 1:N:0:GCCAAT
### NA12892
#@NS500390:8:H2CWJAFXX:1:11101:3349:985 1:N:0:CAGATC

#@NS500390:9:H2CHGAFXX:1:11101:3503:1000 2:N:0:ATCACG


### Read group information.  This will be auto-created from FASTQ ID lines.
RGID=$(head -1 $MATE1 | cut -d ':' -f 1-4)
RGSM=$SAMPLE
RGPL="illumina"  ### Hard code since we will be only dealing with Illumina data through this pipeline.
RGLB=$(head -1 $MATE1 | cut -d ':' -f 3-4)
RGPU=$(head -1 $MATE1 | cut -d ':' -f 10)
RGDS=$REL
RG="@RG\tID:$RGID\tSM:$RGSM\tPL:$RGPL\tLB:$RGLB\tPU:$RGPU\tDS:$RGDS"

# ### Create directory structure.
# mkdir -p $RAW $BWA $SORTED $DUPS $INDEL $BQSR $HC $MUTECT
# #$GG $VR $PHASED $DOC

# ### Index your genome if you haven't done so already with bwa index.
# if [ ! -f $HG.sa ]
# then
#     echo "Indexing $HG."
#     bwa index $HG
# else
#     echo "$HG is already indexed."
# fi

# ### This is currently having issues when run from the command line.  Index files can't be found after job submission.
# ### This can be run by hand without issue.
# ### Align to reference.
# echo "Aligning to $HG."
# $BIN_DIR/bwa mem -M -R $RG -t $DTHREADS $HG $MATE1 $MATE2 > $BWA/$OUT_SAM

# ### Sort your SAM file.  Do not produce an index at this stage.
# echo "Sorting the SAM and outputting a BAM."
# java -jar -Xmx64g $PICARD/SortSam.jar MAX_RECORDS_IN_RAM=10000000 SO=coordinate I=$BWA/$OUT_SAM O=$SORTED/$SAM_BASE\_sorted.bam

# ### Mark duplicates, create and index.
# echo "Marking duplicates."
# java -jar -Xmx64g $PICARD/MarkDuplicates.jar MAX_RECORDS_IN_RAM=10000000 CREATE_INDEX=TRUE I=$SORTED/$SAM_BASE\_sorted.bam O=$DUPS/$SAM_BASE\_sorted_dups.bam M=$DUPS/$SAM_BASE.metrics

# ### Create the interval file, if it doesn't already exist.
# if [ ! -f $PROJ/target_intervals.list ] && [ "$4" == "tumor" ] 
# then
#     echo "Creating intervals for IndelRealigner."
#     java -jar $GATK -R $HG -T RealignerTargetCreator -I $DUPS/$SAM_BASE\_sorted_dups.bam -known $MILLS -o $PROJ/target_intervals.list -nt $DTHREADS
# else
#     while [ ! -f $PROJ/target_intervals.list ]
#     do
# 	sleep 5
#     done
#     echo "Interval file has already been created, proceeding to IndelRealigner."
# fi

# ### Run a realignment around targeted indels.
# echo "Realigning around targeted indel regions."
# java -jar $GATK -R $HG -T IndelRealigner -I $DUPS/$SAM_BASE\_sorted_dups.bam -targetIntervals $PROJ/target_intervals.list -known $MILLS -o $INDEL/$SAM_BASE\_sorted_dups_realigned.bam

# ### Base recalibration.
# echo "Base recalibration round 1."
# java -jar $GATK -R $HG -T BaseRecalibrator -I $INDEL/$SAM_BASE\_sorted_dups_realigned.bam -knownSites $DBSNP -knownSites $MILLS -knownSites $THG -o $BQSR/$SAM_BASE\_recal_data.table -nct $CTHREADS

# echo "Base recalibration round 2."
# java -jar $GATK -R $HG -T BaseRecalibrator -I $INDEL/$SAM_BASE\_sorted_dups_realigned.bam -knownSites $DBSNP -knownSites $MILLS -knownSites $THG -BQSR $BQSR/$SAM_BASE\_recal_data.table -o $BQSR/$SAM_BASE\_post_recal_data.table -nct $CTHREADS

# # echo "Base recalibration plots."
# # java -jar $GATK -R $HG -T AnalyzeCovariates -before $BQSR/$SAM_BASE\_recal_data.table -after $BQSR/$SAM_BASE\_post_recal_data.table -plots $BQSR/$SAM_BASE\_recal_plot.pdf

# echo "Applying recalibration."
# java -jar $GATK -R $HG -T PrintReads -I $INDEL/$SAM_BASE\_sorted_dups_realigned.bam -BQSR $BQSR/$SAM_BASE\_post_recal_data.table -o $BQSR/$SAM_BASE\_sorted_dups_realigned_bqsr.bam -nct $CTHREADS

# ### Calling variants.
# ### Add -D $DBSNP option here?
# echo "Calling variants with HaplotypeCaller."
# java -jar $GATK -R $HG -T HaplotypeCaller -I $1 -o $HC/$REL.vcf -nct $CTHREADS

### Merge and genotype gVCFs.

if [ "$REL" == "tumor" ]
then
    echo "Running MuTect."

    for log in $(ls -1 $LOGS/*.{1}.log)
    do
	echo "Waiting for $log."
	until grep -q "return value 100" $log
	do
	    sleep 5
	done
    done

    NF=$(basename $2 .bam)
    java -jar $MUTECT_BIN -R $HG -T MuTect --dbsnp_normal_lod $1 --tumor_lod 10 --enable_extended_output -I:tumor $2 -I:normal $3 -dbsnp $DBSNP -cosmic $COSMIC -vcf $MUTECT/mutect_$NF\_$1.vcf -o $MUTECT/call_stats_$NF\_$1
else
    echo "DONE"
    exit 100
fi









### Not yet applicable


# ### Build SNP recalibration model.
# echo "Building the SNP recalibration model for VQSR."
# java -jar $GATK -R $HG -T VariantRecalibrator -input $GG/trio.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $THG -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an DP -an QD -an FS -an MQ -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $VR/recalibrate_SNP.recal -tranchesFile $VR/recalibrate_SNP.tranches -rscriptFile $VR/recalibrate_SNP_plots.R -nt $DTHREADS

# ### Apply recalibration.
# echo "Applying VQSR."
# java -jar $GATK -R $HG -T ApplyRecalibration -input $GG/trio.vcf -mode SNP --ts_filter_level 99.0 -recalFile $VR/recalibrate_SNP.recal -tranchesFile $VR/recalibrate_SNP.tranches -o $VR/trio_recal_snps_raw_indels.vcf -nt $DTHREADS

# ### Build INDEL recalibration model
# echo "Building the indel recalibration model for VQSR."
# java -jar $GATK -T VariantRecalibrator -R $HG -input $VR/trio_recal_snps_raw_indels.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS -an QD -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 -recalFile $VR/recalibrate_INDEL.recal -tranchesFile $VR/recalibrate_INDEL.tranches -rscriptFile $VR/recalibrate_INDEL_plots.R -nt $DTHREADS 

# ### Apply recalibration.
# echo "Applying VQSR."
# java -jar $GATK -T ApplyRecalibration -R $HG -input $VR/trio_recal_snps_raw_indels.vcf -mode INDEL --ts_filter_level 99.0 -recalFile $VR/recalibrate_INDEL.recal -tranchesFile $VR/recalibrate_INDEL.tranches -o $VR/trio_recalibrated_variants.vcf 

# ### Apply phasing based on pedigree info.
# echo "Applying phasing."
# java -jar $GATK -T PhaseByTransmission -R $HG -V $VR/trio_recalibrated_variants.vcf -ped $PROJ/raw/pedigree -o $PHASED/final_trio.vcf -mvf $VR/trio_mend_viol

# ### Depth Of Coverage.
# echo "Running DepthOfCoverage."
# java -jar $GATK -R $HG -T DepthOfCoverage -I $BQSR/$SAM_BASE\_sorted_dups_realigned_bqsr.bam -pt sample --outputFormat rtable -L $INTERVALS -isr UNION -baqGOP 40 -DBQ 0 -S STRICT -im ALL --maxBaseQuality 127 -mbq 20 --maxMappingQuality 2147483647 -mmq 20 --nBins 499 -omitIntervals --start 1 --stop 500 -o $DOC/$SAM_BASE\_doc
