#!/bin/bash

# Look through raw file directory to find names matching samples.
for FILE in "${FILES[@]}"; do
    if [[ $FILE == *"$SAMPLE"* ]]; then
	if [[ $FILE =~ $READ1_REGEX ]]; then
	    READ1="$FILE"
	elif [[ $FILE =~ $READ2_REGEX ]]; then
	    READ2="$FILE"
	else
	    echo "No R1 or R2 designation, assuming single-end."
	    READ1="$FILE"
	fi
    fi
done

# Get read groups.
if [ ${READ1: -3} == ".gz" ]; then
    FASTQ_ID=$(gzip -cd $READ1 | head -1)
else
    FASTQ_ID=$(head -1 $READ1)
fi

# Overwrite defaults as defined by read_groups function.
read_groups $FASTQ_ID $SAMPLE
#RGSM=$SAMPLE
#RGDS=$WORKFLOW\_$VERSION

# Run STAR.
if [ ! -e $PROJ/output/$SAMPLE/$SAMPLE\Aligned.sortedByCoord.out.bam ]; then
    if [ ${READ1: -3} == '.gz' ] && [ -z ${READ2} ]; then
	STAR --runThreadN $DTHREADS --genomeDir $INDEX_DIR --readFilesIn "$READ1" "$READ2" --readFilesCommand zcat --outFileNamePrefix $PROJ/output/$SAMPLE/$SAMPLE --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattrRGline $STAR_READGROUPS
    elif [ ${READ1: -3} != '.gz' ] && [ -z ${READ2} ]; then
	STAR --runThreadN $DTHREADS --genomeDir $INDEX_DIR --readFilesIn "$READ1" "$READ2" --outFileNamePrefix $PROJ/output/$SAMPLE/$SAMPLE --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattrRGline $STAR_READGROUPS
    elif [ ${READ1: -3} == '.gz' ]; then
	STAR --runThreadN $DTHREADS --genomeDir $INDEX_DIR --readFilesIn "$READ1" --readFilesCommand zcat --outFileNamePrefix $PROJ/output/$SAMPLE/$SAMPLE --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattrRGline $STAR_READGROUPS
    else
	STAR --runThreadN $DTHREADS --genomeDir $INDEX_DIR --readFilesIn "$READ1" --outFileNamePrefix $PROJ/output/$SAMPLE/$SAMPLE --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattrRGline $STAR_READGROUPS
    fi
fi

# Index the BAM output.
if [ ! -e $PROJ/output/$SAMPLE/$SAMPLE\Aligned.sortedByCoord.out.bam.bai ]
then
    $JAVA8 $JAVAVM -jar $PICARD BuildBamIndex I=$PROJ/output/$SAMPLE/$SAMPLE\Aligned.sortedByCoord.out.bam CREATE_MD5_FILE=true
fi

# SplitNCigarReads step from GATK best practices.
if [ ! -e $PROJ/output/$SAMPLE/splitn.bam.md5 ]
then
    $JAVA8 $JAVAVM -jar $GATK -T SplitNCigarReads -R $GENOME_FASTA -I $PROJ/output/$SAMPLE/$SAMPLE\Aligned.sortedByCoord.out.bam -o $PROJ/output/$SAMPLE/splitn.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --generate_md5
fi

# # Run IndelRealigner step.
# if [ ! -e $PROJ/output/$SAMPLE/ir.bam.md5 ]
# then
#     $JAVA8 $JAVAVM -jar $GATK -T RealignerTargetCreator -R $GENOME_FASTA -I $PROJ/output/$SAMPLE/splitn.bam -o $PROJ/output/$SAMPLE/rtc.interval_list
#     $JAVA8 $JAVAVM -jar $GATK -T IndelRealigner -R $GENOME_FASTA -I $PROJ/output/$SAMPLE/splitn.bam -targetIntervals $PROJ/output/$SAMPLE/rtc.interval_list -o $PROJ/output/$SAMPLE/ir.bam --generate_md5
# fi

# BaseRecalibrator and PrintReads.
if [ ! -e $PROJ/output/$SAMPLE/pr.bam.md5 ]
then
    $JAVA8 $JAVAVM -jar $GATK -T BaseRecalibrator -R $GENOME_FASTA -I $PROJ/output/$SAMPLE/splitn.bam -o $PROJ/output/$SAMPLE/recal.txt -knownSites $DBSNP -nct $CTHREADS

    $JAVA8 $JAVAVM -jar $GATK -T BaseRecalibrator -R $GENOME_FASTA -I $PROJ/output/$SAMPLE/splitn.bam -o $PROJ/output/$SAMPLE/post_recal.txt -knownSites $DBSNP -BQSR $PROJ/output/$SAMPLE/recal.txt -nct $CTHREADS
    
    $JAVA8 $JAVAVM -jar $GATK -T PrintReads -R $GENOME_FASTA -BQSR $PROJ/output/$SAMPLE/post_recal.txt -I $PROJ/output/$SAMPLE/splitn.bam -o $PROJ/output/$SAMPLE/pr.bam -nct $CTHREADS --generate_md5
fi

# Variant calling.
$JAVA8 $JAVAVM -jar $GATK -T HaplotypeCaller -R $GENOME_FASTA -I $PROJ/output/$SAMPLE/pr.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $PROJ/output/$SAMPLE/hc.vcf -D $DBSNP -nct $CTHREADS
#-emitRefConfidence GVCF
