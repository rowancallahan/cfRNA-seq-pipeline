#!/bin/bash

# Call read_groups to create READGROUPS string.
if [ -e $RAW/$SAMPLE\_R1.fastq.gz ] || [ -e $RAW/$SAMPLE\_R1.fq.gz ]; then
    FASTQ_HEADER=$(gzip -cd $RAW/$SAMPLE\_R1.f*.gz | head -1)
else
    FASTQ_HEADER=$(head -1 $RAW/$SAMPLE\_R1.f*)
fi
read_groups $FASTQ_HEADER $SAMPLE

### Align with BWA
if [ ! -e $PROJ/sortsam/$SAMPLE\_sortsam.bam.md5 ]
then
    bwa mem -t $DTHREADS -k 19 -w 100 -d 100 -r 1.5 -c 10000 -A 1 -B 4 -O 6 -E 1 -L 5 -U 17 -T 30 -M -R "$READGROUPS" $HG $RAW/$SAMPLE\_R1.fastq.gz $RAW/$SAMPLE\_R2.fastq.gz > $PROJ/bwa/$SAMPLE.bam

### Sort alignments
    $JAVA8 -jar -Xms36g -Xmx54g $PICARD/SortSam.jar SO=coordinate I=$PROJ/bwa/$SAMPLE.bam O=$PROJ/sortsam/$SAMPLE\_sortsam.bam MAX_RECORDS_IN_RAM=$MAX_READS CREATE_MD5_FILE=true
fi

### Mark Duplicates
if [ ! -e $PROJ/markdups/$SAMPLE\_markdups.bam.md5 ]
then
    $JAVA8 -jar -Xms36g -Xmx54g $PICARD/MarkDuplicates.jar MAX_FILE_HANDLES=200 MAX_RECORDS_IN_RAM=$MAX_READS I=$PROJ/sortsam/$SAMPLE\_sortsam.bam O=$PROJ/markdups/$SAMPLE\_markdups.bam M=$PROJ/markdups/$SAMPLE.metrics CREATE_INDEX=true CREATE_MD5_FILE=true
fi

# GATK has decided they don't like the indel realigner any more, and have suggested not using it.
# if [ ! -e $PROJ/ir/$SAMPLE\_ir.bam.md5 ]
# then
# ### Realign Target Creator
#     $JAVA8 -jar $GATK -R $HG -T RealignerTargetCreator -I $PROJ/markdups/$SAMPLE\_markdups.bam -nt $DTHREADS -o $PROJ/rtc/$SAMPLE.interval_list -known:indels,vcf $MILLS --disable_auto_index_creation_and_locking_when_reading_rods

# ### Indel Realigner
#     $JAVA8 -jar $GATK -R $HG -T IndelRealigner -I $PROJ/markdups/$SAMPLE\_markdups.bam -targetIntervals $PROJ/rtc/$SAMPLE.interval_list -maxInMemory $MAX_READS -known:indels,vcf $MILLS -o $PROJ/ir/$SAMPLE\_ir.bam --generate_md5 --disable_auto_index_creation_and_locking_when_reading_rods
# fi

# Base recalibration.
if [ ! -e $PROJ/pr/$SAMPLE\_pr.bam.md5 ]
then
    $JAVA8 -jar $GATK -R $HG -T BaseRecalibrator -L $PROBES -I $PROJ/markdups/$SAMPLE\_markdups.bam --knownSites:indels,vcf $MILLS --knownSites:dbsnp,vcf $DBSNP --knownSites:snps,vcf $TG -o $PROJ/br/$SAMPLE\_br1.gatk_report --disable_auto_index_creation_and_locking_when_reading_rods
    $JAVA8 -jar $GATK -R $HG -T BaseRecalibrator -L $PROBES -I $PROJ/markdups/$SAMPLE\_markdups.bam -knownSites:indels,vcf $MILLS -knownSites:dbsnp,vcf $DBSNP --knownSites:snps,vcf $TG -BQSR $PROJ/br/$SAMPLE\_br1.gatk_report -o $PROJ/br/$SAMPLE\_br2.gatk_report --disable_auto_index_creation_and_locking_when_reading_rods
    $JAVA8 -jar $GATK -R $HG -T PrintReads -I $PROJ/markdups/$SAMPLE\_markdups.bam -BQSR $PROJ/br/$SAMPLE\_br2.gatk_report -o $PROJ/pr/$SAMPLE\_pr.bam --generate_md5
fi

# Calling variants.
$JAVA8 -jar $GATK -R $HG -T HaplotypeCaller -L $PROBES --emitRefConfidence GVCF -I $PROJ/pr/$SAMPLE\_pr.bam -o $PROJ/hc/$SAMPLE.g.vcf

# May want to rethink this, depending on how the subsequent gVCF genotyping stage is carried out.
# if [ "$?" = "0" ]
# then
#     echo "Removing temporary directory."
#     rm -rf $PROJ/temp/$SAMPLE
# fi
