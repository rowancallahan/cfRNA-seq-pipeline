#!/bin/bash

FILE_LIST=$(for line in $(ls -1 $PROJ/hc/*.gvcf); do printf -- " -V %s " "$line"; done)

if [ ! -e $PROJ/genotype_gvcfs/$VCF_NAME.vcf.idx ]
then
    $JAVA8 -Xms8g -Xmx16g -jar $GATK -T GenotypeGVCFs -R $HG -D $DBSNP -stand_call_conf 20.0 -o $PROJ/genotype_gvcfs/$VCF_NAME.vcf -nt $DTHREADS --disable_auto_index_creation_and_locking_when_reading_rods $FILE_LIST
fi

### This snippet will run VQSR on an appropriate dataset.  Please see the following links for additional information.
### All output files will be placed in the project directory and the output file names will be based off of the
### basename (file without suffix) of the input VCF.
### Note: All resource files should be indexed.  Also, if running on the cluster, it is highly recommended to utilize the
### --disable_auto_index_creation_and_locking_when_reading_rods option as included below.

### http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr
### https://www.broadinstitute.org/gatk/guide/article?id=2805

# Removing MQRankSum annotation for testing, though it could be included for pull dataset.
# Bad input: Found annotations with zero variance. They must be excluded before proceeding.
# -an MQRankSum
# Disabled during testing for following error:
# org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException: Unable to retrieve result
# -nt $DTHREADS 

$JAVA8 -Xmx16g -jar $GATK -R $HG -T VariantRecalibrator --input $PROJ/genotype_gvcfs/$VCF_NAME.vcf -recalFile $PROJ/vqsr/$VCF_NAME\_snps.recal -tranchesFile $PROJ/vqsr/$VCF_NAME\_snps.tranches -rscriptFile $VCF_NAME\_snps_plots.R -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $TG -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an QD -an MQ -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --disable_auto_index_creation_and_locking_when_reading_rods -an MQRankSum -nt $DTHREADS

$JAVA8 -Xmx16g -jar $GATK -T ApplyRecalibration -R $HG -input $PROJ/genotype_gvcfs/$VCF_NAME.vcf -tranchesFile $PROJ/vqsr/$VCF_NAME\_snps.tranches -recalFile $PROJ/vqsr/$VCF_NAME\_snps.recal -o $PROJ/vqsr/$VCF_NAME\_snps_vqsr.vcf --ts_filter_level 99.5 -mode SNP --disable_auto_index_creation_and_locking_when_reading_rods

$JAVA8 -Xmx16g -jar $GATK -R $HG -T VariantRecalibrator --input $PROJ/vqsr/$VCF_NAME\_snps_vqsr.vcf -recalFile $PROJ/vqsr/$VCF_NAME\_indels.recal -tranchesFile $PROJ/vqsr/$VCF_NAME\_indels.tranches -rscriptFile $VCF_NAME\_indels_plots.R -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS -an QD -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 --disable_auto_index_creation_and_locking_when_reading_rods -an MQRankSum -nt $DTHREADS

$JAVA8 -Xmx16g -jar $GATK -T ApplyRecalibration -R $HG -input $PROJ/vqsr/$VCF_NAME\_snps_vqsr.vcf -tranchesFile $PROJ/vqsr/$VCF_NAME\_indels.tranches -recalFile $PROJ/vqsr/$VCF_NAME\_indels.recal -o $PROJ/vqsr/$VCF_NAME\_recalibrated.vcf --ts_filter_level 99.0 -mode INDEL --disable_auto_index_creation_and_locking_when_reading_rods
