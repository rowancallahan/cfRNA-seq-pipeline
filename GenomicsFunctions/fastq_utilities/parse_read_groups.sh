#!/bin/bash

# USAGE: bash parse_read_groups.sh <YOUR_FASTQ_ID>
# Or, include it in your workflow script.
VERSION='0.1.2'

# READ_GROUP="@D00735:131:C9DWPANXX:7:1109:1198:2063 1:N:0:GGCTAC"
# READ_GROUP=$1

read_groups () {
    
    # This will supply read groups that can be determined from the
    # information in a FASTQ record ID line.  Also will add blank values
    # for the other read groups.  These can be overwritten if needed.

# RGID=String
# ID=String                     Read Group ID  Default value: 1. This option can be set to 'null' to clear the default 
#                               value. 

# RGLB=String
# LB=String                     Read Group Library  Required. 

# RGPL=String
# PL=String                     Read Group platform (e.g. illumina, solid)  Required. 

# RGPU=String
# PU=String                     Read Group platform unit (eg. run barcode)  Required. 

# RGSM=String
# SM=String                     Read Group sample name  Required. 

# RGCN=String
# CN=String                     Read Group sequencing center name  Default value: null. 

# RGDS=String
# DS=String                     Read Group description  Default value: null. 

# RGDT=Iso8601Date
# DT=Iso8601Date                Read Group run date  Default value: null. 

# RGPI=Integer
# PI=Integer                    Read Group predicted insert size  Default value: null. 

    FASTQ=$1
    BARCODE=$2
    SAMPLE=$3
    arrFASTQ=(${FASTQ//:/ })
    arrBARCODE=(${BARCODE//:/ })

    # Set RGID and RGPU based on FASTQ ID.
    RGID=$(printf "%s.%s" "${arrFASTQ[0]}" "${arrFASTQ[3]}")
    RGPU=$(printf "%s.%s.%s" "${arrFASTQ[0]}" "${arrFASTQ[2]}" "${arrFASTQ[3]}")
    # Set the date according to when the data is processed.
    # If desired, set the date to when the sequencing was actually performed.
    RGDT=$(date --iso-8601)
    # Set these as default since our data commonly falls in to these categories.
    RGPL="ILLUMINA"
    RGCN="OHSU_IGL"
    # These two can be set as the same thing, unless you have prepared the same
    # sample with multiple library methods, then you would vary RGLB.  We will by default
    # set RGLB and RGSM as the barcode, though RGSM should certainly be altered to be more
    # human readable.
    RGSM=$(printf "%s" "$SAMPLE")
#    RGSM=$(printf "%s" "${arrBARCODE[3]}")
    RGLB=$(printf "%s" "${arrBARCODE[3]}")
    # RGDS should be set according to the assay that is being run.
    RGDS="SMMART"
    # Potentially can set to median, using Picard or Samtools.
    # Must be numeric!
    RGPI=0

    READGROUPS="@RG\tID:$RGID\tPL:$RGPL\tLB:$RGLB\tSM:$RGSM\tPU:$RGPU\tDT:$RGDT\tCN:$RGCN\tDS:$RGDS\tPI:$RGPI"
    # Output a string for STAR also.
    STAR_READGROUPS="ID:$RGID PL:$RGPL LB:$RGLB SM:$RGSM PU:$RGPU DT:$RGDT CN:$RGCN "DS:$RGDS" PI:$RGPI"

}
