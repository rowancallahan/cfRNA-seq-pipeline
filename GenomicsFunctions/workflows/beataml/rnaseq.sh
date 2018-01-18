#!/bin/bash

### RNA-seq DE pipeline script. (GSNAP based)
# Run this on one of our ACC servers (mother, hoolock, fergie2) or on tiger.
# Separate scripts are necessary for exacloud, as jobs need to be split.
###
# Set directories at top of script.  If you are unsure about using this by pushing the enter button, perform
# each step individually.
# Read groups are contained within a CSV placed within the script directory.  These can be produced in Excel
# if you want and exported as a CSV.  Please retain ordering of read groups as stated on the MarkDuplicates
# help page.
# Program and genome directories are appropriate for current installations on hoolock.
# Raw FASTQ files should be concatenated per sample and placed within the raw subdirectory of your project folder.
###
# John Letaw 12/31/14
# letaw@ohsu.edu
###

### Where the genome FASTA's are.
DB_DIR="/u0/letaw/genomes/mmu_grcm38"
### The reference genome.
REF_FA="$DB_DIR/Mus_musculus.GRCm38.dna.primary_assembly.fa"
### Location of PICARD
PICARD="/usr/local/bin/picard-tools-1.119"
### Location of Trimmomatic
TRIMM="/usr/local/bin/trimmomatic-0.30.jar"
### Location of guess_encoding script.
PHRED="/u0/letaw/scripts/guess_encoding.py"
### Prefix used for index files, and splicesites files.
GENOME="mmu_38"
### Gene set.
GTF="$DB_DIR/Mus_musculus.GRCm38.78.gtf"
### Directory containing project structure.
PROJECT_DIR="/u0/letaw/projs/testrun"
PROJ_NAME="testrun"
### Option field for Trimmomatic.  Look at docs for more info.
### http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.30.pdf
TRIMM_OPT="MAXINFO:50:.9 MINLEN:40"
### For HTSeq-Count, is this a strand-specific assay?  It could be yes, no, or reverse.  See HTSeq-count help for more.
STRANDED="yes"
### How many processing threads?
THREADS=18

function add_dir {

    if [ ! -d $PROJECT_DIR/$1 ]
    then
	mkdir $PROJECT_DIR/$1
    fi
}

### Unset variables FWD and REV in case they are set in environment.
unset REV
unset FWD

### Set nullglob to avoid globs with zero elements appearing to be one element.
shopt -s nullglob

### Clean old directory structure if 666 is specified on command line.
if [ $1 == 66 ]
then
    rm -r $PROJECT_DIR/fastqc $PROJECT_DIR/gsnap $PROJECT_DIR/rg $PROJECT_DIR/dups $PROJECT_DIR/htseq $PROJECT_DIR/trimm
    exit 1
fi

### Create project directory structure.

add_dir fastqc
add_dir gsnap
add_dir rg
add_dir dups
add_dir htseq

### Build GSNAP index, if necessary.

if [ ! -d $DB_DIR/$GENOME ]
then
    gmap_build -d $GENOME -D $DB_DIR $REF_FA
fi

### Create splicesites files, assuming a GTF exists.
if [ ! -e $DB_DIR/$GENOME/$GENOME.maps/$GENOME.iit ]
then
    cat $GTF | gtf_splicesites > $DB_DIR/$GENOME/$GENOME.maps/$GENOME.splicesites
    cat $DB_DIR/$GENOME/$GENOME.maps/$GENOME.splicesites | iit_store -o $DB_DIR/$GENOME/$GENOME.maps/$GENOME.iit
    rm $DB_DIR/$GENOME/$GENOME.maps/$GENOME.splicesites
fi
	
### Run FastQC.
###
fastqc --noextract -o $PROJECT_DIR/fastqc -t $THREADS $PROJECT_DIR/raw/*.gz

### Trimmomatic, if desired.  Command line value of "0" indicates no trimming, "1" for single-end, "2" for paired-end.
### Also, check the quality score encoding.

if [ $1 == 0 ]
then
    echo "No trimming will be performed.  Did you look at FastQC output?"
    FOR_GSNAP="$PROJECT_DIR/raw/*.gz"

elif [ $1 == 1 ]
then
    echo "Trimming will be performed on single-end reads."
    add_dir trimm

    for line in $(ls -1 $PROJECT_DIR/raw/*.gz)
    do
	QUAL_SCORE=$(gzip -cd $line | awk 'NR % 4 == 0' | python $PHRED)
	FILENAME=$(basename $line .fastq.gz)

	java -jar $TRIMM SE -threads $THREADS -$QUAL_SCORE -trimlog $PROJECT_DIR/trimm/$FILENAME.log $line $PROJECT_DIR/trimm/$FILENAME.fastq $TRIMM_OPT
	gzip -fq $PROJECT_DIR/trimm/$FILENAME.fastq
    done

    FOR_GSNAP="$PROJECT_DIR/trimm/*.gz"

elif [ $1 == 2 ]
then
    echo "Trimming will be performed on paired-end reads."
    add_dir trimm

    cd $PROJECT_DIR/raw

    FWD=( *R1*gz )
    REV=( *R2*gz )

    START=0
    END=$(expr ${#FWD[@]} - 1)

    for (( i=$START; i<=$END; i++ ))
    do
	QUAL_SCORE=$(gzip -cd ${FWD[$i]} | awk 'NR % 4 == 0' | python $PHRED)
	FILENAME1=$(basename ${FWD[$i]} .fastq.gz)
	FILENAME2=$(basename ${REV[$i]} .fastq.gz)
	echo $FILENAME1
	echo $FILENAME2
	java -jar $TRIMM PE -threads $THREADS -$QUAL_SCORE -trimlog $PROJECT_DIR/trimm/$FILENAME1.log ${FWD[$i]} ${REV[$i]} $PROJECT_DIR/trimm/$FILENAME1.fastq $PROJECT_DIR/trimm/$FILENAME1.unpaired_out $PROJECT_DIR/trimm/$FILENAME2.fastq $PROJECT_DIR/trimm/$FILENAME2.unpaired_out $TRIMM_OPT
	gzip -fq $PROJECT_DIR/trimm/$FILENAME1.fastq
	gzip -fq $PROJECT_DIR/trimm/$FILENAME2.fastq

    done
    
    FOR_GSNAP="$PROJECT_DIR/trimm"

else
    echo "Argument must be 0 (no trimming), 1 (single-end), or 2 (paired-end)."
    exit 1
fi

### Run GSNAP alignments.
### This is currently just for single-end reads, need to adapt for paired.

### If a FWD variable does not exist, run the single-end version.
if [ -z ${FWD+x} ]
then
    for line in $(ls -1 $FOR_GSNAP)
    do
	INPUT_DIR=$(basename $line .fastq.gz)
	echo $INPUT_DIR
	mkdir $PROJECT_DIR/gsnap/$INPUT_DIR
###
	gsnap -D $DB_DIR -d $GENOME -s $GENOME --use-sarray=1 -B 5 -m .05 --gunzip -t $THREADS -N 1 -A sam --failed-input $PROJECT_DIR/gsnap/$INPUT_DIR/$INPUT_DIR.failed --split-output=$PROJECT_DIR/gsnap/$INPUT_DIR/$INPUT_DIR $line
    done

### Otherwise, run the paired-end version.
### This should be more extensively tested, though it has worked for me so far.
else

    FWD=( $FOR_GSNAP/*R1*gz )
    REV=( $FOR_GSNAP/*R2*gz )

    START=0
    END=$(expr ${#FWD[@]} - 1)

    for (( i=$START; i<=$END; i++ ))
    do
	INPUT_DIR=$(basename ${FWD[$i]} .fastq.gz)
	mkdir $PROJECT_DIR/gsnap/$INPUT_DIR
###
	gsnap -D $DB_DIR -d $GENOME -s $GENOME --use-sarray=1 -B 5 -m .05 --gunzip -t $THREADS -N 1 -A sam --failed-input $PROJECT_DIR/gsnap/$INPUT_DIR/$INPUT_DIR.failed --split-output=$PROJECT_DIR/gsnap/$INPUT_DIR/$INPUT_DIR ${FWD[$i]} ${REV[$i]}
	
    done
    
fi

### This chunk will read through a CSV and add read groups.
COUNTER=1
for aligned in $(ls -1d $PROJECT_DIR/gsnap/*)
do

    CURR_RG=$(head -$COUNTER $PROJECT_DIR/read_groups.txt | tail -1)
    IFS=',' read -ra RG <<< "$CURR_RG"

    ### Add read groups and coordinate sort SAM's.

    if [ -z ${FWD+x} ]
    then
	FILE_NAME=$(basename $aligned/*unpaired_uniq .unpaired_uniq)
###
	java -jar $PICARD/AddOrReplaceReadGroups.jar INPUT=$PROJECT_DIR/gsnap/$FILE_NAME/$FILE_NAME.unpaired_uniq OUTPUT=$PROJECT_DIR/rg/$FILE_NAME\_rg.sam SORT_ORDER=coordinate RGID=${RG[0]} RGLB=${RG[1]} RGPL=${RG[2]} RGPU=${RG[3]} RGSM=${RG[4]} RGCN=${RG[5]} RGDS=${RG[6]} RGDT=${RG[7]} RGPI=${RG[8]} 
    else
	FILE_NAME=$(basename $aligned/*concordant_uniq .concordant_uniq)
###
	java -jar $PICARD/AddOrReplaceReadGroups.jar INPUT=$PROJECT_DIR/gsnap/$FILE_NAME/$FILE_NAME.concordant_uniq OUTPUT=$PROJECT_DIR/rg/$FILE_NAME\_rg.sam SORT_ORDER=coordinate RGID=${RG[0]} RGLB=${RG[1]} RGPL=${RG[2]} RGPU=${RG[3]} RGSM=${RG[4]} RGCN=${RG[5]} RGDS=${RG[6]} RGDT=${RG[7]} RGPI=${RG[8]} 

    fi

    let COUNTER+=1

    ### Mark duplicates.
###
    java -jar $PICARD/MarkDuplicates.jar INPUT=$PROJECT_DIR/rg/$FILE_NAME\_rg.sam OUTPUT=$PROJECT_DIR/dups/$FILE_NAME\_rg_dups.sam M=$PROJECT_DIR/dups/$FILE_NAME\_rg_dups.metrics CREATE_INDEX=TRUE
    
    ### Count gene set/alignment overlap.
###
    htseq-count -s $STRANDED -i gene_name $PROJECT_DIR/dups/$FILE_NAME\_rg_dups.sam $GTF > $PROJECT_DIR/htseq/$FILE_NAME.counts

done

### Turn counts files in to matrix of counts for DESeq.
### JOIN_COMM is the join command that will be created for however many files are in the htseq directory.

COUNTS=($PROJECT_DIR/htseq/*.counts)

START=0
END=$(expr ${#COUNTS[@]})

JOIN_COMM="join ${COUNTS[$START]} ${COUNTS[$START+1]} "
echo -n " " > $PROJECT_DIR/htseq/$PROJ_NAME.header
echo -n $(basename ${COUNTS[$START]} .counts) >> $PROJECT_DIR/htseq/$PROJ_NAME.header
echo -n " " >> $PROJECT_DIR/htseq/$PROJ_NAME.header
echo -n $(basename ${COUNTS[$START+1]} .counts) >> $PROJECT_DIR/htseq/$PROJ_NAME.header
let START+=2

for (( i=$START; i<=$END; i++ ))
do

    if [ $i == $END ]
    then
	JOIN_COMM="$JOIN_COMM > $PROJECT_DIR/htseq/$PROJ_NAME.tmp"
    else
	JOIN_COMM="$JOIN_COMM| join - ${COUNTS[$i]} "
	echo -n " " >> $PROJECT_DIR/htseq/$PROJ_NAME.header
	echo -n $(basename ${COUNTS[$i]} .counts) >> $PROJECT_DIR/htseq/$PROJ_NAME.header
    fi
done

echo -n "\n" >> $PROJECT_DIR/htseq/$PROJ_NAME.header

eval $JOIN_COMM
tail -5 $PROJECT_DIR/htseq/$PROJ_NAME.tmp > $PROJECT_DIR/htseq/$PROJ_NAME.stats

cat $PROJECT_DIR/htseq/$PROJ_NAME.header > $PROJECT_DIR/htseq/$PROJ_NAME.txt
head -n -5 $PROJECT_DIR/htseq/$PROJ_NAME.tmp >> $PROJECT_DIR/htseq/$PROJ_NAME.txt

rm $PROJECT_DIR/htseq/$PROJ_NAME.tmp  

### Move $PROJECT_DIR/htseq/$PROJ_NAME.txt to RStudio via WinSCP.
### I prefer to not do this part automatically, but that doesn't mean it couldn't be done.  
