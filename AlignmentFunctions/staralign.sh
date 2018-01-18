#!/bin/bash

### Initialize global variables.
### Can set default values at this time.

gtf_file=""
index_dir=""
project_dir=""
out_dir="./"
sample_files=()
num_files=0
threads=15
paired=false
###  is this a strand-specific assay?  It could be yes, no, or reverse.  See HTSeq-count help for more.
stranded=no
reset=false
read_files_command=zcat
zipped=true
sample_name=""
here="$(pwd)"

#source activate star

function main ()
{
    #echo "in main"

    parseArgs "$@"

    echo "gtf_file: $gtf_file"
    echo "index_dir $index_dir"
    #echo "project_dir: $project_dir"
    echo "out_dir: $out_dir"
    echo "sample_name: $sample_name"
    echo "sample_files: ${sample_files[@]}"
    echo "threads: $threads"
    echo "paired: $paired"
    echo "stranded: $stranded"
    echo "oSsF: $oSsF"
    echo "oFIM: $oFIM" 
    echo "zipped $zipped"
    echo "read_files_command: $read_files_command"
    echo "use_env: $use_env"
    #echo "reset: $reset"
    working_dir=${out_dir}/${sample_name}
    echo "working_dir: $working_dir"


    #exit

    ### set up filenames for paired or unpaired reads
    # variable will contain space-delimited list if more than one file matches
    # exactly 2 files per sample are recommended for paired-end reads
    # STAR takes space-separated pairs of fastq files for paired-end reads; split files are comma-delimited
    
    #SAMP_FASTQ=`ls $PROJECT_DIR/$raw/$SAMP_NAME*$fastq`

    ### Run STAR
    #if [ ! -d ${project_dir}/${output_dir}/${sample_name} ]; then
    #    echo "making project directory" ${project_dir}/${output_dir}/${sample_name}
    #    mkdir ${project_dir}/${output_dir}/${sample_name}
    #    cd ${project_dir}/${output_dir}/${sample_name}
    #fi

    mkdir -p $working_dir
    cd $working_dir

    echo "running STAR from here: $(pwd)"

    ### Comment from original:
    ### consider adding wiggle-track output? raw or norm, read 5' or all, stranded.. STAR can also generate wiggle or bedgraph from position-sorted BAM output
    cmd=" 
    STAR                                        
        --genomeDir $index_dir                  
        --readFilesCommand ${read_files_command}  
        --outFileNamePrefix ${sample_name}        
        --readFilesIn "${sample_files[@]}"      
        --runThreadN $threads                   
        --twopassMode Basic                     
        --quantMode GeneCounts                  
        --outReadsUnmapped Fastq                
        --outSAMstrandField $oSsF               
        --outFilterIntronMotifs $oFIM           
        --outSAMtype BAM SortedByCoordinate     
        --limitBAMsortRAM 34400000000           
        --bamRemoveDuplicatesType UniqueIdentical
    "
    echo "$cmd"
    cmd=$(echo "$cmd" | tr "\n" " " | tr -s " ")
    echo "$cmd"
    eval "$cmd"

    cd $here
}

function showUsage()
{
    #echo "in showUsage"

    local exitcode=$1
    if [[ ! -z $exitcode ]]; then
        exitcode=0
    fi
    echo "usage: starnaseq [ -h ] -g GTF_FILE -d PROJ_DIR -f FASTQ_DIR [ -o OUT_DIR ] [ -t FASTQ_TAG ] [ -r ] [ -p ] [ -s (yes|no|reverse) ]

    Help:
    ----
    -h              Print this help message.
    
    Required arguments:
    ------------------
    -g GTF_FILE     The GTF file for the genome.

    -i INDEX_DIR    Directory containing the genome
                    index file.

    -f FASTQ_FILE(s)
                    Sample fastq files. For paired-end reads, two files
                    are requried. To supply these, use the -f option 
                    twice. If two files are given and the paired ends
                    flag is not provided, then an error will be reported.

    Optional arguments:
    ------------------
    -e ENV          Use the a conda environment ENV. This environment will
                    be activated prior to running STAR.

    -n SAMPLE_NAME  Name to use for the samples. Gets added to the file prefix
                    option of STAR.

    -o OUT_DIR      Output directory for the alignments. This gets
                    passed to the file prefix option of STAR.
                    output directories will be emptied prior to running.

    -p              Flag to indicated paired ends. Default is not paired.
    
    -r              Flag to trigger resetting project directory. All

    -s STRANDED     Is this a strand-specific assay?  Acceptable values are
                    (yes|no|reverse). See HTSeq-count help for more. Default
                    is \"no\".

    -t THREADS      Number of computing threads. Defaults to 15--assuming an 
                    8-core processor with dual threads.

    -z ZIPPED       Boolean value (true|false) stating whether or not
                    the fastq files are gzipped. Default is true.
    "
    exit $exitcode
}

function parseArgs()
{
    #echo "in parseArgs"
    numfiles=0

    while getopts ":hpre:f:g:i:n:o:s:t:z:" opt; do
        #echo "parsing arg $opt"
        case "$opt" in
            h)
                showUsage 0
                ;;
            p)
                paired=true
                #echo "setting paired to true"
                ;;
            r)
                reset=true
                #echo "setting reset to true"
                ;;
            e) 
                use_env="$OPTARG"
                #echo "use_env: $use_env"
                ;;
            f)  ### required
                #num_files=${#sample_files[@]}
                #echo "$num_files"
                num_files=$((num_files+1))
                sample_files[$num_files]="$OPTARG"
                #echo "sample_files: ${sample_files[@]}"
                ;;
            g)  ### required
                gtf_file="$OPTARG"
                #echo "gtf_file: $gtf_file"
                ;;
            i)  ### required
                index_dir="$OPTARG"
                #echo "index_dir: $index_dir"
                ;;
            n)
                sample_name="$OPTARG"
                #echo "sample_name: $name"
                ;;
            o)
                out_dir="$OPTARG"
                #echo "out_dir: $out_dir"
                ;;
            s)
                stranded="$OPTARG"
                #echo "stranded: $stranded"
                ;;
            t)
                threads="$OPTARG"
                #echo "threads: $threads"
                ;;           
            z)
                zipped="$OPTARG"
                #echo "stranded: $stranded"
                ;;
            ?)
                echo "Error: did not recognize option $opt"
                showUsage 1
                ;;
        esac
    done
    shift $((OPTIND-1))

    if  [[ "$gtf_file" == "" ]] \
        || \
        [[ "$index_dir" == "" ]] \
        || \
        [[ $num_files == 0 ]]
    then
        echo
        echo "ERROR: Missing a required argument [g|i|d|f]"
        echo
        showUsage 1
    fi

    if [[ $(echo "yes no reverse" | grep $stranded) == "" ]]; then
        echo
        echo "ERROR: Acceptable stranded values (-s) are (yes|no|reverse)."
        echo 
        showUsage 1
    fi

    case $paired in
        true)
            if [[ $num_files < 2 ]]; then
                echo
                echo "ERROR: Paired selected, but only one sample file provided."
                echo
                showUsage 1
            fi
            ;;
        false)
            if [[ $num_files > 1 ]]; then
                echo
                echo "ERROR: Paired end NOT selected. Only one sample file can be provided."
                echo
                showUsage 1
            fi
            ;;
    esac

    if [ "$stranded" == "no" ]; then
        oSsF="intronMotif"; 
        oFIM="RemoveNoncanonicalUnannotated"
    else
        oSsF="None"; 
        oFIM="None"
    fi

    if [[ "$zipped" == false ]]; then
        read_files_command=cat
    fi
    
    out_dir=${out_dir%%/}
    #echo "done parsing"
}

function clearProject ()
{
    rmdir  ${project_dir}/{"$out",_rg,_dups,_htseq}
    exit 0
}

#echo "$@"
#exit

#showUsage 0

main "$@"
