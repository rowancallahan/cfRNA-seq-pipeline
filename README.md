Omics-QC-pipeline
======================

Pipeline to assess quality of Omics datasets

This is a package of Python and R scripts that enable reading, processing and analysis of Omics' datasets. 
This package implements the Snakemake management workflow system and is currently implemented to work with 
the cluster management and job scheduling system SLURM. 

Questions/issues
======================

Please add an issue to the Omics-QC-pipeline repository. We would appreciate if your issue included sample code/files 
(as appropriate) so that we can reproduce your bug/issue. 


Contributing
======================

We welcome contributors! For your pull requests, please include the following:

* Sample code/file that reproducibly causes the bug/issue
* Documented code providing fix
* Unit tests evaluating added/modified methods. 

Workflow
======================

Locate raw files:
* After your files are sequenced, they are placed in /home/groups/CEDAR/seq/library_name.
* Change into this directory and find out details about the fastq files contained inside.

```
$ cd /path/to/raw/data
$ ls -alh
```

Check md5sum.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Move your files into the archive to be stored.

```
$ mv /path/to/raw/data /path/to/archive
```

Check md5sum again to ensure your sequencing files are not corrupted.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Unzip all fastq files.

```
$ gunzip –d sample.fastq.gz
$ ctrl+z
$ bg
```

Clone the Omics-QC Pipeline into your working directory.

```
$ git clone https://github.com/ohsu-cedar-comp-hub/Omics-QC-pipeline.git
```

Create a sample/raw directory, a logs directory and a data directory (if they do not exist) in your wdir().

```
$ mkdir logs
$ mkdir data
$ mkdir samples
$ cd samples
$ mkdir raw
```

Symbollically link the fastq files of your samples to the wdir/samples/raw directory using a simple bash script loop in your terminal.

```
$ ls -1 /path/to/data/LIB*R1*fastq | while read fastq; do
    R1=$( basename $fastq | cut -d _ -f 2 | awk '{print $1"_R1.fq"}' )
    R2=$( basename $fastq | cut -d _ -f 2 | awk '{print $1"_R2.fq"}' )
    echo $R1 : $R2
    ln -s $fastq ./$R1
    ln -s ${fastq%R1_001.fastq}R2_001.fastq ./$R2
done
```

Upload your metadata file to the /data directory, with the correct formatting:
* Columns should read:
```StudyID   SampleID   Type   Plasma_volume   RNA_volume   RNA_extracted_by   RNA_extraction_date   Lib_prep_by   Lib._Conc.   Sample_Information   Notes```
* Each row should be a sample, with the above information provided
* All values in this file should be tab-separated

Edit the omic_config.yaml in your wdir():
* Change the project_id to a unique project identifier
* Add appropriate contrasts based on your samples under the [diffexp][contrasts] section
* Add the path to your metadata file for the omic_meta_data and samples parameters

Do a dry-run of snakemake to ensure proper execution before submitting it to the cluster (in your wdir).

```
$ snakemake -np --verbose
```

Once your files are symbolically linked, you can submit the job to exacloud via your terminal window.

```
$ sbatch submit_snakemake.sh
```

To see how the job is running, look at your queue.

```
$ squeue -u your_username
```

Detailed Steps and Output
=================================

Alignment
======================
1) Trimming
    * Trimming of paired-end reads was performed using the trimming tool sickle
    * The output is located in `samples/trimmed/`
2) Quality Analysis
    * Trimmed reads were subject to fastqc quality analysis
    * The output is located in `samples/fastqc/{sample}/{samples}_t_fastqc.zip`
3) Alignment
    * Trimmed reads were aligned to the hg19 genome assembly using STAR
        * We included a two pass mode flag in order to increase the number of aligned reads
        * Output is placed in `samples/star/{sample}_bam/`
            * Output directory includes: `Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab`, and `Log.final.out`
    * We extracted the statistics from the STAR run, and placed them in a table, summarizing the results across all samples from the Log.final.out output of STAR
        * Output is `results/tables/{project_id}_STAR_mapping_statistics.txt`
4) Summarizing output
    * htseq and samtools were used to extract the gene counts for each sample
    * We summarize these results into 1 table, which includes the gene counts across all samples
    * The output is located in `data/{project_id}_counts.txt`
