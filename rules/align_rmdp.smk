rule trimming:
    input:
        fwd = "samples/raw/{smp}_R1.fq.gz",
        rev = "samples/raw/{smp}_R2.fq.gz"
    output:
        fwd = "samples/trimmed/{smp}_R1_t.fq",
        rev = "samples/trimmed/{smp}_R2_t.fq",
        single = "samples/trimmed/{smp}_R1_singletons.fq"
    message:
        """--- Trimming."""
    shell:
        """sickle pe -f {input.fwd} -r {input.rev}  -l 40 -q 20 -t sanger  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log"""

rule fastqc:
    input:
        fwd = "samples/trimmed/{smp}_R1_t.fq",
        rev = "samples/trimmed/{smp}_R2_t.fq"
    output:
        fwd = "samples/fastqc/{smp}/{smp}_R1_t_fastqc.zip",
        rev = "samples/fastqc/{smp}/{smp}_R2_t_fastqc.zip"
    log:
        "logs/fastqc/{smp}_fastqc.log"
    conda:
        "envs/omic_qc_wf.yaml"
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """fastqc --outdir  samples/fastqc/{wildcards.smp} --extract  -f fastq {input.fwd} {input.rev}"""


rule STAR:
    input:
        fwd = "samples/raw/{smp}_R1_t.fq",
        rev = "samples/raw/{smp}_R2_t.fq"
    output:
        "samples/star/{smp}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{smp}_bam/ReadsPerGene.out.tab",
    threads: 12
    params:
        gtf=GTF,
    log:
        "logs/star/{smp}_star.log"
    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix samples/star2/{wildcards.smp}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --twopassMode Basic
                """)

rule star_statistics:
    input:
        rules.STAR.output
    params:
        project_id = project_id
    output:
        "results/tables/{}_STAR_mapping_statistics.txt".format(config['project_id'])
    message:
        "Executing star_statistics on {timestamp}"
    log:
        "logs/star_statistics/"
    shell:
        "python StarUtilities.py -d samples/star -p {params.project_id}"


rule picard:
  input:
      "samples/star/{smp}_bam/Aligned.sortedByCoord.out.bam"
  output:
      temp("samples/genecounts_rmdp/{smp}_bam/{smp}.rmd.bam")
  params:
      name="rmd_{smp}",
      mem="5300"
  threads: 1
  run:
    picard=config["picard_tool"]

    shell("java -Xmx3g -jar {picard} \
    INPUT={input} \
    OUTPUT={output} \
    METRICS_FILE=samples/genecounts_rmdp/{wildcards.smp}_bam/{wildcards.smp}.rmd.metrics.text \
    REMOVE_DUPLICATES=true")


rule sort:
  input:
    "samples/genecounts_rmdp/{smp}_bam/{smp}.rmd.bam"
  output:
    "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam"
  params:
    name = "sort_{smp}",
    mem = "6400"
  conda:
    "envs/omic_qc_wf.yaml"
  shell:
    """samtools sort -O bam -n {input} -o {output}"""


rule samtools_stats:
    input:
        "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam"
    output:
        "samples/samtools_stats/{smp}.txt"
    log:
        "logs/samtools_stats/{smp}_samtools_stats.log"
    conda:
        "envs/omic_qc_wf.yaml"
    wrapper:
        "0.17.0/bio/samtools/stats"


rule genecount:
    input:
        "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam"
    output:
        "samples/htseq_count/{smp}_htseq_gene_count.txt",
    log:
        "logs/genecount/{smp}_genecount.log"
    params:
        name = "genecount_{smp}",
        gtf = GTF
    conda:
        "envs/omic_qc_wf.yaml"
    threads: 1
    shell:
        """
          htseq-count \
                -f bam \
                -r name \
                -s reverse \
                -m union \
                {input} \
                {params.gtf} > {output}"""