rule trimming:
    input:
        fwd = "samples/raw/{sample}_R1.fq",
        rev = "samples/raw/{sample}_R2.fq"
    output:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq",
        single = "samples/trimmed/{sample}_R1_singletons.fq"
    message:
        """--- Trimming."""
    shell:
        """sickle pe -f {input.fwd} -r {input.rev}  -l 40 -q 20 -t sanger  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log"""

rule fastqc:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        fwd = "samples/fastqc/{sample}/{sample}_R1_t_fastqc.zip",
        rev = "samples/fastqc/{sample}/{sample}_R2_t_fastqc.zip"
    log:
        "logs/fastqc/{sample}_fastqc.log"
    conda:
        "../envs/omic_qc_wf.yaml"
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """fastqc --outdir  samples/fastqc/{wildcards.sample} --extract  -f fastq {input.fwd} {input.rev}"""


rule STAR:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab"
    threads: 12
    params:
        gtf=config["gtf_file"]
    log:
        "logs/star/{sample}_star.log"
    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix samples/star/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --twopassMode Basic
                """)

rule star_statistics:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab"
    params:
        project_id = config["project_id"]
    output:
        "results/tables/{params.project_id}_STAR_mapping_statistics.txt"
    message:
        "Executing star_statistics on {timestamp}"
    log:
        "logs/star_statistics/"
    shell:
        "python StarUtilities.py -d samples/star -p {params.project_id}"


rule picard:
  input:
      "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
  output:
      temp("samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam")
  params:
      name="rmd_{sample}",
      mem="5300"
  threads: 1
  run:
    picard=config["picard_tool"]

    shell("java -Xmx3g -jar {picard} \
    INPUT={input} \
    OUTPUT={output} \
    METRICS_FILE=samples/genecounts_rmdp/{wildcards.sample}_bam/{wildcards.sample}.rmd.metrics.text \
    REMOVE_DUPLICATES=true")


rule sort:
  input:
    "samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam"
  output:
    "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
  params:
    name = "sort_{sample}",
    mem = "6400"
  conda:
    "../envs/omic_qc_wf.yaml"
  shell:
    """samtools sort -O bam -n {input} -o {output}"""


rule samtools_stats:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        "samples/samtools_stats/{sample}.txt"
    log:
        "logs/samtools_stats/{sample}_samtools_stats.log"
    conda:
        "../envs/omic_qc_wf.yaml"
    wrapper:
        "0.17.0/bio/samtools/stats"


rule genecount:
    input:
        "samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam"
    output:
        "samples/htseq_count/{sample}_htseq_gene_count.txt",
    log:
        "logs/genecount/{sample}_genecount.log"
    params:
        name = "genecount_{sample}",
        gtf = config["gtf_file"]
    conda:
        "../envs/omic_qc_wf.yaml"
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
