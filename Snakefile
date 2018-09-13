__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

configfile:"omic_config.yaml"
import datetime
import sys
import os

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

project_id = config["project_id"]
rseqqc_env = config["rseqc_env"]
samples, = glob_wildcards("/home/exacloud/lustre1/CEDAR/cfrna/analysis/pancan/genecounts_rmdp/{sample}_bam/")

#TODO Remove abspath for glob_wildcards function

ext = ['r','R1.pdf','R2.pdf','xls']
insertion_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']

rseq_func=['insertion_profile','read_distribution','clipping_profile']
rseq_dir=['insertion_profile','read_distribution','clipping_profile']

rule_dirs = ['insertion_profile','inner_distance','clipping_profile','read_GC','star_statistics','compile_counts','generate_qc_qa','run_qc_qa','star_statistics','deseq2','fastqc']
for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
  sys.stderr.write("|--- " + mes + "\n")


def format_plot_columns():
    factors = config['meta_columns_to_plot'].keys()
    reformat_factors = '"' + '","'.join(factors) + '"'
    return 'c({})'.format(reformat_factors)


def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]


rule all:
    input:
        "results/tables/{}_STAR_mapping_statistics.txt".format(config['project_id']),
        "data/{}_counts.txt".format(config['project_id']),
        "analysis_code/{}_analysis.R".format(config['project_id']),
        "results/tables/{}_Normed_with_Ratio_and_Abundance.txt".format(config['project_id']),
        expand("samples/trimmed/{smp}_R1_t.fq", smp = samples),
        expand("samples/star/{smp}_bam/Aligned.sortedByCoord.out.bam", smp = samples),
        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=samples, ext=insertion_prof_ext),
        expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = samples, ext = inner_distance_ext),
        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = samples, ext = clipping_prof_ext),
        expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = samples, ext = read_dist_ext),
        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = samples, ext = read_gc_ext),
        expand(["results/diffexp/{contrast}.diffexp.tsv", "results/diffexp/{contrast}.ma-plot.pdf"],contrast = config["diffexp"]["contrasts"]),
        "results/diffexp/pca.pdf"

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


rule STAR:
        input:
            fwd = "samples/trimmed/{smp}_R1_t.fq",
            rev = "samples/trimmed/{smp}_R2_t.fq"
        output:
            temp("samples/star/{smp}_bam/Aligned.sortedByCoord.out.bam")
        threads: 12
        params:
            name="STAR_{smp}",
            mem="64000",
            threads="12"
        run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""mkdir -p samples/star/{wildcards.smp}
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input.fwd} {input.rev} \
                --outFileNamePrefix samples/star/{wildcards.smp}_bam \
                --outSAMtype BAM SortedByCoordinate
                """)

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
    "envs/intergrin.yaml"
  shell:
    """samtools sort -O bam -n {input} -o {output}"""


rule genecount:
  input:
    "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam"
  output:
    "samples/htseq_count/{smp}_htseq_gene_count.txt"
  params:
    name = "genecount_{smp}",
    mem = "5300"
  conda:
    "envs/intergrin.yaml"
  threads: 1
  shell:
    """
      htseq-count \
            -f bam \
            -r name \
            -s reverse \
            -m union \
            --samout=genecounts_rmdp/htseq_samout/Output_{wildcards.smp}.sam \
            {input} \
            /home/exacloud/lustre1/CEDAR/anurpa/genomes/gencode.v27.annotation.gtf > {output}"""


rule insertion_profile:
    input:
        "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam",
    params:
        seq_layout=config['seq_layout'],
    output:
        "rseqc/insertion_profile/{sample}/{sample}.insertion_profile.r",
        "rseqc/insertion_profile/{sample}/{sample}.insertion_profile.R1.pdf",
        "rseqc/insertion_profile/{sample}/{sample}.insertion_profile.R2.pdf",
        "rseqc/insertion_profile/{sample}/{sample}.insertion_profile.xls",
    conda:
        "envs/rseqc.yaml"
    shell:
        "insertion_profile.py -s '{params.seq_layout}' -i {input} -o rseqc/insertion_profile/{wildcards.sample}/{wildcards.sample}"


rule inner_distance:
    input:
        "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam",
    params:
        bed=config['bed_file']
    output:
        "rseqc/inner_distance/{sample}/{sample}.inner_distance.txt",
        "rseqc/inner_distance/{sample}/{sample}.inner_distance_plot.r",
        "rseqc/inner_distance/{sample}/{sample}.inner_distance_plot.pdf",
        "rseqc/inner_distance/{sample}/{sample}.inner_distance_freq.txt",
    conda:
        "envs/rseqc.yaml"
    shell:
        "inner_distance.py -i {input} -o rseqc/inner_distance/{wildcards.sample}/{wildcards.sample} -r {params.bed}"


rule clipping_profile:
    input:
        "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam",
    params:
        seq_layout=config['seq_layout'],
    output:
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.r",
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.R1.pdf",
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.R2.pdf",
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.xls",
    conda:
        "envs/rseqc.yaml"
    shell:
        "clipping_profile.py -i {input} -s '{params.seq_layout}' -o rseqc/clipping_profile/{wildcards.sample}/{wildcards.sample}"


rule read_distribution:
    input:
        "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam",
    params:
        bed=config['bed_file']
    output:
        "rseqc/read_distribution/{sample}/{sample}.read_distribution.txt",
    conda:
        "envs/rseqc.yaml"
    shell:
        "read_distribution.py -i {input} -r {params.bed} > {output}"

rule read_GC:
    input:
        "samples/genecounts_rmdp/{smp}_bam/{smp}_sort.rmd.bam",
    output:
        "rseqc/read_GC/{sample}/{sample}.GC.xls",
        "rseqc/read_GC/{sample}/{sample}.GC_plot.r",
        "rseqc/read_GC/{sample}/{sample}.GC_plot.pdf",
    conda:
        "envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o rseqc/read_GC/{wildcards.sample}/{wildcards.sample}"

rule star_statistics:
    input:
        dir="samples/genecounts_rmdp"
    params:
        project_id = project_id
    output:
        "results/tables/{}_STAR_mapping_statistics.txt".format(config['project_id'])
    message:
        "Executing star_statistics on {timestamp}"
    log:
        "logs/star_statistics/"
    shell:
        "python StarUtilities.py -d {input.dir} -p {params.project_id}"

rule compile_counts:
    input:
        sample_counts="samples/htseq_count/"
    params:
        project_id = project_id
    output:
        "data/{project_id}_counts.txt".format(project_id=project_id)
    shell:
        "python StarUtilities.py -d {input.sample_counts} -p {params.project_id} -c"

rule generate_qc_qa:
 input:
    counts = expand("data/{project_id}_counts.txt", project_id=config['project_id'])
 params:
    project_id = config["project_id"],
    datadir = config['base_dir'],
    meta = config["omic_meta_data"],
    baseline = config["baseline"],
    linear_model = config["linear_model"],
    sample_id = config["sample_id"],
    gtf_file = config["gtf_file"],
    meta_viz = format_plot_columns(),
 output:
    "analysis_code/{project_id}_analysis.R".format(project_id=config['project_id'])
 log:
    "logs/generate_qc_qa/"

 shell:
    "python GenerateAbundanceFile.py -d {params.datadir} -mf {params.meta} -p {params.project_id} -b {params.baseline} -lm {params.linear_model} -id '{params.sample_id}' -pl '{params.meta_viz}' -g '{params.gtf_file}' -df -da {input.counts}"

rule run_qc_qa:
    input:
        script = "analysis_code/{}_analysis.R".format(config['project_id'])
    output:
        "results/tables/{}_Normed_with_Ratio_and_Abundance.txt".format(config['project_id'])
    conda:
        "envs/omic_qc_wf.yaml"
    log:
        "logs/run_qc_qa/"
    shell:
        "Rscript analysis_code/{}_analysis.R".format(config['project_id'])


rule deseq2_init:
    input:
        counts="data/{project_id}_counts.txt".format(project_id = project_id)
    output:
        rds="results/diffexp/{project_id}_all.rds".format(project_id=project_id),
        normed_counts="results/tables/{project_id}_normed_counts.txt".format(project_id = project_id),
        rld_out = "results/diffexp/{project_id}_rlog_dds.rds".format(project_id = project_id),
    params:
        samples=config["samples"],
        design=config["pca"]["labels"],
        row_names=config["sample_id"],
    conda:
        "envs/deseq2.yaml",
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        "scripts/deseq2-init.R"


rule deseq2_plots:
    input:
        rds_object = "results/diffexp/{project_id}_rlog_dds.rds".format(project_id = project_id),
        dds_object = "results/diffexp/{project_id}_all.rds".format(project_id = project_id),

    output:
        pca="results/diffexp/pca.pdf",
        sd_mean_plot="results/diffexp/sd_mean_plot.pdf",
        heatmap_plot = "results/diffexp/heatmap_plot.pdf",
        distance_plot = "results/diffexp/distance_plot.pdf",
        panel_ma = "results/diffexp/panel_ma.pdf",
        var_heat = "results/diffexp/variance_heatmap.pdf",
        ggplot_pca_factor = "results/diffexp/ggplot_factor_pca.pdf",
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "envs/deseq2_plots.yaml"
    log:
        "logs/deseq2/pca.log"
    script:
        "scripts/plot-pca.R"


rule deseq2:
    input:
        rds="results/diffexp/{project_id}_all.rds".format(project_id=project_id)
    output:
        table="results/diffexp/{contrast}.diffexp.tsv",
        ma_plot="results/diffexp/{contrast}.ma_plot.pdf",
        p_hist="results/diffexp/{contrast}.phist_plot.pdf",
    params:
        contrast=get_contrast,
        condition = config["linear_model"]
    conda:
        "envs/deseq2.yaml",
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads()
    script:
        "scripts/deseq2.R"


rule fastqc:
    input:
        fwd = "/home/exacloud/lustre1/CEDAR/cfrna/data/LIB180515JW/{project}_{sample}_{s}_{lane}_R1_001.fastq",
        rev = "/home/exacloud/lustre1/CEDAR/cfrna/data/LIB180515JW/{project}_{sample}_{s}_{lane}_R2_001.fastq"

    output:
        fwd = "fastqc/{sample}/{project}_{sample}_{s}_{lane}_R1_001.fastq_fastqc.zip",
        rev = "fastqc/{sample}/{project}_{sample}_{s}_{lane}_R2_001.fastq_fastqc.zip"

    conda:
        "envs/fastqc-omic-wf.yaml"

    message:
        """--- Quality check of raw data with Fastqc."""

    shell:
        """fastqc --outdir  fastqc/{wildcards.smp} --extract  -f fastq {input.fwd} {input.rev}"""
