__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""


import datetime
import sys
import os

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]
rseqqc_env = config["rseqc_env"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq")

ext = ['r','R1.pdf','R2.pdf','xls']
fastq_ext = ['R1_P','R2_P']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']


# TODO generate initializing rule to automatically generate log out for all rules

rule_dirs = ['fastqc','fastqscreen','star','index','bam_statistics','get_bam_coverage','picard','sort','samtools_stats','genecount','count_exons','compile_counts','compile_exon_counts','trimming','insertion_profile','read_distribution','inner_distance','clipping_profile','read_GC','star_statistics','generate_qc_qa','run_qc_qa','star_statistics','deseq2','bwa','ciri','ciri_junction_counts','GO','volcano']
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
    return 1 if len(config["samples"]) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")


rule all:
    input:
        expand("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        expand("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai", sample = SAMPLES),
        expand("samples/star/{sample}_bam/ReadsPerGene.out.tab", sample = SAMPLES),
        expand("samples/star/{sample}_bam/Log.final.out",sample=SAMPLES),
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        expand("samples/bamstats/{sample}/genome_coverage.json", sample = SAMPLES),
        "data/{project_id}_coverage.txt".format(project_id=config["project_id"]),
        expand("samples/fastqc/{sample}/{sample}_{fastq_ext}_t_fastqc.zip", sample = SAMPLES, fastq_ext = fastq_ext),
        expand("samples/fastqscreen/{sample}/{sample}_{fastq_ext}_t_screen.{fastqscreen_ext}", sample=SAMPLES, fastq_ext=fastq_ext, fastqscreen_ext=fastqscreen_ext),
        expand("samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam", sample = SAMPLES),
        expand("samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam", sample = SAMPLES),
        "data/{project_id}_counts.txt".format(project_id=config['project_id']),
        "data/{project_id}_counts_w_stats.txt".format(project_id=config['project_id']),
        "data/{project_id}_exon_counts.txt".format(project_id = config["project_id"]),
        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt", sample=SAMPLES),
        expand("samples/htseq_exon_count/{sample}_htseq_exon_count.txt", sample=SAMPLES),
        "results/tables/{}_Normed_with_Ratio_and_Abundance.txt".format(config['project_id']),
        "results/diffexp/pca.pdf",
        expand("results/diffexp/{project_id}_all.rds",project_id = config['project_id']),
        expand(["results/diffexp/{contrast}.diffexp.tsv", "results/diffexp/{contrast}.ma_plot.pdf","results/diffexp/{contrast}.phist_plot.pdf"],contrast = config["diffexp"]["contrasts"]),
        expand("samples/circexplorer/{sample}_chim_bam/Chimeric.out.junction", sample = SAMPLES),
        expand("samples/ciri/{sample}.sam",sample = SAMPLES),
        expand("results/ciri_out/{sample}_ciriout.txt",sample = SAMPLES),
        "results/tables/{project_id}_ciri_junctioncounts.txt".format(project_id=project_id),
        "results/tables/{project_id}_ciri_frequency.txt".format(project_id=project_id),
        expand(["results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_GO.txt", "results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_GO.txt", "results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01.BP.pdf", "results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01.BP.pdf","results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_classic_5_all.pdf","results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_classic_5_all.pdf"], contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/{contrast}.diffexp.01.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"])

include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"
include: "rules/ciri.smk"
include: "rules/circ_star.smk"
