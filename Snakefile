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
insertion_and_clipping_prof_ext = ['r','R1.pdf','R2.pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']
#EXT=['circrna/circularRNA_known.txt']


# TODO generate initializing rule to automatically generate log out for all rules

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
    return 1 if len(config["samples"]) < 100 or few_coeffs else 6


def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]


for sample in SAMPLES:
    message("Sample " + sample + " will be processed")


rule all:
    input:
        expand("samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        expand("samples/star/{sample}_bam/ReadsPerGene.out.tab", sample = SAMPLES),
        "results/tables/{}_STAR_mapping_statistics.txt".format(config['project_id']),
        expand("samples/fastqc/{sample}/{sample}_{ext}_t_fastqc.zip", sample = SAMPLES, ext = ext),
        expand("samples/genecounts_rmdp/{sample}_bam/{sample}.rmd.bam", sample = SAMPLES),
        expand("samples/genecounts_rmdp/{sample}_bam/{sample}_sort.rmd.bam", sample = SAMPLES),
        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        expand("rseqc/inner_distance/{sample}/{sample}.inner_distance{ext}", sample = SAMPLES, ext = inner_distance_ext),
        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt", sample=SAMPLES),
        "results/tables/{}_Normed_with_Ratio_and_Abundance.txt".format(config['project_id']),
        "results/diffexp/pca.pdf",
        expand("results/diffexp/{project_id}_all.rds",project_id = config['project_id']),


include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"
#include: "rules/circ_explorer.smk"
