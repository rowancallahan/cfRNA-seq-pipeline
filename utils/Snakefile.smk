#vim: set syntax=python

__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

configfile: "config.yaml"

import sys
import yaml
	
rule all:
 input:
  "data/{project_id}_counts.txt"

rule qc_qa:
 input:


rule trim_reads:
 input:
  "reads/{sample}.fastq.gz"
 output:
  "reads_trimmed/{sample}.trimmed.fastq"
 log:
  "logs/ctadpt_{sample}.log"
 shell:
  "cutadapt -n 2 -g {config[adapter5]} -a {config[adapter3]} {input} > {output}"
  
rule align_reads:
 input:
  "reads_trimmed/{sample}.trimmed.fastq"
 output:
  "reads_aligned/{sample}.sam"
 log:
  "logs/bwt2_{sample}.log"
 shell:
  "bowtie2 -p 8 -x {config[bowtie2libidx]} --norc -U {input} -S {output} 2> {log}"

rule sam2bam_reads:
 input:
  "reads_aligned/{sample}.sam"
 output:
  "reads_sam2bam/{sample}.bam"
 log:
  "logs/sam2bam_{sample}.log"
 shell:
  "samtools view -bS {input} > {output} 2> {log}"

rule filter_reads:
 input:
  "reads_sam2bam/{sample}.bam"
 output:
  "reads_filtered/{sample}.bam"
 log:
  "logs/filtered_{sample}.log"
 shell:
  "samtools view -bq 10 {input} > {output} 2> {log}"
