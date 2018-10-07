
rule circ_star:
  input:
    fwd = "samples/raw/{sample}_R1.fq",
    rev = "samples/raw/{sample}_R2.fq"
  output:
    "samples/circexplorer/{sample}_chim_bam/Chimeric.out.junction"
  threads: 12
  params:
    STAR=config["star_tool"],
    pathToGenomeIndex = config["star_index"],
    name="star_{sample}", mem="64000"
  run:
    STAR=config["star_tool"],
    pathToGenomeIndex = config["star_index"]
    
    shell("""
    {STAR} --runThreadN {threads} \
           --runMode alignReads \
           --genomeDir {pathToGenomeIndex} \
           --readFilesIn {input.fwd} {input.rev} \
           --outFileNamePrefix samples/circexplorer/{wildcards.sample}_chim_bam/  \
           --outSAMtype BAM SortedByCoordinate \
           --chimSegmentMin 5 \
           --chimJunctionOverhangMin 5
    """)

rule circ_explorer_parse:
  input:
    "samples/circexplorer/{sample}_chim_bam/Chimeric.out.junction"
  output:
    "samples/circexplorer/{sample}_circrna/back_spliced_junction.bed"
  conda:
    "../envs/circexplorer.yaml"

  shell:"""
    CIRCexplorer2 parse -b {output} -t STAR {input}
    """

rule circ_explorer_annotate:
  input:
    "samples/circexplorer/{sample}_circrna/back_spliced_junction.bed"
  output:
    "samples/circexplorer/{sample}_circrna/circularRNA_known.txt"
  params:
    mem="64000",
    refflat= config["refflat2"],
    genome= config["genome"]
  conda:
    "../envs/circexplorer.yaml"

  shell:"""
        CIRCexplorer2 annotate -r {params.refflat} -g {params.genome} -b {input} -o {output}
    """

rule circexplorer_junction_counts:
  input:
    expand("samples/circexplorer/{sample}_circrna/circularRNA_known.txt",sample=SAMPLES)
  output:
    "results/tables/{project_id}_circexplorer_junctioncounts.txt".format(project_id=project_id)
  conda:
    "../envs/junction_counts.yaml"
  script:
    "../scripts/circexplorer_junction_counts.R"
