
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
    name="star_{sample}",
    mem="64000"
  run:
    STAR=config["star_tool"]
    pathToGenomeIndex = config["star_index"]

    shell("""
    {STAR} --runThreadN {threads} \
           --runMode alignReads \
           --genomeDir {pathToGenomeIndex} \
           --readFilesIn {input[0]} {input[1]} \
           --outFileNamePrefix samples/circexplorer/{wildcards.sample}_chim_bam/  \
           --outSAMtype BAM SortedByCoordinate \
           --chimSegmentMin 5 \
           --chimJunctionOverhangMin 5
    """)

rule circ_explorer:
  input:
    "samples/circexplorer/{sample}_chim_bam/Chimeric.out.junction"
  output:
    "results/circexplorer/{sample}_circrna/circularRNA_known.txt"

  params:
    mem="64000",
    refflat= config["refflat2"],
    genome= config["genome"]
  conda:
    "../envs/circexplorer.yaml"

  shell:"""
    fast_circ.py parse \
        -r {params.refflat} \
        -g {params.genome} \
        -t STAR \
        -o results/circexplorer/{wildcards.sample}_circrna \
        {input}
    """

rule junction_counts:
  input:
    expand("results/circexplorer/{sample}_circrna/circularRNA_known.txt",sample=SAMPLES)
  output:
    "results/tables/{project_id}_circexplorer_junctioncounts.txt".format(project_id=config['project_id'])
  params:
    project_id=config["project_id"]
  conda:
    "../envs/junction_counts.yaml"
  script:
    "../scripts/compile_splice_junctions.R"

