rule circ_star:
  input:
    fwd = "samples/trimmed/{sample}_R1_t.fq",
    rev = "samples/trimmed/{sample}_R2_t.fq"
  output:
    "samples/circexplorer/{sample}_chim_bam/Chimeric.out.junction"
  threads: 12
  params:
    name="star_{sample}", mem="64000"
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

