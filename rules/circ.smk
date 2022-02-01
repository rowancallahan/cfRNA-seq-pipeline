rule bwa: 
        input: "samples/raw/{sample}_R1.fq" 
        output: "data/sam/bwa_map_{sample}.sam"  

        params: name="bwa_{sample}", partition="exacloud", mem="64000" 
        threads: 12  

        run: 
                bwa="/home/groups/CEDAR/tools/bwa" 
                ref_fa="/home/groups/CEDAR/anurpa/genomes/GRCh38.primary_assembly.genome.fa" 

                shell (""" 
                        {bwa}  mem -T 12 {ref_fa} {input[0]}  > {output}

                        """)

rule ciri2:
        input:"data/sam/bwa_map_{sample}.sam"
        output:"ciri2_output/ciri2out_bwa_map_{sample}.txt"

        params:name="ciri2_{sample}", partition="long_jobs", mem="6000"
        threads:4

        run:
                ciri2="/home/groups/CEDAR/tools/CIRI2.pl"
                ref_fa="/home/groups/CEDAR/anurpa/genomes/GRCh38.primary_assembly.genome.fa"
                ref_anno="/home/groups/CEDAR/callahro/cfrna_references/combined_gencode_v27_ercc.gtf"  

                shell ("""
                        perl {ciri2} -T 12 -I {input} -O {output} -F {ref_fa} -A {ref_anno}

                        """) 
