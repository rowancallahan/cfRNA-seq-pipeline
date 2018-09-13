

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