#!/usr/bin/bash
#SBATCH --time 35:00:00
#SBATCH --partition exacloud
#SBATCH --job-name omic-wf-submit

snakemake --use-conda -j --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem {cluster.mem}" -s Snakefile