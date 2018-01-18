#!/bin/bash

PROJ="/home/users/letaw/lustre1/projs/beataml2"
export PATH="/opt/installed"

echo $PATH
/mnt/lustre1/CompBio/bin/rsem-prepare-reference --transcript-to-gene-map $PROJ/gene_trans.tsv $PROJ/known_ensembl_trans.fasta $PROJ/hg19_2
