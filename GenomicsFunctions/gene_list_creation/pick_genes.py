#!/usr/bin/python

# NOTE: Run crush_exons.py on exons from BioMart.

# In BioMart (grch37.ensembl.org/biomart)
# Select hsa, under filter:
# chroms 1-22,X,Y
# with HGNC symbols only
# protein_coding
# under attributes:
# chromsome name
# gene start
# gene end
# ensembl gene id
# ensembl transcript id
# transcript length
# strand
# ccds id
# RefSeq mRNA
# HGNC symbol

# After running through script, look for duplicate id's in resulting list.
# The following modifications are made:

# Duplicate ENST
# There are many duplicate ENST id's when looking at the larger list, mainly in MIR's.  We will ignore these duplicates for the time being.

# SEBOX and VTN same thing (remove SEBOX, not in kc)

# Duplicate CCDS

# CCDS10162 (GCOM1 and MYZAP, remove MYZAP, no RefSeq)
# CCDS11304 (CCL15-CCL14 and CCL15, remove CCL15-CCL14, no RefSeq)
# CCDS114 (APITD1 and APITD1-CORT, remove APITD1-CORT, no RefSeq)
# CCDS12649 (APOC4 and APOC4-APOC2, remove APOC4-APOC2, no RefSeq)
# CCDS13424 (TMEM189-UBE2V1 and TMEM189, remove TMEM-UBE2V1, shorter)
# CCDS32682 (NME1-NME2 and NME2, remove NME2, no RefSeq)
# CCDS33030 (RAB4B and RAB4B-EGLN2, remove RAB4B-EGLN2, no RefSeq)
# CCDS4224 (ANKHD1 and ANKHD1-EIF4EBP3, remove ANKHD1-EIF4EBP3, no RefSeq)
# CCDS54887 (TMED7-TICAM2 and TICAM2, remove TMED7-TICAM2, no RefSeq)
# CCDS60529 (PGBD3 and ERCC6-PGBD3, keep ERCC6-PGBD3, it is longer, remove from kc)
# CCDS7542 (C10orf32 and C10orf32-ASMT, remove C10orf32-ASMT, no RefSeq)

# Duplicate RefSeq

# NM_001258399 (CCDC103 and FAM187A, remove FAM187A, no CCDS, not in kc)
# NM_001242493 (FAM156B and FAM156A, different coordinates, not sure how to proceed here, perhaps remove RefSeq ID? (Removed RefSeq ID for FAM156B, but kept gene in list)

# Special case: Duplicates only in larger list including not just protein coding genes.
# NM_004678 (BPY2, BPY2C, BPY2B) (Remove RefSeq ID from 2C and 2B, the CCDS ID's are all unique for these 3)
# NR_003334 (SNORD116-20, SNORD116-21) (Remove SNORD116-21, these entries are completely identical except HGNC symbol)

# CCDS196 (neither NBL1 or MINOS1-NBL1 are in curated list, remove CCDS id from MINOS1-NBL1)
# CCDS2041 (neither MRPL30 or C2orf15 in curated list, remove CCDS id from C2orf15)
# CCDS6195 (neither C8orf44-SGK3 or SGK3 in curated list, remove CCDS id from C8orf44-SGK3)

# Look for divergence between the kdl_curated_list and this list.
# Take divergent list and run back through BioMart, but this time with no protein_coding option set.
# Add this list to the previous list, run through this script again.
# Look for divergence once again and duplicate id's.  Send divergent list to Andrew, or someone at KDL to deal with.

# We are creating a larger list as well.  This list is made by not restricting to protein coding in the first step.  Duplicates should be resolved here as well, for consistency.

import sys

handle_genes = open(sys.argv[1], 'rU')

fgenes = {}
BLACKLIST = ["SEBOX", "MYZAP", "CCL15-CCL14", "APITD1-CORT", "APOC4-APOC2", "TMEM189-UBE2V1", "NME2", "RAB4B-EGLN2", "ANKHD1-EIF4EBP3", "TMED7-TICAM2", "C10ORF32-ASMT", "FAM187A", "PGBD3", "SNORD116-21", "FAM226B", "MIR1184-3", "MIR1302-10", "MIR1302-11", "MIR1302-9", "MIR29B1", "MIR3116-2", "MIR3118-2", "MIR3118-3", "MIR3119-2", "MIR3130-2", "MIR3150B", "MIR3158-2", "MIR3160-2", "MIR3191", "MIR3199-2", "MIR3202-2", "MIR3622B", "MIR3680-2", "MIR3688-2", "MIR3910-2", "MIR3913-2", "MIR3926-2", "MIR4283-2", "MIR4436B2", "MIR4444-2", "MIR4477B", "MIR4509-2", "MIR4509-3", "MIR451A", "MIR451B", "MIR4732", "MIR4520B", "MIR4524B", "MIR4659B", "MIR4662B", "MIR4679-2", "MIR4771-2", "MIR4773-2", "MIR4776-2", "MIR5692A2", "MIR5701-2", "RNA5SP389", "MIR3118-6", "MIR548D2", "MIR338", "MIR133A1", "MIR4435-2", "MIR548D1"]


with handle_genes as genes:
    next(genes)
    for line in genes:
        line = line.rstrip('\n').split('\t')
        
        chrom = line[0]
        start = line[1]
        stop = line[2]
        ensg = line[3]
        enst = line[4]
        tlen = int(line[5])
        strand = line[6]
        ccds = line[7]
        refseq = line[8]
        hgnc = line[9]

        if ccds and refseq:
            if hgnc not in fgenes:
                fgenes[hgnc] = [chrom, start, stop, ensg, enst, tlen, strand, ccds, refseq, hgnc]
            else:
                if tlen > fgenes[hgnc][5]:
                    fgenes[hgnc] = [chrom, start, stop, ensg, enst, tlen, strand, ccds, refseq, hgnc]

handle_genes = open(sys.argv[1], 'rU')

with handle_genes as genes:
    next(genes)
    for line in genes:
        line = line.rstrip('\n').split('\t')
        
        chrom = line[0]
        start = line[1]
        stop = line[2]
        ensg = line[3]
        enst = line[4]
        tlen = int(line[5])
        strand = line[6]
        ccds = line[7]
        refseq = line[8]
        hgnc = line[9]
        
        if ccds or refseq:
            if hgnc not in fgenes:
                fgenes[hgnc] = [chrom, start, stop, ensg, enst, tlen, strand, ccds, refseq, hgnc]
            else:
                if (tlen > fgenes[hgnc][5] and fgenes[hgnc][7] == "") or (tlen > fgenes[hgnc][5] and fgenes[hgnc][8] == ""):
                    fgenes[hgnc] = [chrom, start, stop, ensg, enst, tlen, strand, ccds, refseq, hgnc]

handle_genes = open(sys.argv[1], 'rU')

with handle_genes as genes:
    next(genes)
    for line in genes:
        line = line.rstrip('\n').split('\t')

        chrom = line[0]
        start = line[1]
        stop = line[2]
        ensg = line[3]
        enst = line[4]
        tlen = int(line[5])
        strand = line[6]
        ccds = line[7]
        refseq = line[8]
        hgnc = line[9]
        
        if not ccds and not refseq:
            if hgnc not in fgenes:
                fgenes[hgnc] = [chrom, start, stop, ensg, enst, tlen, strand, ccds, refseq, hgnc]
            else:
                if tlen > fgenes[hgnc][5] and (fgenes[hgnc][7] == "") and (fgenes[hgnc][8] == ""):
                    fgenes[hgnc] = [chrom, start, stop, ensg, enst, tlen, strand, ccds, refseq, hgnc]

handle_out = open(sys.argv[2], 'w')

for key in sorted(fgenes):
    if key.upper() not in BLACKLIST:
        if fgenes[key][6] == "1":
            handle_out.write(fgenes[key][0] + '\t' + fgenes[key][1] + '\t' + fgenes[key][2] + '\t' + fgenes[key][3] + \
                                 '\t' + fgenes[key][4] + '\t+\t' + fgenes[key][7] + '\t' + fgenes[key][8] + '\t' + fgenes[key][9] + '\n')
        elif fgenes[key][6] == "-1":
            handle_out.write(fgenes[key][0] + '\t' + fgenes[key][1] + '\t' + fgenes[key][2] + '\t' + fgenes[key][3] + \
                                 '\t' + fgenes[key][4] + '\t-\t' + fgenes[key][7] + '\t' + fgenes[key][8] + '\t' + fgenes[key][9] + '\n')
        else:
            raise Exception("Strand should be either 1 or -1")

handle_out.close()
