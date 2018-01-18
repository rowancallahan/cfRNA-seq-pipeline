#!/usr/bin/env python

### From the BED file describing amplicons given to us by Carol, create the amplicon files
### that will be fed in to the R CNV script.  We will need to utilize samtools faidx to pull GC counts
### for the particular genomic coordinates covered by amplicons.

# BED files supplied by Carol look like: 
# track type=bedDetail ionVersion=4.0 name="Solid_Tumor_IAD32788_Designed" description="Amplicon_Insert"
# chr1	115256501	115256590	AMPL1624980	0	+	.	GENE_ID=NRAS

# amplicons.txt to be fed in to CNV script needs to look like:
# AmpliconId	AmpliconIndex	ChromNum	StartPos	EndPos	Gene	NumGC	Length	GC
# 440774676	1	1	36931917	36932093	CSF3R	115	177	0.650
# 440787931	2	1	36932082	36932236	CSF3R	102	155	0.658
# 440509141	3	1	36932232	36932319	CSF3R	56	88	0.636
# 440524407	4	1	36932307	36932433	CSF3R	79	127	0.622

# AMPL number will be AmpliconId
# Arbitrarily assign index from 1 -> x
# Pull down chromosome, remove chr prefix
# amplicons.txt is not in BED format, while the input file claims to be BED.
# gene can be taken directly from the input file, though we need to verify these gene names all match with CGD.
# use samtools to count GC's, then make the simple calculation

### USAGE: python reformat_amplicons_for_cnv.py <Input BED from Carol> <Output to be used with CNV Rscript> <Location of Ref Genome>

import sys
import subprocess
import string
import re


def parseBedInput(bed_handle):
    """
    Input a BED file as shown above.
    """

    amplicons = []

    with bed_handle as bed:
        next(bed) # Skip header row
        amp = 1
        for amplicon in bed:
            amplicon = amplicon.rstrip('\n').split('\t')
            chrom = sex_chrom_switch(amplicon[0][3:])
            start = int(amplicon[1]) + 1 # BED files are 0-based at the start position, so we need to add 1 here.
            stop = int(amplicon[2])
            amplicon_id = filter(str.isalnum, amplicon[3])
          #   amplicon_id = amplicon[3]
          # #  filter(amplicon_id.isalnum, string.printable)   ### GO HERE NEXT TO DEAL WITh AMPLICON IDS
          #   if '.' in amplicon_id:
          #       amplicon_id = amplicon_id.split('.')[-1]
          #   else:
          #       amplicon_id = amplicon_id[4:]
            gene = amplicon[7].split('=')[1]
            if ';' in gene:  ### This stuff needs to happen because the input files aren't consistent.
                ###Should build a different function to validate format.
                gene = gene.split(';')[0]
                if '_' in gene:
                    gene = gene.split('_')[0]

            amp += 1
            amplicons.append([chrom, start, stop, amplicon_id, gene])
            
    return amplicons


def sex_chrom_switch(chrom):
    """
    Convert X and Y chroms to 23 and 24.
    """
    if chrom == "X":
        return "23"
    elif chrom == "Y":
        return "24"
    else:
        return chrom


def createAmpliconOutput(amplicons, handle_out, ref_path):
    """
    This is the file that is used as input to the CNV Rscript.
    amplicons = [chrom, start, stop, amplicon_id, gene]
    """

    handle_out.write('\t'.join(["AmpliconId", "AmpliconIndex", "ChromNum", "StartPos", "EndPos", "Gene", "NumGC", "Length", "GC\n"]))
    amplicon_index = 0

    for amplicon in amplicons:
        chrom = amplicon[0]
        amp_length = amplicon[2] - amplicon[1] + 1.0
        start = str(amplicon[1])
        stop = str(amplicon[2])
        amplicon_id = amplicon[3]
        gene = amplicon[4]

        amplicon_index += 1

        numgc = countGC(chrom, start, stop, ref_path)
        gc =  "{:.3f}".format(numgc/amp_length)
        handle_out.write('\t'.join([amplicon_id, str(amplicon_index), chrom, start, stop, gene, str(numgc), str(int(amp_length)), str(gc)]))
#        handle_out.write('\t'.join([str(amplicon_index), str(amplicon_index), chrom, start, stop, gene, str(numgc), str(int(amp_length)), str(gc)]))
        handle_out.write('\n')
    

def countGC(chrom, start, stop, ref_path):
    """
    Use "samtools faidx" to find the genomics sequence associated with a particular set of coordinates.
    Return the count of G's and C's in that sequence.
    """

    if chrom == "23":
        chrom = "X"
    if chrom == "24":
        chrom = "Y"
    interval = chrom + ':' + start + '-' + stop
    cmd = ['samtools', 'faidx', ref_path, interval]
    seq = subprocess.check_output(cmd)
    gc_count = (seq.count('G') + seq.count('C'))
    
    return gc_count


def main():

    bed_handle = open(sys.argv[1], 'rU')
    handle_out = open(sys.argv[2], 'w')
    ref_path = sys.argv[3]

    amplicons = parseBedInput(bed_handle)
    createAmpliconOutput(amplicons, handle_out, ref_path)


if __name__ == "__main__":
    main()


