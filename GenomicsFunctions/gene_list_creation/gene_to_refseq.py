#!/usr/bin/env python

# INPUTS: [RefSeq(GFF3), GENE LIST TO CONVERT(TXT), ACCEPTABLE REFSEQ TRANSCRIPTS (TXT)]
# OUTPUTS: [Transcripts(TXT), INTERVALS(GATK_INTERVAL, BED)]
# USAGE: gene_to_refseq.py <input_genes> <input_transcripts> <input_gff3> <write_txs> <write_ints> <write_bed>
# CODED BY: John Letaw

import argparse

VERSION = '0.1.1'

# Maps RefSeq chromosome IDs with common.
CHROM_MAP = {'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9', '6': 'NC_000006.11', '7': 'NC_000007.13', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10', '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9', '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10', '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_genes', type=file, help='Input list of HGNC genes.')
    parser.add_argument('input_transcripts', type=file, help='Input list of accepted transcripts.')
    parser.add_argument('input_gff3', type=file, help='RefSeq GFF3 input.')
    parser.add_argument('write_txs', help='Write the transcript list to this file.')
    parser.add_argument('write_ints', help='Write the GATK interval list to this file.')
    parser.add_argument('write_bed', help='Write the BED interval list to this file.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def genes_to_find(in_genes):
    """
    Create the list of genes that we will be looking for in the RefSeq GFF3 file.
    Also will be used to take list of transcripts and place in to list.
    """
    gene_list = []
    with in_genes as genes:
        for line in genes:
            gene_list.append(line.rstrip('\n'))

    return gene_list


def parse_gff3(args, genes, txs):
    """
    Parse GFF according ot specific rules, grabbing only transcripts of interest.

    NC_000001.10	BestRefSeq	mRNA	156268415	156269428	.	-	.	ID=rna5463;Name=NM_001004319.2;Parent=gene2403;Dbxref=GeneID:391104,Genbank:NM_001004319.2,HGNC:30666,HPRD:15648;gbkey=mRNA;gene=VHLL;product=von Hippel-Lindau tumor suppressor-like;transcript_id=NM_001004319.2

    We just want the entries that match the gene name (gene=xxx) and where the RefSeq transcript id matches the list of
    approved transcripts (Name=xxx).  Only look at "mRNA" entries.

    """
    final_gff_list = []
    write_txs = open(args.write_txs, 'w')
    write_ints = open(args.write_ints, 'w')
    write_bed = open(args.write_bed, 'w')
    with args.input_gff3 as gff:
        for line in gff:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                entry_type = line[2]
                if entry_type == 'mRNA':
                    chrom = line[0]
                    start = line[3]
                    stop = line[4]
                    gent = find_term(line[8], 'gene')
                    txt = find_term(line[8], 'Name')
                    hgnc = line[8].split(';')[gent].split('=')[1]
                    refseq = line[8].split(';')[txt].split('=')[1]
                    if hgnc in genes and refseq in txs and chrom.startswith('NC_'):
                        # Write transcript IDs to file.
                        print("Writing list of transcripts")
                        write_txs.write(refseq)
                        write_txs.write('\n')
                        # Write intervals in GATK format to file.
                        print("Writing intervals in gatk_interval format")
                        chrom = find_value(chrom)
                        interval = chrom + ':' + start + '-' + stop
                        write_ints.write(interval)
                        write_ints.write('\n')
                        # Write BED intervals to file.'
                        print("Writing intervals in BED format")
                        start = int(start) - 1
                        to_write = [chrom, str(start), stop, '\n']
                        write_bed.write('\t'.join(to_write))

                        final_gff_list.append(line)

    write_bed.close()
    write_txs.close()
    write_ints.close()
    return final_gff_list


def find_value(chrom):
    """
    Dig up the key value corresponding to the RefSeq formatted chrom id.
    """
    for key, value in CHROM_MAP.iteritems():
        if value == chrom:
            return key
 

def find_term(blob, term):
    """
    Find the index of a particular term in a GFF INFO section.
    """
    blob = blob.split(';')
    for entry in blob:
        if entry.split('=')[0] == term:
            return blob.index(entry)


def main():

    args = supply_args()

    gene_list = genes_to_find(args.input_genes)
    tx_list = genes_to_find(args.input_transcripts)
    gff_list = parse_gff3(args, gene_list, tx_list)

if __name__ == "__main__":
    main()


