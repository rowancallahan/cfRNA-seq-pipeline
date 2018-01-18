#!/usr/bin/env python

# USAGE: remove_ref_alt_hyphens_bed.py <input_bed> <output_bed> --ref_genome <ref_genome>
# CODED BY: John Letaw

import argparse
import pysam

VERSION = '0.1.0'

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='Remove entries from a BED file that contain hyphens, as it doesn\'t match the VCF specification.')
    parser.add_argument('input_bed', type=file, help='Input BED file to remove hyphens from.')
    parser.add_argument('output_bed', help='Output BED.')
    parser.add_argument('--ref_genome', help='Reference genome to pull flanking bases from.', default="/opt/installed/galaxy_genomes/hg19/Homo_sapiens_assembly19.fasta")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def return_ref_base(ref_genome, chrom, coord):
    """
    Find the base flanking the ref and alt bases we are looking at.
    """
    to_find = chrom + ':' + coord + '-' + coord
    return pysam.faidx(ref_genome, to_find)[1].rstrip('\n')


def write_entry(handle_out, variant):
    """
    Write a variant entry to file.
    """
    handle_out.write('\t'.join(variant))
    handle_out.write('\n')


def main():

    args = supply_args()
    handle_out = open(args.output_bed, 'w')

    with args.input_bed as bed:
        for variant in bed:
            variant = variant.rstrip('\n').split('\t')
            chrom = variant[0]
            start = variant[1]
            stop = variant[2]
            ref = variant[3]
            alt = variant[4]
            if '-' in alt:
                base_to_add = return_ref_base(args.ref_genome, chrom, start)
                variant[1] = str(int(start) - 1)
                variant[2] = str(int(start) - 1)
                variant[3] = base_to_add + ref
                variant[4] = base_to_add

            elif '-' in ref:
                base_to_add = return_ref_base(args.ref_genome, chrom, stop)
                variant[3] = base_to_add
                variant[4] = base_to_add + alt

            write_entry(handle_out, variant)

    handle_out.close()

if __name__ == "__main__":
    main()


