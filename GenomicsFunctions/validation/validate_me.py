#!/usr/bin/env python

# USAGE: validate_me.py <input_vcf> <compare_bed>
# CODED BY: John Letaw

import argparse
from natsort import natsorted

VERSION = '0.1.1'

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', type=file, help='Input VCF to compare against a truth set.')
    parser.add_argument('compare_bed', type=file, help='Compare the VCF against this BED-formatted truth set.')
    parser.add_argument('outfile', help='Metrics output file.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def create_vcf_dict(vcf):
    """
    Create the VCF data structure.
    """
    vcf_dict = {}
    for variant in vcf:
        variant = variant.rstrip('\n').split('\t')
        if not variant[0].startswith('#'):
            chrom = variant[0]
            coord = variant[1]
            ref = variant[3]
            alt = variant[4]
            uniq_key = (chrom, coord, ref, alt)
            if sample not in vcf_dict:
                vcf_dict[sample] = {}
            if uniq_key not in vcf_dict[sample]:
                vcf_dict[sample][uniq_key] = variant
        else:
            if variant[0] == '#CHROM':
                sample = variant[9]

    return vcf_dict
            

def create_bed_dict(bed):
    """
    Create a dictionary to hold variant from a BED-formatted truth set.
    """
    bed_dict = {}
    for variant in bed:
        variant = variant.rstrip('\n').split('\t')
        chrom = variant[0]
        start = str(int(variant[1]) + 1)
        ref = variant[3]
        alt = variant[4]
        sample = variant[5]
        uniq_key = (chrom, start, ref, alt)
        if sample not in bed_dict:
            bed_dict[sample] = {}
        if uniq_key not in bed_dict[sample]:
            bed_dict[sample][uniq_key] = variant

    return bed_dict


def tally_matches(vcf_dict, bed_dict):
    """
    Iterate through VCF variants and determine if there is a match in the BED file.
    Return True is a match, and then do something else.
    """
    counts = {'tp': 0.0, 'fp': 0.0, 'fn': 0.0}
    sample_counts = {}
    match = {}
    no_match = {}

    for sample in vcf_dict:

        if sample not in sample_counts:
            sample_counts[sample] = counts

        for variant in vcf_dict[sample]:
            if variant in bed_dict[sample]:
                sample_counts[sample]['tp'] += 1
                if sample not in match:
                    match[sample] = []
                match[sample].append(vcf_dict[sample][variant])

            else:
                sample_counts[sample]['fp'] += 1
        
        for variant in bed_dict[sample]:
            if variant not in vcf_dict[sample]:
                sample_counts[sample]['fn'] += 1
                if sample not in no_match:
                    no_match[sample] = []
                no_match[sample].append(bed_dict[sample][variant])

    return sample_counts, match, no_match


def calc_metrics(counts):
    """
    Calculate sensitivity and PPV metrics for each sample in the truth set.
    """
    metrics = {}
    for sample in counts:
        if sample not in metrics:
            metrics[sample] = {}

        tp = counts[sample]['tp']
        fp = counts[sample]['fp']
        fn = counts[sample]['fn']
        sensitivity = tp / (tp + fn)
        ppv = tp / (tp + fp)

        metrics[sample] = [sensitivity, ppv]

    return metrics


def to_write(vcf, bed, metrics, match, no_match, handle_out):
    """
    """
    header = ['CHROM', 'COORD', 'REF', 'ALT', 'REF_COUNT', 'ALT_COUNT', 'VAF', 'VAF_EXP']
    for sample in bed:
        if sample in metrics:
            handle_out.write("SAMPLE: " + sample + '\n')
            handle_out.write("Senstivity: " + str(metrics[sample][0]) + '\n')
            handle_out.write("PPV: " + str(metrics[sample][1]) + '\n')
            handle_out.write('\n')
            handle_out.write('FOUND CALLS\n')
            handle_out.write('\t'.join(header))
            handle_out.write('\n')

        for variant in natsorted(bed[sample], key=lambda x: (x[0], x[1])): 
            if sample in match:
                if variant in vcf[sample]:
                    if vcf[sample][variant] in match[sample]:
                        chrom = vcf[sample][variant][0]
                        coord = vcf[sample][variant][1]
                        ref = vcf[sample][variant][3]
                        alt = vcf[sample][variant][4]
                        ref_count = vcf[sample][variant][9].split(':')[1].split(',')[0]
                        alt_count = vcf[sample][variant][9].split(':')[1].split(',')[1]
                        vaf = vcf[sample][variant][9].split(':')[2]
                        vaf_exp = bed[sample][variant][6]
                        to_write = [chrom, coord, ref, alt, ref_count, alt_count, vaf, vaf_exp]
                        handle_out.write('\t'.join(to_write))
                        handle_out.write('\n')

        if sample in metrics:
            handle_out.write("\nCALLS NOT FOUND\n")
        for variant in natsorted(bed[sample], key=lambda x: (x[0], x[1])):
            if sample in no_match:
                if bed[sample][variant] in no_match[sample]:
                    chrom = bed[sample][variant][0]
                    coord = bed[sample][variant][2]
                    ref = bed[sample][variant][3]
                    alt = bed[sample][variant][4]
                    vaf_exp = bed[sample][variant][6]
                    to_write = [chrom, coord, ref, alt, '', '', '', vaf_exp]
                    handle_out.write('\t'.join(to_write))
                    handle_out.write('\n')
                    
def main():

    args = supply_args()
    handle_out = open(args.outfile, 'w')
    
    # Create a dict for the input VCF.
    with args.input_vcf as vcf:
        vcf_dict = create_vcf_dict(vcf)        
    with args.compare_bed as bed:
        bed_dict = create_bed_dict(bed)

    sample_counts, match, no_match = tally_matches(vcf_dict, bed_dict)
    metrics = calc_metrics(sample_counts)

    to_write(vcf_dict, bed_dict, metrics, match, no_match, handle_out)

    handle_out.close()

if __name__ == "__main__":
    main()


