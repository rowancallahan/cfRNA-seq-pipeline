#!/usr/bin/env python

# Filter out entries that contain poor allele balance, as defined by args.allele_balance.
# USAGE: python filter_strand_balance_vcf.py [input_vcf] > [output_vcf]
# VERSION: 0.1.0
# John Letaw

import argparse
import vcf

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', type=file, help='Input VCF to filter based on allele balance annotations.')
#    parser.add_argument('output_vcf', help='Input VCF to filter based on allele balance annotations.')
    parser.add_argument('--allele_balance', type=float, default=.15, help='Filter sites with allele balance +/- <allele-balance>')

    args = parser.parse_args()
    return args


class AlleleBalance(vcf.filters.Base):
    """
    Filter sites by AlleleBalance
    """

    name = 'ab'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--allele_balance', type=float, default=.15,
                help='Filter sites with allele balance +/- <allele-balance>')

    def __init__(self, args):
        self.threshold = args.allele_balance

    def allele_balance(self, sample):
        ref = float(sample['AD'][0])
        alt = float(sample['AD'][1])

        try:
            ab = ref/(ref + alt)
        except:
            ab = 0.0

        return ab

    def __call__(self, record):
        for sample in record.samples:
            if sample['GT'] == '0/1':
                ab = self.allele_balance(sample)
                if ab < (0.5 + self.threshold) or ab > (0.5 - self.threshold):
                    return record.samples


def check_mend_x(variant):
    """
    Assuming we have a male proband...
    """
    chrom = variant.rstrip('\n').split('\t')[0]
    if chrom == 'X':
        father = variant.rstrip('\n').split('\t')[9]
        mother = variant.rstrip('\n').split('\t')[10]
        proband = variant.rstrip('\n').split('\t')[11]
        if (mother == '0/0' and proband == '1/1') or (mother == '1/1' and proband == '0/0'):
           return True
        else:
            return False
    else:
        return True


def main():

    args = supply_args()
    high_af = 0.5 + args.allele_balance
    low_af = 0.5 - args.allele_balance

    with args.input_vcf as vcf:
        for variant in vcf:
            if variant[0] != '#':
                if check_mend_x(variant):
                    for i in range(9,12):
                        sample = variant.rstrip('\n').split('\t')[i]
                        ad = sample.split(':')[1]
                        gt = sample.split(':')[0]
                        if gt == '0/1':
                            if len(ad.split(',')) == 2:
                                ref_ad = float(ad.split(',')[0])
                                alt_ad = float(ad.split(',')[1])
                                try:
                                    ab = ref_ad/(ref_ad + alt_ad)
                                except:
                                    ab = 0.0

                                if (ab < high_af and ab > low_af) or ab == 0.0:
                                    print(variant),
            else:
                print(variant),

if __name__ == "__main__":
    main()


