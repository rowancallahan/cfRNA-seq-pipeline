#!/usr/bin/env python

# Script to perform a few custom filtering steps on Pindel output.
# USAGE: python pindel_filtering.py <input_vcf> <output_vcf>
# CODED BY: John Letaw

import argparse

VERSION = '0.1.0'

def compareAF(af1, af2, margin, thresh):
    """
    Make a comparison between the allele frequencies, and remove variants that don't
    differ by at least margin and aren't greater than thresh.
    The idea is to keep those with af close to 0, so we don't lose low af somatic variants.
    """

    af1 = float(af1)
    af2 = float(af2)

    if af1 >= af2:
        diff = af1-af2
    else:
        diff = af2-af1    

    if diff <= margin:
        if af1 <= thresh or af2 <= thresh or af1 >= 1-thresh or af2 >= 1-thresh:
            return False
    elif diff > margin:
        return False

    return True


def checkGenos(geno1, geno2):
    """
    If the genotypes are the same, do not write to output.
    """
    if geno1 == geno2:
        return True

    return False


def calcAF(count, total):
    """
    Calculate the allele frequency and format the output.
    """
    
    count = int(count) + 0.0
    total = int(total) + 0.0

    if total != 0:
        af = float(count/total)
        if af != 0.0 and af != 1.0:
            af = "{:.7f}".format(count/total).rstrip('0')
    else:
        af = 0.0

    return str(af)


def calcDP(counts):
    """
    Pindel VCF's do not contain a DP value in the SAMPLE field.  Calculate it here.
    'counts' should be a list of values.
    """

    total = 0
    for count in counts:
        total += int(count)
    return str(total)


def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='Perform custom filtering on Pindel output.')
    parser.add_argument('infile', help='Input VCF')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('--margin', type=float, default=0.2, help='Allele frequencies with spread lower than this value will be considered the same if they are also not below the threshold defined below.')
    parser.add_argument('--thresh', type=float, default=0.1, help='If allele frequencies are below this value, this record will be included regardless.')
    parser.add_argument('--normal_af', type=float, default=0.05, help='VAF cutoff for normal sample.')
    parser.add_argument('--tumor_af', type=float, default=0.01, help='VAF cutoff for tumor sample.')
    parser.add_argument('--depth', type=int, default=5, help='Total depth cutoff for filtering.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def main():

    # 1	884041	.	C	CCCTGGCTGCACCCTGGTCCCCCTGGTCCCTTTGGCCCTGCA	.	PASS	END=884041;HOMLEN=59;HOMSEQ=CCTGGCTGCACCCTGGTCCCCCTGGTCCCTTTGGCCCTGCACCTGGCTGCACCCTGGTC;SVLEN=41;SVTYPE=INS	GT:AD	0/0:92,4	0/0:41,1
    args = supply_args()
    handle = open(args.infile, 'rU')
    handle_out = open(args.outfile, 'w')

    with handle as myvcf:
        for variant in myvcf:
            if variant[0] != "#":

                split_var = variant.rstrip('\n').split('\t')
                normal = split_var[9]
                tumor = split_var[10]
                normal_geno = normal.split(':')[0]                                                                                                                    
                tumor_geno = tumor.split(':')[0]                                                                                                                      
                normal_ref = int(normal.split(':')[1].split(',')[0])
                normal_alt = int(normal.split(':')[1].split(',')[1])
                normal_total = normal_ref + normal_alt
                normal_af = calcAF(normal_alt, normal_total)
                tumor_ref = int(tumor.split(':')[1].split(',')[0])
                tumor_alt = int(tumor.split(':')[1].split(',')[1])
                tumor_total = tumor_ref + tumor_alt
                tumor_af = calcAF(tumor_alt, tumor_total)

                # Check if genotypes match.
                if checkGenos(normal_geno, tumor_geno) == False:
                    # Check if allele frequencies are different by at least an amount defined by 
                    # margin and thresh.
                    if compareAF(float(normal_af), float(tumor_af), args.margin, args.thresh) == False:
                        # Check to make sure the depth is at least args.depth.
                        if int(calcDP([tumor_ref, tumor_alt])) > args.depth:
                            # Check if the normal VAF is less than a certain amount.
                            if float(normal_af) < args.normal_af:
                                # Check if the tumor VAF is greater than a certain amount.
                                if float(tumor_af) >= args.tumor_af:
                                    if (float(tumor_af) - float(normal_af)) > args.thresh:
                                        handle_out.write(variant) 
                                    else:
                                        print("diff kicked out: \n" + variant)
                                else:
                                    print("tumor_af kicked out: \n" + variant)
                            else:
                                print("normal_af kicked out: \n" + variant)
                        else:
                            print("calcDP kicked out: \n" + variant)
                    else:
                        print("compareAF kicked out: \n" + variant)
                else:
                    print("checkGenos(normal_geno, tumor_geno) kicked out: \n" + variant) 

            else:
                handle_out.write(variant)

    handle_out.close()

if __name__ == "__main__":
    main()


