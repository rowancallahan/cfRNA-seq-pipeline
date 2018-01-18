#!/usr/bin/env python

### Add an allele frequency and total depth entry, useful for Pindel VCF's.
### USAGE: python add_sample_format_info.py <input VCF> <output VCF>
# TODO: Split filtering functionality out of this script in to with another tool
# or a separate script.  This script was originally put together for one specific
# use case (Pindel) but now may be more widely applicable.

import sys
import argparse

VERSION = '0.5.0'

def createHeaderEntry(attrib):
    """
    This goes in the header of the VCF, and looks like this:
    ##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth, how many reads support this allele">
    Modify to allow for different Number, Type, and Description field to be passed.
    """

    if attrib == "AF":
        header = "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Variant allele frequency\">\n"
    elif attrib == "DP":
        header = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at this locus.\">\n"
        
    return header


def findIndex(splitting, delim, to_find):
    """
    Find the index of the field in the FORMAT column you care about.  Usually will be AD.
    """

    for entry in splitting.split(delim):
        if entry == to_find:
            curr_index = splitting.split(delim).index(to_find)
            return curr_index

    return None


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


def check_format(format, field, delim):
    """
    Check to see if the field you want to add is already there.
    """
    if field in format.split(delim):
        return True
    return False


def gather_header_terms(header_line):
    """
    Gather everything from the header, and determine which fields are already there.
    Can look like:
    ##INFO=<ID=NTLEN,Number=.,Type=Integer,Description="Number of bases inserted in place of deleted code">
    ##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
    """

    return header_line.split('=')[2].split(',')[0]


def replace_field(info):
    """
    Replace one of the fields in the INFO column of a VCF.
    For instance:
    0/1:1253,5:2.037e-03:5:0:1.00:41072,171:659:594:1258
    Want to replace AF in scientific notation so we can compare it.
    """
    pass


def create_format(format, delim, field):
    """
    Take the current FORMAT string and return a new one based on which fields
    we would like to add to the VCF.
    """
    return delim.join([format, field])


def supply_args():
    """
    Populate args.                                         
    https://docs.python.org/2.7/library/argparse.html
    """

    parser = argparse.ArgumentParser(description='Add information to VCF FORMAT and SAMPLE columns.')
    parser.add_argument('infile', help='Input VCF')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('format_label', type=str, choices=['AD', 'DPR'], help='Label to search for in the format field, such as AD or DPR.')
    parser.add_argument('--add_af', action='store_true', help='Should the AF field be added to the FORMAT and SAMPLE columns?')
    parser.add_argument('--add_dp', action='store_true', help='Should the DP field be added to the FORMAT and SAMPLE columns?')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    return args


def main():

    args = supply_args()

    handle_vcf = open(args.infile, 'rU')
    handle_out = open(args.outfile, 'w')
    delim = ':'  # We will usually target ':' delimiters, as that is what the VCF uses.
    to_find = args.format_label  # AD or DPR
    header_terms = []

    if not args.add_af and not args.add_dp:
        raise Exception("You have chosen to not write an AF or DP field, therefore this script has no purpose.")

    with handle_vcf as myvcf:
        for variant in myvcf:
            if variant[0] != "#":

                curr_index = []
                split_variant = variant.rstrip('\n').split('\t')
                format = split_variant[8]
                to_write = split_variant[:8]

                if curr_index == []:
                    if to_find == 'DPR':
                        curr_index.append(findIndex(format, delim, 'RO'))
                        curr_index.append(findIndex(format, delim, 'AO'))
                    elif to_find == 'AD':
                        curr_index.append(findIndex(format, delim, to_find))
                if None in set(curr_index):
                    continue

                if args.add_af:
                    field = "AF"
                    if check_format(format, field, delim) == False:
                        format = create_format(format, delim, field)
                        add_af = True
                    else:
                        add_af = False

                if args.add_dp:
                    field = "DP"
                    if check_format(format, field, delim) == False:
                        format = create_format(format, delim, field)
                        add_dp = True
                    else:
                        add_dp = False

                print(format)

                # Append FORMAT string to to_write.
                to_write.append(format)

                for samp in range(9, len(split_variant)):

                    sample = split_variant[samp]
                    
                    if to_find == 'AD':
                        ref = str(sample.split(delim)[curr_index[0]].split(',')[0])
                        alt = sample.split(delim)[curr_index[0]].split(',')[1:]
                    elif to_find == 'DPR':
                        ref = str(sample.split(delim)[curr_index[0]])
                        alt = sample.split(delim)[curr_index[1]].split(',')

                    all_counts = alt + [ref]
                    total = calcDP(all_counts)

                    if args.add_af:
                        af = []
                        for count in alt:
                            af.append(calcAF(count, total))
                        if add_af == True:
                            sample = delim.join([sample, str(','.join(af))])
                    if args.add_dp:
                        if add_dp == True:
                            sample = delim.join([sample, total])

                    to_write.append(sample)

                handle_out.write('\t'.join(to_write))
                handle_out.write('\n')

            else:
                
                if "CHROM" not in variant:
                    handle_out.write(variant)
                    if "INFO" in variant or "FORMAT" in variant:
                        header_terms.append(gather_header_terms(variant))
                else:
                    if "AF" not in header_terms:
                        handle_out.write(createHeaderEntry("AF"))
                    if "DP" not in header_terms:
                            handle_out.write(createHeaderEntry("DP"))
                    handle_out.write(variant)

    handle_out.close()


if __name__ == "__main__":
    main()
