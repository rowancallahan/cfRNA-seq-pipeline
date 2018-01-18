#!/usr/bin/env python

# Currently this is set up to take output from Billy.  With a couple of field changes this
# will readily take VCF input.
# USAGE: hgvs_test.py <input> <output>
# CODED BY: John Letaw
VERSION = '0.1.1'

import argparse

import hgvs.location
import hgvs.posedit
import hgvs.edit
import hgvs.variant
import hgvs.variantmapper
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.normalizer

CHROM_MAP = {'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9', '6': 'NC_000006.11', '7': 'NC_000007.14', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10', '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9', '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10', '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}


def supply_args():
    """                                                                                                                            Populate args.                                                                                                                 https://docs.python.org/2.7/library/argparse.html                                                                              """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', type=file, help='Input VCF.')
    parser.add_argument('output_vcf', help='Output VCF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class PrepareVariants(object):
    """
    Take coordinate information from the input file and return the structure that can be used
    to pull HGVS p. and c.
    """
    def __init__(self, chrom, start, end, ref, alt):
        self.chrom = chrom
        self.start = hgvs.location.BaseOffsetPosition(base=start, datum=hgvs.location.CDS_START)
        self.end = hgvs.location.BaseOffsetPosition(base=end, datum=hgvs.location.CDS_END)
        self.iv = hgvs.location.Interval(start=self.start, end=self.end)
        # if alt != None:
        #     if 'A' in alt or 'C' in alt or 'T' in alt or 'G' in alt:
        #         self.edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
        #     else:
        #         self.edit = hgvs.edit.NARefAlt(ref=ref, alt=None)
        # else:
        self.edit = hgvs.edit.NARefAlt(ref=ref, alt=alt)
        self.posedit = hgvs.posedit.PosEdit(pos=self.iv, edit=self.edit)
        self.gvar = hgvs.variant.SequenceVariant(ac=self.chrom, type='g', posedit=self.posedit)
        

def find_start_and_stop(coord, ref):
    """
    Based on a coordinate and a reference allele string, determine
    what the start and stop coordinates should be.
    """
    start = coord
    stop = coord + len(ref) - 1
    return start, stop


def fix_refseq_no_dot(chrom, start, end, refseq, hdp):
    """
    Access UTA data source and determine proper version number of 
    RefSeq identifiers without the version number.
    """

    tx_choices = hdp.get_tx_for_region(chrom, 'splign', start, end)

    for entry in tx_choices:
        if refseq in entry[0]:
            return entry[0]

    return None


def one_or_none(in_str):
    """
    If we can slice out the everything but first character, return that string.
    Otherwise, return None.
    """
    try:
        return in_str[1:]
    except:
        return None


def fix_vcf_coords(ref, alt, coord):
    """
    Fix deletion/insertion variants in a VCF to allow for proper c. notation.
    https://groups.google.com/forum/#!topic/hgvs-discuss/psf97BkLqYw
    """
    try:
        if ref[0] == alt[0]:
            ref = one_or_none(ref)
            alt = one_or_none(alt)
            if alt == None:
                coord += 1
    except:
        if len(ref) == 0:
            ref = None
        elif len(alt) == 0:
            alt = None

    return ref, alt, coord


def write_header():
    """
    Return the header on the output file.
    """
    header = ['description', 'position_start', 'reference_base', 'variant_base', 'description', 'genotype_cdna', 'description', 'genotype_amino_acid_onel', 'genotype_amino_acid_threel', 'hgvspy_gvar', 'hgvspy_cvar', 'hgvspy_pvar', '\n']

    return '\t'.join(header)


def main():

    args = supply_args()
    hdp = hgvs.dataproviders.uta.connect()
    variantmapper = hgvs.variantmapper.VariantMapper(hdp)
    handle_out = open(args.output_vcf, 'w')

    # Write header
    header = write_header()
    handle_out.write(header)

    with args.input_vcf as myvcf:
        for line in myvcf:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')

                # Remove and replace with regular expression matching.
                if line[0] != 'description':

                    try:
                        chrom = CHROM_MAP[line[0]]
                    except:
                        chrom = CHROM_MAP[line[0][3:]]

                    coord = int(line[1])
                    ref = line[2]
                    alt = line[3]
#                    if len(ref) > 1 or len(alt) > 1:
#                        ref, alt, coord = fix_vcf_coords(ref, alt, coord) 
                    start, end = find_start_and_stop(coord, ref)
                    refseq = line[4]
                    
                    if '.' not in refseq:
                        refseq = fix_refseq_no_dot(chrom, start, end, refseq, hdp)

                    if refseq != None:
                        my_hgvs = PrepareVariants(chrom, start, end, ref, alt)
                        gvar = my_hgvs.gvar
                        print(hdp.get_tx_for_region(chrom, 'splign', start, end))
                        cvar = variantmapper.g_to_c(gvar, refseq)
                        print(str(cvar))
                        print(hgvs.normalizer.Normalizer(hdp).normalize(cvar))


                        pvar = variantmapper.c_to_p(cvar)
                        line.extend([str(gvar), str(cvar), str(pvar)])


                    handle_out.write('\t'.join(line))
                    handle_out.write('\n')

    handle_out.close()


if __name__ == "__main__":
    main()
