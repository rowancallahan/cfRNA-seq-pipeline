#!/usr/bin/env python

# Smoosh overlapping intervals, and create a new interval with id defined by <AMPL1_AMPL2>.
# BED files for amplicon coords look roughly like this.  Assuming a 0-based start coordinate, as file
# ext is BED.
# chr1	162688765	162688853	AMPL289803	0	+	.	GENE_ID=DDR2
# chr1	162688848	162688950	AMPL291704	0	+	.	GENE_ID=DDR2

# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse
import os
import sys

# Including this because is it quite useful...
# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess


VERSION = '0.2.0'

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_bed', help='Input BED file to remove overlapping genomic intervals from.')
    parser.add_argument('output_bed', help='Output BED')
    parser.add_argument('-m', '--manifest', action='store_true', help='The BED file is in manifest style, as provided by PTRL.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def main():

    args = supply_args()

    handle_bed = open(args.input_bed, 'rU')
    handle_out = open(args.output_bed, 'w')

    bed_list = []
    first = True

    with handle_bed as bed:
        
        if args.manifest:
            # Skip the header.
            next(bed)
            for line in bed:
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                start = line[1]
                stop = line[2]

                ampl_id = line[3]
                strand = line[5]
                gene = line[7]

                if first == True:
                    bed_list = [line]
                    first = False
                    i = 0
                elif first == False:

                    if int(start) < int(bed_list[i][2]) and chrom == bed_list[i][0]:
                        new_ampl_id = bed_list[i][3] + '_' + ampl_id
                        bed_list.append([chrom, start, bed_list[i][2], new_ampl_id, '0', strand, '.', gene])
                        bed_list.append([chrom, bed_list[i][2], stop, ampl_id, '0', strand, '.', gene])
                        bed_list[i][2] = start
                        i += 2
                    else:
                        bed_list.append(line)
                        i += 1
        else:

            for line in bed:
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                start = line[1]
                stop = line[2]

                if first == True:
                    bed_list = [line]
                    first = False
                    i = 0
                elif first == False:

                    if int(start) < int(bed_list[i][2]) and chrom == bed_list[i][0]:
                        bed_list.append([chrom, start, bed_list[i][2]])
                        bed_list.append([chrom, bed_list[i][2], stop])
                        bed_list[i][2] = start
                        i += 2
                    else:
                        bed_list.append(line)
                        i += 1
            

    for entry in bed_list:
        handle_out.write('\t'.join(entry))
        handle_out.write('\n')

main()
