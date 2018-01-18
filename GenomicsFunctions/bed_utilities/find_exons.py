#!/usr/bin/env python

# USAGE: find_exons.py <GFF3 File> <RefSeq transcript list> <BED file of intervals>
# CODED BY: John Letaw
# Jira: BCORE-224

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

refseq_chrom = {'NC_000001.10':'1', 'NC_000002.11':'2', 'NC_000003.11':'3', 'NC_000004.11':'4', 'NC_000005.9':'5', 'NC_000006.11':'6', 'NC_000007.13':'7', 'NC_000008.10':'8', 'NC_000009.11':'9', 'NC_000010.10':'10', 'NC_000011.9':'11', 'NC_000012.11':'12', 'NC_000013.10':'13', 'NC_000014.8':'14', 'NC_000015.9':'15', 'NC_000016.9':'16', 'NC_000017.10':'17', 'NC_000018.9':'18', 'NC_000019.9':'19', 'NC_000020.10':'20', 'NC_000021.8':'21', 'NC_000022.10':'22', 'NC_000023.10':'X', 'NC_000024.9':'Y', 'NC_012920.1':'MT'}


def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_gff', type=file, help='')
    parser.add_argument('input_refseq', type=file, help='')
    parser.add_argument('input_bed', type=file, help='')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def find_entry(info, field):
    """
    Find an entry from the last column of the GFF3 file.
    inputs: info string, field you want
    """
    for entry in info.split(';'):
        if entry.split('=')[0] == field:
            return entry.split('=')[1]


def make_list(fileh):
    """
    Create a list out of a one column input file.
    """
    my_list = []
    with fileh as my_file:
        for line in my_file:
            my_list.append(line.rstrip('\n'))

    return my_list


def grab_intervals(fileh):
    """
    chr1	2488077	2488302	TNFRSF14	0	+
    """
    to_check = []
    with fileh as my_file:
        for line in my_file:
            if line.startswith('chr'):
                line = line.rstrip('\n').split('\t')
                chrom = line[0][3:]
                start = str(int(line[1]) + 1)
                stop = line[2]
                to_check.append((chrom, start, stop))
    return to_check


def grab_exons(handle_gff, my_refseq):
    """
    Find the exon entries from the GFF3 file based on transcript list.
    """
    coords = {}
    with handle_gff as gff:
        for line in gff:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                entry_type = line[2]
                info = line[8]
                if entry_type == 'exon' and line[0] in refseq_chrom:
                    refseq = find_entry(info, 'transcript_id')
                    chrom = refseq_chrom[line[0]]
                    start = line[3]
                    stop = line[4]
                    if refseq in my_refseq:
                        if refseq not in coords:
                            coords[refseq] = [(chrom, start, stop)]
                        else:
                            coords[refseq].append((chrom, start, stop))
    return coords


def get_match(bed_coord, my_coords):
    """
    Find match between the interal input list and the exon list from the GFF3.
    """
    for coord in my_coords:
        for trio in my_coords[coord]:
            if bed_coord[0] == trio[0]:
                if not ((bed_coord[1] < trio[1] and bed_coord[2] < trio[1]) or (bed_coord[1] > trio[2] and bed_coord[2] > trio[2])): 
                    return '\t'.join([bed_coord[0], str(int(bed_coord[1])-1), bed_coord[2], coord, str(my_coords[coord].index(trio)+1)])
                
    return '\t'.join([bed_coord[0], str(int(bed_coord[1])-1), bed_coord[2]])


def main():

    args = supply_args()
    to_check = grab_intervals(args.input_bed)
    my_refseq = make_list(args.input_refseq)
    my_coords = grab_exons(args.input_gff, my_refseq)

    # Send to stdout and redirect.
    for bed_coord in to_check:
        print(get_match(bed_coord, my_coords))


if __name__ == "__main__":
    main()


