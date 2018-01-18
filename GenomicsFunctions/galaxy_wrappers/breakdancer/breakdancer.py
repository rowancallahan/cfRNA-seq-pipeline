#!/usr/bin/env python

# Wrapper for BreakDancer.
# All options are currently supported.
# CODED BY: John H. Letaw

from __future__ import print_function
from argparse import ArgumentParser
import logging
import sys
import pysam
import os

# Including this because is it quite useful...
# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess


VERSION = '0.2.2'
PROGRAM_VERSION = '1.4.5'

def bd_argparse():
    parser = ArgumentParser(description='Run BreakDancer')
    ### Required
    parser.add_argument("--config_file", action="store", help="BreakDancer Analysis Config File", default="bd.config")

    ### Optional
    parser.add_argument("-o", help="operate on a single chromosome [all chromosomes]")
    parser.add_argument("-s", help="minimum length of a region")
    parser.add_argument("-c", help="cutoff in unit of standard deviation")
    parser.add_argument("-m", help="maximum SV size")
    parser.add_argument("-q", help="minimum alternative mapping quality")
    parser.add_argument("-r", help="minimum number of read pairs required to establish a connection")
    parser.add_argument("-x", help="maximum threshold of haploid sequence coverage for regions to be ignored")
    parser.add_argument("-b", help="buffer size for building connection")
    parser.add_argument("-t", help="only detect transchromosomal rearrangement", action="store_true")
    parser.add_argument("-d", help="prefix of fastq files that SV supporting reads will be saved by library")
    parser.add_argument("-g", help="dump SVs and supporting reads in BED format for GBrowse")
    parser.add_argument("-l", help="analyze Illumina long insert (mate-pair) library", action="store_true")
    parser.add_argument("-a", help="print out copy number and support reads per library rather than per bam", action="store_true")
    parser.add_argument("-f", help="print out Allele Frequency column", action="store_true")
    parser.add_argument("-y", help="output score filter")

    ### Wrapper arguments
    parser.add_argument("--input_tumor", help="Input Tumor BAM", required=True)
    parser.add_argument("--input_normal", help="Input Normal BAM", required=True)
    parser.add_argument("--output_file", help="Output Breakpoints", required=True)
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    # Process arguments
    args = parser.parse_args()
    
    return args


def symlink_and_index_bams(args, tfile, nfile):
    """
    Index the BAM file.
    """
    os.symlink(args.input_tumor, tfile)
    os.symlink(args.input_normal, nfile)
    pysam.index(tfile)
    pysam.index(nfile)


def head_bam(infile, lines=100000):
    """
    Take first x lines from a BAM file then write to a new file.
    TODO: Get this to work as a stream instead.
    """
    samfile = pysam.AlignmentFile(infile, "rb")
    sam_out = pysam.AlignmentFile('head_bam_temp.bam', 'wb', template=samfile)
    
    i = 0
    for entry in samfile:
        if i <= lines:
            sam_out.write(entry)
            i += 1
        else:
            return 0

    sam_out.close()


def find_stats():
    """
   Find the necessary stats to plug in to the BreakDancer config file.

    SN	insert size average:229.7
    SN	insert size standard deviation:55.1
    """

    bam_stats = []
    stats = pysam.stats('head_bam_temp.bam')
    for line in stats.split('\n'):
        if line.startswith('SN') and 'insert size' in line or 'maximum length' in line:
            line = line.split('\t')
            bam_stats.append(line[2])
    
    return bam_stats


def build_cmd(options, config):
    """
    Build this command in to a list for execution.
    """
    ignore = ('input_normal', 'input_tumor', 'config_file', 'output_file')
    cmd = ['breakdancer-max', options['config_file']]
    for arg, value in options.items():
        if value:
            if arg not in ignore:
                cmd.append('-{0}'.format(arg))
                cmd.append(value)

    return cmd


def run_cmd(cmd):
    """
    Run command.
    """
    logging.info("RUNNING: %s" % (cmd))
    print('Running the following command:')
    print('\t'.join(cmd))

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    print("From Output: " + stdout)
    if stderr:
        eprint("From Error: " + stderr)

    return (p.wait(), stdout, stderr)


def eprint(*args, **kwargs):
    """
    http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


def writeTemp(outfile, normal_stats, tumor_stats, input_normal, input_tumor):

    outfile.write("map:" + input_tumor + "\tmean:" + tumor_stats[1] + "\tstd:" + tumor_stats[2] + "\treadlen:" + tumor_stats[0] + "\tsample:tumor\texe:samtools view\n")
    outfile.write("map:" + input_normal + "\tmean:" + normal_stats[1] + "\tstd:" + normal_stats[2] + "\treadlen:" + normal_stats[0] + "\tsample:normal\texe:samtools view\n")
    outfile.close()
    

def main():
    args = bd_argparse()
    tumor_file = "tumor.bam"
    normal_file = "normal.bam"
    symlink_and_index_bams(args, tumor_file, normal_file)

    print("Calculating normal stats.")
    head_bam(normal_file)
    normal_stats = find_stats()

    print("Calculating tumor stats.")
    head_bam(normal_file)
    tumor_stats = find_stats()

    print("Writing config file.")
    writeTemp(open(args.config_file, 'w'), normal_stats, tumor_stats, normal_file, tumor_file)
    this_cmd = build_cmd(vars(args), args.config_file)
    to_write = run_cmd(this_cmd)[1]
    handle_out = open(args.output_file, 'w')
    handle_out.write(to_write)
    handle_out.close()

if __name__ == "__main__":
    main()
