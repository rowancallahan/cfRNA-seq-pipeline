#!/usr/bin/env python

### Create a counts file to be input for the Corless CNV R script.
### From ProbeQC:
### 1	115256501	115256590	NM_002524.4	NRAS	6183.2	93.68%	98.89%	98.89%	98.89%	98.89%
### From Amplicons BED:
### chr1	115256501	115256590	AMPL1624980	0	+	.	GENE_ID=NRAS
### Both files have headers that should be skipped.
### Output should look like:
### AmpliconId,PG2.58.013.normal,PG2.58.014.normal,PG2.59.010.normal,PG2.64.004.normal,PG2.65.005.normal,PG2.65.006.normal,PG2.69.008.normal,PG2.69.012.normal,PG2.69.015.normal,PG2.69.016.normal,RDM.15.01189.tumor
### 440774676,1465,2499,1862,1466,972,983,997,1624,836,1125,2190
### Usage: python create_counts_file.py <amplicons BED> <list of files to process> <output counts table csv>


import sys


def parseProbeQC(probeqc_file, probes):
    """
    Remove depth field from ProbeQC coverage metrics.  Match depth by genomic coordinates.
    """

    with probeqc_file as probeqc:
        next(probeqc)
        for line in probeqc:
            line = line.rstrip('\n').split('\t')
            chrom = sex_chrom_switch(line[0])
            start = line[1]
            stop = line[2]
            avgd = line[5].split('.')[0]
        
            coord = chrom + ':' + start + '-' + stop
            if coord not in probes:
                probes[coord] = [avgd]
            else:
                probes[coord].append(avgd)

    return probes


def parseAmplicons(amp_file):
    """
    Grab genomic coordinates and amplicon id.
    """

    amps = {}

    with amp_file as amplicons:
        next(amplicons)
#        amp_id = 1
        for line in amplicons:
            line = line.rstrip('\n').split('\t')
            chrom = sex_chrom_switch(line[0][3:])
            start = line[1]
            stop = line[2]

            amp_id = filter(str.isalnum, line[3])
            coord = chrom + ':' + start + '-' + stop
            amps[coord] = str(amp_id)
#            amp_id += 1
    
    return amps


def sex_chrom_switch(chrom):
    """
    Convert X and Y chroms to 23 and 24.
    """
    if chrom == "X":
        return "23"
    elif chrom == "Y":
        return "24"
    else:
        return chrom


def parse_filelist(handle_filelist):
    """
    Go through a list of files and return the file name for opening.
    """
    with handle_filelist as flist:
        for entry in flist:
            entry = entry.rstrip('\n')
            yield entry


def writeCounts(probes, amplicons, handle_out, header):
    """
    Write a counts file for the CNV Rscript.
    """

    handle_out.write(','.join(header))
    handle_out.write('\n')

    for key in amplicons:
        write_probes = ','.join(probes[key])
        handle_out.write(','.join([amplicons[key], write_probes]))
        handle_out.write('\n')

    handle_out.close


def main():

    handle_amplicons = open(sys.argv[1], 'rU')
    handle_filelist = open(sys.argv[2], 'rU')
    handle_out = open(sys.argv[3], 'w')

    probes = {}
    header = ["AmpliconId"]
    for entry in parse_filelist(handle_filelist):
        handle_probeqc = open(entry, 'rU')
        probes = parseProbeQC(handle_probeqc, probes)
        handle_probeqc.close()
        entry = filter(str.isalnum, entry)
        new_entry = entry + ".normal"
        header.append(new_entry)

    amps = parseAmplicons(handle_amplicons)

    writeCounts(probes, amps, handle_out, header)


if __name__ == "__main__":
    main()
