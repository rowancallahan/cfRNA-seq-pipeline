#!/usr/bin/env python
# Merge intervals in a BED file that is formatted as such:
# X	100627109	100629986	ENSG00000000003
# USAGE: python new_merge.py <input_bed> <output_bed>

import sys

def merge_intervals(intervals):
    """
    Input: intervals, a list of tuples [(x,y), (x,y), (x,y), ...]
    Output: merged list of tuples
    """

    # Sort the incoming interval list according to start position.
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []

    for interval in intervals:
        if not merged:
            merged.append(interval)
        else:
            last_coord = merged.pop()
            # Check for overlap, if not, append both the previous coordinate and new_coordinate to merged.
            if last_coord[1] >= interval[0]:#we know interval[0]>last_coord[0]
                new_ival = (min(last_coord[0],interval[0]), max(last_coord[1],interval[1]) )
                merged.append(new_ival)
            else:
                merged.append(last_coord)
                merged.append(interval)

    return merged


def write_bed(gene_coords, outfile):
    """
    Input is a dictionary that looks like:
    {"ensg_chrom": (start, stop)} 
    Close the handle after processing.
    """

    for entry in sorted(gene_coords.items(), key=lambda(k,v): (k.split('_')[1], int(v[0][0]))):
        
        key = entry[0]
        item = entry[1]
        gene = key.split('_')[0]
        chrom = key.split('_')[1]

        for coord in item:
            start = coord[0]
            stop = coord[1]
            outfile.write('\t'.join([chrom, start, stop, gene]))
            outfile.write('\n')

    outfile.close()


def main():
    input_bed = open(sys.argv[1], 'rU')
    output_bed = open(sys.argv[2], 'w')

    all_ints = {}

    with input_bed as bed:
        for coord in bed:

            try:
                this_coord = coord.rstrip('\n').split('\t')[:4]
            except:
                raise Expection("Your BED column probably does not have a column of gene identifiers.")

            chrom = this_coord[0]
            start = this_coord[1]
            stop = this_coord[2]
            ensg = this_coord[3]

            # Create a dictionary with unique gene and chrom keys.
            uniq_key = ensg + '_' + chrom
            if uniq_key not in all_ints:
                all_ints[uniq_key] = [(start, stop)]
            else:
                all_ints[uniq_key].append((start, stop))

    # Perform merging operation.
    for key in all_ints:
        all_ints[key] = merge_intervals(all_ints[key])

    # Write the BED file.
    write_bed(all_ints, output_bed)

if __name__ == "__main__":
    main()
