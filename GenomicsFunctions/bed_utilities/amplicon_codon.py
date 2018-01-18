#!/usr/bin/env python

# Write a file that contains codon coordinates for all reginos in a BED
# coordinate file.
# chromosome, start position, end position, gene, transcript, first exon \in
#  the region, start codon for the region, end codon for the region
# USAGE: python amplicon_codon.py <Interval BED file> <RefSeq GFF3> <RefSeq Transcript List>
# CODED BY: John Letaw

VERSION = '0.2.0'

import sys
from gff3 import GffReader
#from gff3 import GffLine
from bed import BedReader


def coding_to_codon(coords):

    """
    Convert the coding coordinates to actual codon coordinates.  If the
    number is not divisible by 3, we will treat any remainder as a whole
    number coordinate value.
    :param coords:
    :return codon_coords:
    """

    codon_coords = []

    for coord in coords:
        start_coord = coord[0]
        end_coord = coord[1]

        if start_coord == 1:
            new_start = 1
        elif start_coord % 3 != 0:
            new_start = start_coord / 3 + 1
        else:
            new_start = start_coord / 3

        if end_coord % 3 != 0:
            new_end = end_coord / 3 + 1
        else:
            new_end = end_coord / 3

        codon_coords.append([new_start, new_end])

    return codon_coords


def place_coord(bed, gff):

    """
    Connect coordinates from the BED file to the map structure.  This will
    help ensure we only develop codon coordinates for the portion of BED
    coordinate that overlaps with GFF coding features.
    :param bed:
    :param gff:
    :return bed_map, len(temp_map):
    """

    bed_map = []
    temp_map = map_coords(gff)
    for bcoord in bed:
        start = bcoord[0]
        stop = bcoord[1]

        if start in temp_map:
            new_start = temp_map[start]
        else:
            while start not in temp_map:
                if min(temp_map, key=temp_map.get) < start:
                    start = start - 1
                else:
                    start = start + 1

            new_start = temp_map[start]

        if stop in temp_map:
            new_stop = temp_map[stop]
        else:
            while stop not in temp_map:
                if min(temp_map, key=temp_map.get) > stop:
                    stop = stop + 1
                else:
                    stop = stop - 1

            new_stop = temp_map[stop]

        bed_map.append([new_start, new_stop])
    return bed_map, len(temp_map)


def map_coords(coords):

    """
    Map the coordinates to integer values, making it easier to convert from
    genomic coordinates to coding coordinates.
    :param coords:
    :return temp_map:
    """

    temp_map = {}
    j = 1
    for coord in coords:
        for i in range(int(coord[0]), int(coord[1]) + 1):
            temp_map[i] = j
            j += 1
    return temp_map


def locate_coords(parent_coords, bcoord, chrom):

    """
    Find locations for the BED coordinates, so that we can connect them to a
    parent id from the GFF.

    :param parent_coords:
    :param bcoord:
    :param chrom:
    :return parent or None:
    """

    if chrom in parent_coords:
        for parent in parent_coords[chrom]:
            if (int(bcoord[0]) > int(parent_coords[chrom][parent][0][1])) or \
                    (int(bcoord[1]) < int(parent_coords[chrom][parent][0][0])):
                pass
            else:
                return parent

    return None


def flip_coords(codon_coords, length):

    """
    Flip codon coordinates, for the case where we are working on the
    opposite strand.
    :param codon_coords:
    :param length:
    :return new_coords:
    """

    new_coords = [[length/3 - coord + 1 for coord in coords] for coords in \
            codon_coords]
    new_coords.reverse()

    for coords in new_coords:
        coords.reverse()

    return new_coords


def write_codon_coords(outfile, chrom, bed_parent, parent,
                       refseq_parent, final_coords, exon_number,
                       hgnc_parent):

    """
    Write the codon coordinates to an output file.
    chromosome, start position, end position, gene, transcript, first exon
    in the region, start codon for the region, end codon for the region

    :param outfile:
    :param chrom:
    :param bed_parent:
    :param parent:
    :param refseq_parent:
    :param final_coords:
    :param exon_number:
    :param hgnc_parent:
    :return output_file:
    """

    bed_index = 0

    for coord in final_coords:

        #bed_index = final_coords.index(coord)
        to_write = [str(chrom), str(bed_parent[parent][bed_index][0]),
                    str(bed_parent[parent][bed_index][1]),
                    hgnc_parent[chrom][parent],
                    refseq_parent[chrom][parent],
                    str(exon_number[bed_index]),
                    str(coord[0]),
                    str(coord[1]), '\n']

        outfile.write('\t'.join(to_write))
        bed_index += 1


def find_exon_number(bed_parent, exon_parent, parent, chrom):

    """
    Given a CDS sequence, find which exon it exists within.  This way,
    we can write the exon number to the output file.
    :return temp_exons:
    """

    temp_exons = []
    curr_index = 0
    if parent not in bed_parent:
        print("Not in bed_parent " + parent)
    else:
        for coord in bed_parent[parent]:
            start = coord[0]
            stop = coord[1]
            for ecoord in exon_parent[chrom][parent]:
                if start > int(ecoord[1]) or stop < int(ecoord[0]):
                    if exon_parent[chrom][parent].index(ecoord) == len(
                            exon_parent[chrom][parent])-1:
                        temp_exons.append(curr_index)
                else:
                    curr_index = exon_parent[chrom][parent].index(ecoord)+1
                    temp_exons.append(exon_parent[chrom][parent].index(ecoord)+1)
                    break

    return temp_exons


def decide_flip(my_gff, chrom, parent):

    """
    If the CDS coordinates are listed in reverse order, flip them.
    :return True or False:
    """

    if my_gff.cds_parent[chrom][parent][0][0] > \
            my_gff.cds_parent[chrom][parent][1][0]:
        return True
    else:
        return False


def create_bed_parent(my_gff, my_bed):

    """
    Make a dictinoary of GFF parent id's to BED regions.
    :return bed_parent:
    """

    bed_parent = {}
    parent_coords = my_gff.coords_parent
    for chrom in my_bed.bed_ints:
        for bcoord in my_bed.bed_ints[chrom]:
            bparent = locate_coords(parent_coords, bcoord, chrom)
            if bparent != None:
                if bparent not in bed_parent:
                    bed_parent[bparent] = [bcoord]
                else:
                    bed_parent[bparent].append(bcoord)

    return bed_parent


def main():

    my_bed = BedReader(sys.argv[1])
    my_gff = GffReader(sys.argv[2], sys.argv[3])
    outfile = open(sys.argv[4], 'w')

    # Write the header to file.
    header = ["Chromosome", "Start", "Stop", "HGNC",
              "RefSeq", "Exon Number", "Start Codon",
              "Stop Codon", '\n']
    outfile.write('\t'.join(header))

    bed_parent = create_bed_parent(my_gff, my_bed)

    i = 0
    for chrom in my_gff.cds_parent:
        for parent in my_gff.cds_parent[chrom]:

            if i % 10 == 0:
                print(i)
            i += 1

            if len(my_gff.cds_parent[chrom][parent]) > 1:
                flip = decide_flip(my_gff, chrom, parent)

            if parent in bed_parent:
                placing, length = place_coord(bed_parent[parent], sorted(
                    my_gff.cds_parent[chrom][parent]))
                final_coords = coding_to_codon(placing)
                if flip == True:
                    final_coords = flip_coords(final_coords, length)

            curr_exons = find_exon_number(bed_parent, my_gff.exon_parent,
                                        parent, chrom)
            if flip:
                curr_exons = curr_exons[::-1]
                bed_parent[parent] = bed_parent[parent][::-1]

            write_codon_coords(outfile, chrom, bed_parent, parent,
                               my_gff.refseq_parent, final_coords,
                               curr_exons, my_gff.hgnc_parent)


if __name__ == "__main__":
    main()
