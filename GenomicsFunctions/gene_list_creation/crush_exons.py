#!/usr/bin/python

### Input a TSV from BioMart that looks like:
### <chrom><start><stop><ense><ensg><hgnc>
### If exon coordinates are identical, pick an identifier.
### Using coding coordinates only right now, as that is what KDL cares about.
### if an hgnc/ense pair can correspond to different coordinates, look for larger coordinate range.

# In BioMart (grch37.ensembl.org/biomart)
# Select hsa, under filter:
# chroms 1-22,X,Y
# with HGNC symbols only
# protein_coding
# under attributes/structures:
# chromsome name
# genomic coding start
# genomic coding end
# ensembl exon id
# ensembl gene id
# associated gene name

import sys

def createKDL(handle_kdl):

    kdl_list = []

    with handle_kdl as kdl:
        for line in kdl:
            kdl_list.append(line.rstrip('\n'))

    return kdl_list

def checkDupExon(exons):

    """
    exons[chrom:start-stop] = [ense, [hgnc]]
    """

    genes = {}

    for exon in exons:

        hgnc = exons[exon][1]

        for i in range(len(hgnc)):

            chrom = exon.split(':')[0]
            start = exon.split(':')[1].split('-')[0]
            stop = exon.split(':')[1].split('-')[1]
            ense = exons[exon][0]
            key = hgnc[i] + "_" + ense

            if key not in genes:
                genes[key] = [chrom, start, stop, hgnc[i], ense]
            else:
                if chrom == genes[key][0] and start == genes[key][1] and stop >= genes[key][2]:
                    genes[key] = [chrom, start, stop, hgnc[i], ense]

    return genes


def writeExons(dups, handle_out):

    """
    dups[hgnc_ense] = [chrom, start, stop, hgnc, ense]
    """

    handle_out.write("CHROM\tSTART\tSTOP\tHGNC\tENSE\n")

    for dup in sorted(dups.values(), key=lambda x: (x[0], int(x[1]))):
#        handle_out.write('\t'.join([dups[dup][0], dups[dup][1], dups[dup][2], dups[dup][3], dups[dup][4]]))
        handle_out.write('\t'.join([dup[0], dup[1], dup[2], dup[3], dup[4]]))
        handle_out.write('\n')

    handle_out.close()


def createExons(handle):

    exon_dict = {}

    with handle as exons:
        next(exons)
        for line in exons:
            line = line.rstrip('\n').split('\t')

            chrom = line[0]
            start = line[1]
            stop = line[2]
            ense = line[3]
            ensg = line[4]
            hgnc = line[5]

            if start == "" or stop == "" or chrom == "":
                continue

            coord = chrom + ':' + start + '-' + stop
            if coord not in exon_dict:
                exon_dict[coord] = [ense, [hgnc]]
            else:
                if hgnc not in exon_dict[coord][1]:
                    exon_dict[coord][1].append(hgnc)
                    if int(ense[4:]) > int(exon_dict[coord][0][4:]):
                        exon_dict[coord][0] = ense

    return exon_dict


def main():

    handle = open(sys.argv[1], 'rU')
    handle_kdl = open(sys.argv[2], 'rU')
    handle_out = open(sys.argv[3], 'w')

    exon_dict = createExons(handle)
#    kdl_list = createKDL(handle_kdl)
    remove_dups = checkDupExon(exon_dict)
    writeExons(remove_dups, handle_out)

   
if __name__ == "__main__":
    main()
