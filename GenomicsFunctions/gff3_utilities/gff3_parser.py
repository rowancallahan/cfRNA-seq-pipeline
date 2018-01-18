# VERSION: 0.1.0

"""
Contains:
GFFReader
GFFLine
"""

class GffReader(object):
    """
    Reads a GFF3 file and creates several data structures relating portions
    of the GFF info field to the parent id of the particular feature.

    Currently this takes a list of RefSeq transcript id's as well, to limit
    the search space for BED coordinates.  This will change when I have a
    chance to put more effort in to this.

    Input: filename, refseq_file

    * Describe seqnames
    region

    * No parent
    gene

    *Parent is "gene"
    primary_transcript
    transcript
    mRNA
    ncRNA
    rRNA
    tRNA

    *Parent is in above group.
    CDS
    exon

    *Singletons
    cDNA_match
    match
    C_gene_segment
    D_gene_segment
    J_gene_segment
    D_loop
    V_gene_segment

    """
    def __init__(self, filename, refseq_file):

        """
        Set following data structure:

        refseq - List of RefSeq id's, from file
        header - List of lines containing # as line.startswith()
        refseq_parent - Relates RefSeq NM values to the feature parent id
        hgnc_parent - Relates HGNC symbols to the feature parent id
        cds_parent - Relates CDS coordinates to the feature parent id
        exon_parent - Relates exon coordinates to the feature parent id
        coords_parent - Relates mRNA transcript coordinates to the feature
        parent id

        :param filename:
        :param refseq_file:
        """

        self.filename = filename
        self.refseq = self._parse_refseq(refseq_file)
        self.header = self._get_header()
        self.refseq_parent = self._refseq_parent(self.refseq, feature='mRNA',
                                                 ident='Name')
        self.hgnc_parent = self._info_parent(self.refseq_parent,
                                               feature='mRNA',
                                               ident='gene')
        self.cds_parent = self._cds_parent(self.refseq_parent, feature='CDS')
        self.exon_parent = self._cds_parent(self.refseq_parent, feature='exon')
        self.coords_parent = self._cds_parent(self.refseq_parent,
                                              feature='mRNA', ident='ID')
        self.strand_parent = self._strand_parent(self.refseq_parent,
                                              feature='mRNA')


    def _strand_parent(self, refseq_parent, feature='mRNA'):
        """
        Connect strand and parent id for feature.
        """
        strand_parent = {}
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if not this_line.startswith('#'):
                    this_line = this_line.rstrip('\n').split('\t')
                    chrom = self._refseq_to_common_chrom(this_line[0])
                    if chrom != None:
                        this_feature = this_line[2]
                        strand = this_line[6]
                        info = this_line[8]
                        if chrom not in strand_parent:
                            strand_parent[chrom] = {}
                        if this_feature == feature:
                            this_id = self._find_info_fields(info, 'ID')
                            if this_id in refseq_parent[chrom]:
                                if this_id not in strand_parent[chrom]:
                                    strand_parent[chrom][this_id] = strand
        return strand_parent
                        

    def _info_parent(self, refseq_parent, feature='mRNA', ident='Name'):

        """
        Connect stuff in info field to Parent ID's based on feature and the
        field in the info entry.

        :param refseq_parent:
        :param feature:
        :param ident:
        :return info_parent:
        """

        info_parent = {}
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if not this_line.startswith('#'):
                    this_line = this_line.rstrip('\n').split('\t')
                    chrom = self._refseq_to_common_chrom(this_line[0])
                    if chrom != None:
                        info = this_line[8]
                        this_feature = this_line[2]
                        if chrom not in info_parent:
                            info_parent[chrom] = {}
                        if this_feature == feature:
                            parent = self._find_info_fields(info, 'ID')
                            this_id = self._find_info_fields(info, ident)
                            if parent in refseq_parent[chrom]:
                                if parent not in info_parent[chrom]:
                                    info_parent[chrom][parent] = this_id

        return info_parent


    def _parse_refseq(self, filename):

        """
        Put RefSeq ID's from file in to a list.
        :param filename:
        :return refseq_list:
        """

        refseq_list = []
        with open(filename, 'rU') as refseq:
            for line in refseq:
                line = line.rstrip('\n')
                refseq_list.append(line)

        return refseq_list


    def _cds_parent(self, refseq_parent, feature='CDS', ident='Parent'):

        """
        Get 'feature' sequences, and put them in a dictionary with 'ident'.
        This is in position [8] and follows "Parent=(rna[0-9]+);"

        :param refseq_parent:
        :param feature:
        :param ident:
        :return cds_parent:
        """

        cds_parent = {}
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if not this_line.startswith('#'):
                    this_line = this_line.rstrip('\n').split('\t')
                    chrom = self._refseq_to_common_chrom(this_line[0])
                    if chrom == None:
                        continue
                    info = this_line[8]
                    this_feature = this_line[2]
                    if chrom not in cds_parent:
                        cds_parent[chrom] = {}
                    if this_feature == feature:
                        parent = self._find_info_fields(info, ident)
                        start = this_line[3]
                        stop = this_line[4]
                        if parent in refseq_parent[chrom]:
                            if parent not in cds_parent[chrom]:
                                cds_parent[chrom][parent] = []
                            cds_parent[chrom][parent].append([start, stop])

        return cds_parent


    def _refseq_parent(self, refseq, feature='mRNA', ident='Name'):

        """
        Connect 'feature' id's to 'ident' id's.  Differs from _cds_parent in
        that we are not tracking coordinates here.
        NOTE: Come back to this and refactor based on overlap with _cds_parent.

        :param refseq:
        :param feature:
        :param ident:
        :return refseq_parent:
        """

        refseq_parent = {}
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if not this_line.startswith('#'):
                    this_line = this_line.rstrip('\n').split('\t')
                    chrom = self._refseq_to_common_chrom(this_line[0])
                    if chrom != None:
                        info = this_line[8]
                        this_feature = this_line[2]
                        if chrom not in refseq_parent:
                            refseq_parent[chrom] = {}
                        if this_feature == feature:
                            parent = self._find_info_fields(info, 'ID')
                            refseq_id = self._find_info_fields(info, ident)
                            if refseq_id in refseq:
                                if parent not in refseq_parent[chrom]:
                                    refseq_parent[chrom][parent] = refseq_id

        return refseq_parent


    def _get_header(self):

        """
        Isolate the header section of the GFF3.  This is a list of entries
        beginning with the '#' symbol.
        :return header:
        """

        header = []
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if this_line.startswith('#'):
                    header.append(this_line.rstrip('\n'))
        return header


    def _find_info_fields(self, info, field):

        """
        Find a field from the INFO column.  Return the value based on a
        particular field title.

        :return entry[1]:
        """

        for entry in info.split(';'):
            entry = entry.split('=')
            if entry[0] == field:
                return entry[1]


    def _refseq_to_common_chrom(self, chrom):

        """
        Set dictionary of RefSeq to chromosome identifiers.  Hard code this
        for the time being until I figure out a better way to deal with this
        situation.
        :return refseq_to_chrom[chrom] if in refseq_to_chrom:
        :else return None:
        """

        refseq_to_chrom = {'NC_012920.1': 'MT', 'NC_000001.10': '1',
                           'NC_000002.11': '2', 'NC_000003.11': '3',
                           'NC_000004.11': '4', 'NC_000005.9': '5',
                           'NC_000006.11': '6', 'NC_000007.13': '7',
                           'NC_000008.10': '8', 'NC_000009.11': '9',
                           'NC_000010.10': '10', 'NC_000011.9': '11',
                           'NC_000012.11': '12', 'NC_000013.10': '13',
                           'NC_000014.8': '14', 'NC_000015.9': '15',
                           'NC_000016.9': '16', 'NC_000017.10': '17',
                           'NC_000018.9': '18', 'NC_000019.9': '19',
                           'NC_000020.10': '20', 'NC_000021.8': '21',
                           'NC_000022.10': '22', 'NC_000023.10': 'X',
                           'NC_000024.9': 'Y'}
        try:
            return refseq_to_chrom[chrom]
        except:
            return None


class GffLine(object):
    """
    This is what a typical line looks like.
    NC_000011.9	Gnomon	CDS	533766	533944	.	-	0	ID=cds34112;Name=XP_005252941.1;Parent=rna42823;Dbxref=G
eneID:3265,Genbank:XP_005252941.1,HGNC:5173,HPRD:01813,MIM:190020;gbkey=CDS;gene=HRAS;product=GTPase HRas isoform X1

    NOTE: In progress.

    """

    def __init__(self, gff_entry):

        self.gff_entry = gff_entry.rstrip('\n').split('\t')
        self.refseq_to_chrom = self._refseq_to_common_chrom()
        self.chrom = self.gff_entry[0]
        if self.chrom in self.refseq_to_chrom:
            self.seqname = self.refseq_to_chrom[self.chrom]
        else:
            self.seqname = self.chrom
        self.source = self.gff_entry[1]
        self.feature = self.gff_entry[2]
        self.start = self.gff_entry[3]
        self.end = self.gff_entry[4]
        self.score = self.gff_entry[5]
        self.strand = self.gff_entry[6]
        self.frame = self.gff_entry[7]
        self.grouping = self.gff_entry[8]
