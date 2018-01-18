class BedReader(object):
    """
    Simple class to ingest BED files and return a data structure as such:
    {chrom: [start1, stop1], [start2, stop2], ...}

    Input: filename
    """

    def __init__(self, filename):

        self.filename = open(filename, 'rU')
        self.bed_ints = self._create_bed()

    def _create_bed(self):

        """
        Create the structure of BED coordinates connected to chromosome
        identifiers.
        :return bed_ints:
        """

        bed_ints = {}
        with self.filename as bed:
            for interval in bed:
                interval = interval.rstrip('\n').split('\t')
                chrom = str(interval[0])
                start = int(interval[1])+1 # 0-based
                stop = int(interval[2]) # 1-based

                if chrom not in bed_ints:
                    bed_ints[chrom] = [[start, stop]]
                else:
                    bed_ints[chrom].append([start, stop])
        return bed_ints
