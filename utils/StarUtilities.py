import os, sys, getopt
from numpy import *


def compileFinalLog(StarDir, out, experimentalDesign = None, key = None):
    """Function accepts a STAR output directory and an optional metadata file and row identifier
    in the metadata file

    Args:
        StarDir (str/path): Path to STAR output
        out (str): Compiled STAR data file
        experimentalDesign (str/path): Path to optional metadata file. **The first line of the file must list 'Sample' 
		as a row identifier and the samples in that row correspond to sample directories in StarDir.
        key (str/path): ['all' or rowID for experimentalDesign] Optional key parameter for metadata file.
            If experimental not None, key can be set to 'all' or specific row delimiter in experimental design.

    Examples:
        metadata file format:

        Sample  RNA160304JO_C4_S13_L004_R1_001  RNA160304JO_D2_S12_L004_R1_001  RNA160304JO_G5_S10_L004_R1_001
        AnimalID    700 701 702
        Pressure    Naive   Naive   Naive

     if key == 'all':
        All rows of experimentalDesign will be incorporated into consolidated Star data file
     else:
        key == 'Sample' | 'AnimalID' | 'Pressure'
        Specified row of experimentalDesign will be incorporated into consolidated Star data file

    Returns:
        Compiled STAR log.final.out with optional metadata columns as tab delimited file.
    """
    sample_stats = {}
    column_order = []
    meta_header = []
    directories = [d for d in os.listdir(StarDir) if os.path.isdir(os.path.join(StarDir, d))]
    meta_dic = {}
    if experimentalDesign:
        print('experimentalDesign flag utilized! Using {} as metadata file'.format(experimentalDesign))
        if key != 'all':
            print('Incorportating {} of experimentalDesign into consolidated STAR log.final table: {}'.format(key, out))
            for line in open(experimentalDesign):
                line = line.strip().split()
                if line[0] == 'Sample':
                    samples = line[1:]
                if line[0] == key:
                    keys = samples.copy()
                    values = line[1:]
                    meta_dic = dict(zip(keys, values))
        else:
            for line in open(experimentalDesign):
                line = line.strip().split()
                if line[0] == 'Sample':
                    samples = line[1:]
                else:
                    meta_header.append(line[0])
                    keys = samples.copy()
                    values = line[1:]
                    for i, e in enumerate(keys):
                        if 'Sample' not in meta_dic:
                            meta_dic['Sample'] = {}
                        if e not in meta_dic['Sample']:
                            meta_dic['Sample'][e] = []
                        meta_dic['Sample'][e].append(values[i])

    if experimentalDesign:
        if key == 'all':
            assert (meta_dic[
                        'Sample'].keys().sort() == directories.sort()), 'Samples listed in experimentalDesign are not concordant with STAR directories'
        else:
            assert (
                        meta_dic.keys().sort() == directories.sort()), 'Samples listed in experimentalDesign are not concordant with STAR directories'

    for d in directories:
        FL = os.path.join(StarDir, d + '/' + d + '_Log.final.out')
        for line in open(FL):
            if '|' in line:
                line = line.strip().split('|')
                line = [x.strip() for x in line]
                if d not in sample_stats:
                    sample_stats[d] = {}
                if line[0] not in sample_stats[d]:
                    sample_stats[d][line[0]] = line[1]
                if line[0] not in column_order:
                    column_order.append(line[0])
            else:
                continue

    samples = sample_stats.keys()
    mat = chararray((len(sample_stats), len(column_order)), itemsize = 50)

    for i, e in enumerate(samples):
        for j, m in enumerate(column_order):
            try:
                mat[i, j] = sample_stats[e][m]
            except:
                mat[i, j] = 'NA'

    meta_list = []
    if experimentalDesign:
        if key != 'all':
            for i, e in enumerate(samples):
                meta_list.append(meta_dic[e])
        else:
            for i, e in enumerate(samples):
                meta_list.append(meta_dic['Sample'][e])

    if experimentalDesign:
        if key != 'all':
            meta_list.insert(0, key)
            column_order.insert(0, 'Samples')
            row_mat = column_stack((samples, mat))
            col_mat = row_stack((column_order, row_mat))
            col_mat = column_stack((meta_list, col_mat))
            savetxt(out, col_mat, fmt = '%s', delimiter = '\t')
        else:
            meta_list.insert(0, meta_header)
            meta_array = array(meta_list)
            column_order.insert(0, 'Samples')
            row_mat = column_stack((samples, mat))
            col_mat = row_stack((column_order, row_mat))
            col_mat = column_stack((meta_array, col_mat))
            savetxt(out, col_mat, fmt = '%s', delimiter = '\t')
    else:
        column_order.insert(0, 'Samples')
        row_mat = column_stack((samples, mat))
        col_mat = row_stack((column_order, row_mat))
        savetxt(out, col_mat, fmt = '%s', delimiter = '\t')


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv = None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print(compileFinalLog.__doc__)
                sys.exit(0)
    except Usage as err:
        print(sys.stderr, err.msg)
        print(sys.stderr, "for help use --help")
        return 2
    if len(args) == 2:
        print(args, '2 arguments')
        compileFinalLog(args[0], args[1])
    else:
        print(args)
        try:
            print(args, '4 arguments')
            compileFinalLog(args[0], args[1], args[2], args[3])
        except:
            if len(args) == 4:
                print('Arguments passed incorrectly: ', zip(['StarDir', 'out', 'experimentalDesign', 'key'], args))
            else:
                print('Additional Arguments required')


if __name__ == "__main__":
    main()
