import argparse
import os
import numpy as np


# noinspection PyTypeChecker
def compile_final_log(star_dir, project_title):
    """Function accepts a STAR output directory and an optional metadata file and row identifier
    in the metadata file

    Args:
        star_dir (str/path): Path to STAR output
        project_title (str): Project title for compiled STAR mapping statistics

    Returns:
        Compiled STAR log.final.out as tab delimited file.
    """

    sample_stats = {}
    column_order = []
    directories = [d for d in os.listdir(star_dir) if os.path.isdir(os.path.join(star_dir, d))]

    for d in directories:
        fh = os.path.join(star_dir, d + '/' + d + '_Log.final.out')
        for line in open(fh):
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
    mat = np.chararray((len(sample_stats), len(column_order)), itemsize = 50)

    for i, e in enumerate(samples):
        for j, m in enumerate(column_order):
            try:
                mat[i, j] = sample_stats[e][m]
            except KeyError:
                mat[i, j] = 'NA'

    column_order.insert(0, 'Samples')
    row_mat = np.column_stack((samples, mat))
    col_mat = np.row_stack((column_order, row_mat))
    out = "{}_STAR_mapping_statistics.txt".format(project_title)
    np.savetxt(out, col_mat, fmt = '%s', delimiter = '\t')


def main():
    parser = argparse.ArgumentParser(description = 'Generate STAR mapping statistics summary table',
                                     usage = 'use "python StarUtilities.py --help" for more information',
                                     formatter_class = argparse.RawTextHelpFormatter)
    metagroup = parser.add_argument_group('Mandatory Compile STAR summary table arguments')
    metagroup.add_argument("-d", "--read_dir", type = str,
                           help = "Absolute path to abundance data directory i.e ../STAR/", metavar = "")
    metagroup.add_argument("-p", "--project_title", type = str, metavar = "",
                           help = """Project title associated with STAR summary table""")

    args = parser.parse_args()
    compile_final_log(args.read_dir, args.project_title)


if __name__ == "__main__":
    main()
