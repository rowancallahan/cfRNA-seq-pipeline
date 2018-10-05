import argparse
import os
import pandas as pd


def compile_star_log(files, out):
    """Function accepts a STAR output directory and compiles all sample information from Log.final.out

    Args:
        data_dir (str/path): Path to STAR output
        project_title (str): Project title for compiled STAR mapping statistics

    Returns:
        Compiled STAR log.final.out as tab delimited file.
    """

    tables = [pd.read_csv(fh, sep = '\t', index_col = 0, names = [fh.split('/')[-2]]) for fh in files]
    joined_table = pd.concat(tables, axis=1)
    joined_table_sorted = joined_table.reindex(sorted(joined_table.columns), axis = 1)
    joined_table_sorted.to_csv(out, sep='\t')


def compile_counts_table(files, project_title):
    """Function accepts a STAR output directory and compiles {sample}_gene_count.txt into a joined tab sep file

    Args:
        files (list): list of globbed wildcards
        project_title (str): Project title for compiled STAR gene counts table

    Returns:
        Compiled STAR gene counts table as tab delimited file.
    """
    tables = [pd.read_csv(fh, sep='\t', index_col=0, names=[fh.split('/')[-1].split('_')[0]]) for fh in files]
    joined_table = pd.concat(tables, axis=1)
    filtered_joined = joined_table.iloc[:-5, :]
    filtered_joined_sorted = filtered_joined.reindex(sorted(filtered_joined.columns), axis = 1)
    out = "data/{}_counts.txt".format(project_title)
    filtered_joined_sorted.to_csv(out, sep='\t')


def main():
    parser = argparse.ArgumentParser(description = 'Generate STAR mapping statistics summary table',
                                     usage = 'use "python StarUtilities.py --help" for more information',
                                     formatter_class = argparse.RawTextHelpFormatter)
    metagroup = parser.add_argument_group('Mandatory Compile STAR summary table arguments')
    metagroup.add_argument("-d", "--read_dir", type = str,
                           help = "Absolute path to abundance data directory i.e ../STAR/", metavar = "")
    metagroup.add_argument("-p", "--project_title", type = str, metavar = "",
                           help = """Project title associated with STAR summary table""")
    metagroup.add_argument("-c", "--compile_counts",
                           help = "Compile counts table", action = "store_true", default = False)
    args = parser.parse_args()

    if args.compile_counts:
        compile_counts_table(args.read_dir, args.project_title)
    else:
        compile_star_log(args.read_dir, args.project_title)


if __name__ == "__main__":
    main()
