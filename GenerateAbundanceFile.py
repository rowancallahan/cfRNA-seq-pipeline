# coding=utf-8
"""Frameworks for applying QC functions to -omics datasets.

This file contains functions used to generate SLURM submission scripts,
generate Rscripts to assess quality of omic datasets.

Author: Joey Estabrook <estabroj@exacloud.ohsu.edu>

"""

import os
import sys
import textwrap
import argparse
import subprocess
import stat
from operator import itemgetter
from abundance_models import *


def ensure_dir(relnm):
    """ Accept relative filepath string, create it if it doesnt already exist
        return filepath string

    Args:
        relnm (str) : Relative name/path

    Returns:
        relnm (str)

    """

    if not os.path.exists(relnm):
        os.makedirs(relnm)

    return relnm


def generate_meta_file(read_dir, sample_meta_data_list, select_meta_data_list, split_hyphen = None):
    """Generate meta file from sample directory by column orientation based on IGL data formatting

    Args:
        split_hyphen (NoneType): Split IGL header and create meta file by respective field hyphen split represents
        read_dir (str): Abs path to data directory
        sample_meta_data_list (str): list representing all fields of IGL sample label
        select_meta_data_list (str): list representing fields of IGL sample label to subset by
    """
    sample_meta_data_list = sample_meta_data_list.strip().split(',')
    select_meta_data_list = select_meta_data_list.strip().split(',')

    check = set(sample_meta_data_list).issuperset(set(select_meta_data_list))
    if check:

        files = [f for f in os.listdir(read_dir) if f.startswith('RNA')]
        files.sort()
        project = files[0].split('_')[0]
        out_f = open(project + '_metadata.txt', 'w')

        if split_hyphen:
            out_f.write('\t'.join(select_meta_data_list) + '\t' + 'ID_Group' + '\n')
        else:
            out_f.write('\t'.join(select_meta_data_list) + '\n')

        idx = [sample_meta_data_list.index(i) for i in select_meta_data_list]
        for f in files:
            metadata = f.split('_')
            pertinent = list(itemgetter(*idx)(metadata))
            if split_hyphen:
                hyphen_idx = [x for x in pertinent if '-' in x]
                if len(hyphen_idx) > 1:
                    print('Multiple meta fields containing hyphens!')
                else:
                    pertinent.append(hyphen_idx[0].split('-')[0])
                    out_f.write('\t'.join(pertinent) + '\n')
            else:
                out_f.write('\t'.join(pertinent) + '\n')
        out_f.close()
    else:
        for i in select_meta_data_list:
            if i not in sample_meta_data_list:
                print(i)
        sys.exit(1)


def generate_slurm_submit_script(project_title, read_dir):
    """Generates SLURM submission script
    
    Args:
        project_title (str): Label for experiment directory and plot titles
        read_dir (str): Abs path to data directory
    """
    log_files = os.path.join(read_dir, '/logs')
    analysis_sub = os.path.join(read_dir, '/analysis_code')
    out_f = open(os.path.join(analysis_sub, project_title + '_analysis.submit'), 'w')
    cmd = """
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --time=34:00:00
#SBATCH --job-name={project_title}
#SBATCH --mem=30000
#SBATCH --error={log_files}/{project_title}.%j.AutoGenerate.stderr
#SBATCH --output={log_files}/{project_title}.%j.AutoGenerate.stdout

source activate Crayon

Rscript {analysis_sub}/{project_title}_analysis.R
"""
    reformatted_cmd = textwrap.dedent(cmd).strip()
    context = {"project_title": project_title, "log_files": log_files, "analysis_sub": analysis_sub}
    out_f.write(reformatted_cmd.format(**context))
    out_f.close()


def launch_slurm_submit_script(project_title, read_dir):
    """Launch generated SLURM script

    Args:
        project_title (str): Label for experiment directory and plot titles
        read_dir (str): Abs path to data directory
    """
    analysis_sub = os.path.join(read_dir, '/analysis_code')
    sub_script = "{}_analysis.submit".format(project_title)
    analysis_script = "{}_analysis.R".format(project_title)
    scripts = [sub_script, analysis_script]
    cmd = "sbatch {}/{}_analysis.submit".format(analysis_sub, project_title)
    # Set permissions on scripts to be launched to scheduler

    current_permissions = stat.S_IMODE(os.lstat(analysis_sub).st_mode)
    for s in scripts:
        os.chmod(os.path.join(analysis_sub, s), current_permissions)

    subprocess.call(cmd, shell = True)


def generate_abundance_script(read_dir, meta_file, code_dir, tax_id, gtf_file, project_title, baseline,
                              sample_id = "sample_id", mart_dataset = "hsapiens_gene_ensembl", lm_by = "ID_Group",
                              gtf_feature = "gene", read_pattern = "*", useme_cols = "*", label_from_colname = "*",
                              path_type = "gene.counts", path_norms = "loess", covariate = False, ann_colplotme = None,
                              load_table = False, dataset_path = None, deseq = False, logged_b = 'F'):
    """Generate QC -omics dataset Rscript

    Args:
        read_dir (str): Abs path to out directory
        meta_file (str): Abs path to associated metafile for -omics dataset
        code_dir (str): Abs path to code directory
        tax_id (int): Taxa ID
        gtf_file (str): Abs path to gtf file for associated omics dataset samples
        project_title (str): Label for experiment directory and plot titles
        baseline (str): Baseline factor level to generate linear model
        sample_id (str): Column header in [meta_file] with sample ids to be reported in results table
        mart_dataset (str): Query for Biomart to map genomic feature id to ensemble or Hugo gene id
        lm_by (str): Column header in [meta_file] with factor level [baseline]
        gtf_feature (str): Feature to filter gtf file by
        path_norms (str): Normalization technique implemented
        covariate (bool): Boolean to implement covariates in model
        ann_colplotme (NoneType|str): defaults to {lm_by} or specified as str of format "c('column[i]','column[i+1]')"
        load_table (bool): Boolean to load dataset matrix [n_features, n_samples]
        dataset_path (str): Absolute path to dataset matrix
        deseq (bool): Boolean to generate DESeq2 fitted values to be plotted in QC/QA analysis
        logged_b (bool): Boolean to log2 scale data before processing
        path_type (str): STAR counts to quantify - defaults to gene.counts
        label_from_colname (str): regex expression to identify samples to perform QC/QA analysis on
        useme_cols (str): regex expression to incorporate meta data into expression object
        read_pattern (str): regex expression to identify samples to perform QC/QA analysis on

    """
    if ann_colplotme is None:
        ann_colplotme = lm_by

    gtf_read_dir = '/'.join(gtf_file.split('/')[:-1])
    analysis_sub = os.path.join(read_dir,'analysis_code')
    results = os.path.join(read_dir, 'results')
    log_files = os.path.join(read_dir, 'logs')

    print(analysis_sub)
    print(results)
    print(log_files)

    _ = ensure_dir(analysis_sub)
    _ = ensure_dir(results)
    _ = ensure_dir(log_files)

    out_f = open(os.path.join(analysis_sub, project_title + '_analysis.R'), 'w')

    code_context = {"code_dir": code_dir, "meta_file": meta_file, "sample_id": sample_id, "tax_id": tax_id,
                    "gtf_file": gtf_file, "gtf_feature": gtf_feature, "project_title": project_title,
                    "gtf_read_dir": gtf_read_dir, "read_dir": read_dir, "read_pattern": read_pattern,
                    "useme_cols": useme_cols, "lm_by": lm_by, "baseline": baseline, "path_type": path_type,
                    "path_norms": path_norms, "mart_dataset": mart_dataset, "label_from_colname": label_from_colname,
                    "ann_colplotme": ann_colplotme, "results": results, "dataset": dataset_path, "logged_b": logged_b}

    if not covariate:
        contrast_str = """contr_ls = list("{lm_by}" = list(baseline="{baseline}", contr.FUN = "contr.treatment"))"""\
            .format(**code_context)
        lm_expr = "y ~ {lm_by}".format(**code_context)
        code_context['lm_expr'] = lm_expr
        code_context['contr_ls'] = contrast_str
        code_context['ann_colplotme'] = '{}'.format(ann_colplotme)
        code_context['annCollm_by'] = '"{}"'.format(lm_by)
        code_context['oneclass'] = '"{}"'.format(lm_by)

    else:

        if len(covariate.split(',')) != len(baseline.split(',')):
            print("""Provide co-variates and their associated baselines in ordered comma delimited list .ie

            covariate = 'Time,Pressure'
            baseline = '0Hr,0mmHg'

            resulting in:

            contr_ls = list(Time = list(baseline="0Hr", contr.FUN = "contr.treatment"), 
                Pressure = list(baseline="0mmHg", contr.FUN = "contr.treatment"))
            
            & 
            
            lm_expr = 'y ~ Time + Pressure

            ***Note***
            Both Time and Pressure need to be column headers in the table provided as an argument 
            to expt.design in regressMatrix covariate list provided:
            """)
            print(covariate.split(','))
            print("baseline list provided:")
            print(baseline.split(','))
        else:
            code_context['oneclass'] = '"{}"'.format(lm_by)
            covariate_list = covariate.split(',')
            baseline_list = baseline.split(',')
            lm_expr = "y ~ " + ' + '.join(covariate_list)
            code_context['lm_expr'] = lm_expr
            cov_str = ""
            for i, j in zip(covariate_list, baseline_list):
                cov_str += """'{}' = list(baseline='{}', contr.FUN = 'contr.treatment'),""".format(i, j)
            covariate_str = "list(" + cov_str[:-1] + ")"
            code_context['contr_ls'] = covariate_str
            reformat_covariate = '"' + '","'.join(covariate_list) + '"'
            code_context['ann_colplotme'] = 'c({})'.format(reformat_covariate)
            code_context['annCollm_by'] = 'c({})'.format(reformat_covariate)

    if not load_table:
        code = qc_model
    else:
        if deseq:
            code = qc_matrix_w_deseq_model
        else:
            code = qc_matrix_model
    reformatted_code = textwrap.dedent(code).strip()
    out_f.write(reformatted_code.format(**code_context))
    out_f.close()


def main():
    parser = argparse.ArgumentParser(description = 'Generate Abundance Workflow Wrapper',
                                 usage='use "python GeneratAbundanceFile.py --help" for more information',
                                 formatter_class = argparse.RawTextHelpFormatter)
    metagroup = parser.add_argument_group('Mandatory Metadata arguments')
    metagroup.add_argument("-mb", "--meta",
                           action = "store_true",
                           help = "Boolean to determine whether to make a meta file",
                           default = False)
    metagroup.add_argument("-md", "--sample_meta_data_list",
                           type = str,
                           metavar='',
                           help = """
Comma delimited list of meta information provided by underscore delimited sample directory name 
i.e RNA160606DM_294-1_S2_L001_R1_001 = 'Project','ID','sample','lane','r1','01'""")
    metagroup.add_argument("-ms", "--select_meta_data_list",
                           type = str,
                           metavar='',
                           help = """
subset of the comma delimited list of meta information provided by underscore delimited sample directory name
i.e 'ID','sample': These fields will be incorporated into metadata.txt""")
    metagroup.add_argument("-sh", "--split_hyphen",
                           action = "store_true",
                           default = False,
                           help = """
Split and incorporate hyphenated meta data field into meta data i.e. 294-1 in 
RNA160606DM_294-1_S2_L001_R1_001 will be incorporated into metadata.txt under 
the column ID_Group as 294""")
    abundancegroup = parser.add_argument_group('Mandatory Abundance script generation arguments')
    abundancegroup.add_argument("-d", "--read_dir",
                            type = str,
                            help = "Absolute path to abundance data directory i.e ../STAR/",
                            metavar="")
    abundancegroup.add_argument("-mf", "--meta_file",
                            type = str,
                            metavar="",
                            help = "Absolute path to metafile.txt generated via .generate_meta_file")
    abundancegroup.add_argument("-p", "--project_title",
                            type = str,
                            metavar="",
                            help = """
Project title associated with abundance dataset.(Will be incorporated into base file name of plots)""")
    abundancegroup.add_argument("-b", "--baseline",
                            type = str,
                            metavar="",
                            default = False,
                            help = "Baseline to generate contrasts against when generating lm")
    optabundancegroup = parser.add_argument_group('Optional Abundance script generation arguments')
    optabundancegroup.add_argument("-lj", "--launch_job",
                            action = "store_true",
                            help = "Generate SLURM submission script and launch job",
                            default = False)
    optabundancegroup.add_argument("-df", "--load_table",
                            action = "store_true",
                            help = "Loads appropriate abundance model to read in matrices",
                            default = False)
    optabundancegroup.add_argument("-ds", "--deseq",
                            action = "store_true",
                            help = "Loads appropriate abundance model to read in matrices with deseq",
                            default = False)
    optabundancegroup.add_argument("-c", "--code_dir",
                            type = str,
                            metavar = '',
                            help = "Absolute path to code directory i.e. ProjectDirectory/code/",
                            default = "/home/exacloud/lustre1/CEDAR/roskamsh/tools/omic_resources/")
    optabundancegroup.add_argument("-t", "--tax_id",
                            type = str,
                            metavar = '',
                            help = """
tax_id can be found www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi""",
                            default = "9606")
    optabundancegroup.add_argument("-g", "--gtf_file",
                            type = str,
                            metavar = '',
                            help = "Absolute path to gtf used for alignment",
                            default = """/home/exacloud/lustre1/BioCoders/DataResources/Genomes/hg19/release-75/gtf/Homo_sapiens.GRCh37.75.gtf""")
    optabundancegroup.add_argument("-lb", "--logged_B",
                            type = str,
                            metavar='',
                            help = "R Boolean flag to log2 raw matrix F|T",
                            default = "F")
    optabundancegroup.add_argument("-da", "--data_file_path",
                            type = str,
                            metavar='',
                            help = "Absolute path to dataset matrix")
    optabundancegroup.add_argument("-co", "--covariate",
                            metavar='',
                            type = str,
                            help = "Covariates to generate contrasts against when generating lm")
    optabundancegroup.add_argument("-id", "--sample_id",
                            metavar='',
                            type = str,
                            help = """
Column in metadata.txt that identifies samples read in by label_from_colname""",
                            default = "sample_id")
    optabundancegroup.add_argument("-bm", "--mart_dataset",
                            metavar='',
                            type = str,
                            help =
"""BiomaRt dataset to use. Datasets include: mmusculus_gene_ensembl | hsapiens_gene_ensembl | 
aplatyrhynchos_gene_ensembl | drerio_gene_ensembl | ggallus_gene_ensembl | oaries_gene_ensembl | 
rnorvegicus_gene_ensembl | sscrofa_gene_ensembl
""",
                            default = "hsapiens_gene_ensembl")
    optabundancegroup.add_argument("-lm", "--lmBy",
                            metavar='',
                            type = str,
                            help = "Column name in metadata.txt that contains [baseline] value",
                            default = "ID_Group")
    optabundancegroup.add_argument("-gf", "--gtf_feature",
                            metavar='',
                            type = str,
                            help = "GTF feature to annotate abundance and ratio tables with",
                            default = "gene")
    optabundancegroup.add_argument("-rp", "--read_pattern",
                            metavar='',
                            type = str,
                            help = """
Read pattern expression provided to R to read in sample associated abundance information""",
                            default = "*")
    optabundancegroup.add_argument("-uc", "--useme_cols",
                            metavar='',
                            type = str,
                            help = """
Read pattern expression provided to R to select data to be incorporated in STAR.data""",
                            default = "*")
    optabundancegroup.add_argument("-lc", "--label_from_colname",
                            metavar='',
                            type = str,
                            help = """
Read pattern expression provided to R to select for unique sample label identifiers""",
                            default = "*")
    optabundancegroup.add_argument("-pt", "--path_type",
                            metavar='',
                            type = str,
                            help = "gene.counts | SJ.counts",
                            default = "gene.counts")
    optabundancegroup.add_argument("-pn", "--path_norms",
                            metavar='',
                            type = str,
                            help = "alograw | loess | lowess | qspln | quant",
                            default = "loess")
    optabundancegroup.add_argument("-pl", "--plot_me",
                            metavar='',
                            type = str,
                            help = "Column headers in meta file",
                            default = None)

    args = parser.parse_args()

    if args.meta:
        generate_meta_file(args.read_dir, args.sample_meta_data_list, args.select_meta_data_list, args.split_hyphen)

    else:
        generate_abundance_script(args.read_dir, args.meta_file, args.code_dir, args.tax_id, args.gtf_file,
                              args.project_title, args.baseline, args.sample_id, args.mart_dataset, args.lmBy,
                              args.gtf_feature, args.read_pattern, args.useme_cols, args.label_from_colname,
                              args.path_type, args.path_norms, args.covariate, args.plot_me, args.load_table,
                              args.data_file_path, args.deseq, args.logged_B)
        if args.launch_job:
            generate_slurm_submit_script(args.project_title, args.read_dir)
            launch_slurm_submit_script(args.project_title, args.read_dir)


if __name__ == '__main__':
    main()
