# coding=utf-8
import os
import sys
import textwrap
import argparse
import subprocess
import stat
from operator import itemgetter
from abundance_models import *


def generate_meta_file(read_dir, sample_meta_data_list, select_meta_data_list, split_hyphen=None):
    """
    Generates a meta_file
    by column orientation
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
                print (i)
        sys.exit(1)


def generate_slurm_submit_script(project_title, read_dir):
    log = '/'.join(read_dir.split('/')[:-3])+'/logs'
    submit = '/'.join(read_dir.split('/')[:-2]) + '/analysis_code'
    out_f = open(os.path.join(submit, project_title + '_analysis.submit'), 'w')
    cmd = """
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --time=34:00:00
#SBATCH --job-name={project_title}
#SBATCH --mem=30000
#SBATCH --error={log}/{project_title}.%j.AutoGenerate.stderr
#SBATCH --output={log}/{project_title}.%j.AutoGenerate.stdout

source activate python2

Rscript {submit}/{project_title}_analysis.R
"""
    reformatted_cmd = textwrap.dedent(cmd).strip()
    context = {"project_title": project_title, "log" : log, "submit":submit}
    out_f.write(reformatted_cmd.format(**context))
    out_f.close()


def launch_slurm_submit_script(project_title, read_dir):
    log = '/'.join(read_dir.split('/')[:-2])+'/analysis_code'
    sub_script = "%s_analysis.submit" % (project_title)
    analysis_script = "%s_analysis.R" % (project_title)
    scripts = [sub_script, analysis_script]
    cmd = "sbatch %s/%s_analysis.submit" % (log, project_title)
    # Set permissions on scripts to be launched to scheduler

    current_permissions = stat.S_IMODE(os.lstat(log).st_mode)
    for s in scripts:
        os.chmod(os.path.join(log, s), current_permissions)

    subprocess.call(cmd, shell=True)


def generate_abundance_script(read_dir, meta_file, code_dir, taxID, gtfFile, project_title,
                              baseline, SampleID="SampleID", mart_dataset="hsapiens_gene_ensembl",
                              lmBy="ID_Group", gtf_feature="gene", read_pattern="*", useme_cols="*",
                              label_from_colname="*", path_type="gene.counts",
                              path_norms="loess", covariate=False, annColplotme=None, load_table=False, dataset_path=None,
                              deseq=False):
    if annColplotme == None:
        annColplotme = lmBy

    gtf_read_dir = '/'.join(gtfFile.split('/')[:-1])
    log = '/'.join(read_dir.split('/')[:-2])+'/analysis_code'
    out_f = open(os.path.join(log, project_title + '_analysis.R'), 'w')
    results = '/'.join(read_dir.split('/')[:-2])+'/results'
    code_context = {"code_dir": code_dir, "meta_file": meta_file, "SampleID": SampleID, "taxID": taxID,
                    "gtfFile": gtfFile, "gtf_feature": gtf_feature, "project_title": project_title,
                    "gtf_read_dir": gtf_read_dir, "read_dir": read_dir, "read_pattern": read_pattern,
                    "useme_cols": useme_cols, "lmBy": lmBy, "baseline": baseline, "path_type": path_type,
                    "path_norms": path_norms, "mart_dataset": mart_dataset, "label_from_colname": label_from_colname,
                    "annColplotme": annColplotme, "results":results, "dataset":dataset_path}

    if not covariate:
        contrast_str = """contr_ls = list("{lmBy}" = list(baseline="{baseline}", contr.FUN = "contr.treatment"))""".format(**code_context)
        lm_expr = "y ~ {lmBy}".format(**code_context)
        code_context['lm_expr'] = lm_expr
        code_context['contr_ls'] = contrast_str
        code_context['annColplotme'] = '%s' % annColplotme
        code_context['annCollmBy'] = '"%s"' % lmBy
        code_context['oneclass'] = '"%s"' % lmBy

    else:

        if len(covariate.split(',')) != len(baseline.split(',')):
            print("""Provide co-variates and their associated baselines in ordered comma delimited list .ie

            covariate = 'Time,Pressure'
            baseline = '0Hr,0mmHg'

            resulting in:

            contr_ls = list(Time = list(baseline="0Hr", contr.FUN = "contr.treatment"), Pressure = list(baseline="0mmHg", contr.FUN = "contr.treatment"))

            &

            lm_expr = 'y ~ Time + Pressure

            ***Note***
            Both Time and Pressure need to be column headers in the table provided as an argument to expt.design in regressMatrix

            covariate list provided:
            """)
            print(covariate.split(','))
            print("baseline list provided:")
            print(baseline.split(','))
        else:
            code_context['oneclass'] = '"%s"' % lmBy
            covariate_list = covariate.split(',')
            baseline_list = baseline.split(',')
            lm_expr = "y ~ " +' + '.join(covariate_list)
            code_context['lm_expr'] = lm_expr
            cov_str = ""
            for i in range(len(zip(covariate_list,baseline_list))):
                cov_str += """'%s' = list(baseline='%s', contr.FUN = 'contr.treatment'),"""%(covariate_list[i],baseline_list[i])
            covariate_str = "list(" + cov_str[:-1] + ")"
            code_context['contr_ls'] = covariate_str
            reformat_covariate = '"'+'","'.join(covariate_list)+'"'
            code_context['annColplotme'] = 'c(%s)'%reformat_covariate
            code_context['annCollmBy'] = 'c(%s)'%reformat_covariate

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
    parser = argparse.ArgumentParser(description='Generate Abundance Workflow Wrapper')
    parser.add_argument("-d", "--read_dir", type=str, help="Absolute path to abundance data directory i.e ../STAR/")
    metagroup = parser.add_argument_group('Mandatory Metadata arguments')
    metagroup.add_argument("-mb", "--meta", action="store_true", dest="meta", help="Boolean to determine whether to make a meta file", default=False)
    metagroup.add_argument("-md", "--sample_meta_data_list", type=str, help="comma delimited list of meta information provided by underscore delimited sample directory name i.e RNA160606DM_294-1_S2_L001_R1_001 = 'Project','ID','sample','lane','r1','01'")
    metagroup.add_argument("-ms", "--select_meta_data_list", type=str, help="subset of the comma delimited list of meta information provided by underscore delimited sample directory name i.e 'ID','sample': These fields will be incorporated into metadata.txt")
    metagroup.add_argument("-sh", "--split_hyphen", action="store_true", dest="split_hyphen", help="split and incorporate hyphenated meta data field into meta data i.e. 294-1 in RNA160606DM_294-1_S2_L001_R1_001 will be incorporated into metadata.txt under the column ID_Group as 294", default=False)

    abundancegroup = parser.add_argument_group('Mandatory Abundance script generation arguments')
    abundancegroup.add_argument("-lj", "--launch_job", action="store_true", dest="launch_job", help="Generate SLURM submission script and launch job", default=False)
    abundancegroup.add_argument("-df", "--load_table", action="store_true", dest="load_table", help="Loads appropriate abundance model to read in matrices", default=False)
    abundancegroup.add_argument("-ds", "--deseq", action="store_true", dest="deseq", help="Loads appropriate abundance model to read in matrices with deseq", default=False)
    abundancegroup.add_argument("-mf", "--meta_file", type=str, help="Absolute path to metafile.txt generated via .generate_meta_file")
    abundancegroup.add_argument("-c", "--code_dir", type=str, help="Absolute path to code directory i.e. ProjectDirectory/code/", default="/home/users/estabroj/scratch/CEDAR/new_repo/")
    abundancegroup.add_argument("-t", "--taxID", type=str, help="TaxID can be found www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi", default="9606")
    abundancegroup.add_argument("-g", "--gtfFile", type=str, help="Absolute path to gtf used for alignment", default="/home/exacloud/lustre1/BioCoders/DataResources/Genomes/hg19/release-75/gtf/Homo_sapiens.GRCh37.75.gtf")
    abundancegroup.add_argument("-p", "--project_title", type=str, help="Project title associated with abundance dataset. (Will be incorporated into base file name of plots)")
    abundancegroup.add_argument("-b", "--baseline", type=str, help="Baseline to generate contrasts against when generating lm", default=False)

    optabundancegroup = parser.add_argument_group('Optional Abundance script generation arguments')
    optabundancegroup.add_argument("-da", "--data_file_path", type=str, help="Absolute path to dataset matrix")
    optabundancegroup.add_argument("-co", "--covariate", type=str, help="Covariates to generate contrasts against when generating lm")
    optabundancegroup.add_argument("-id", "--SampleID", type=str, help="Column in metadata.txt that identifies samples read in by label_from_colname", default="SampleID")
    optabundancegroup.add_argument("-bm", "--mart_dataset", type=str, help="biomaRt dataset to use. Datasets include: mmusculus_gene_ensembl | hsapiens_gene_ensembl | aplatyrhynchos_gene_ensembl | drerio_gene_ensembl | ggallus_gene_ensembl | oaries_gene_ensembl | rnorvegicus_gene_ensembl | sscrofa_gene_ensembl", default="hsapiens_gene_ensembl")
    optabundancegroup.add_argument("-lm", "--lmBy", type=str, help="Column name in metadata.txt that contains [baseline] value", default="ID_Group")
    optabundancegroup.add_argument("-gf", "--gtf_feature", type=str, help="GTF feature to annotate abundance and ratio tables with", default="gene")
    optabundancegroup.add_argument("-rp", "--read_pattern", type=str, help="Read pattern expression provided to R to read in sample associated abundance information", default="*")
    optabundancegroup.add_argument("-uc", "--useme_cols", type=str, help="Read pattern expression provided to R to select data to be incorporated in STAR.data ", default="*")
    optabundancegroup.add_argument("-lc", "--label_from_colname", type=str, help="Read pattern expression provided to R to select for unique sample label identifiers", default="*")
    optabundancegroup.add_argument("-pt", "--path_type", type=str, help="gene.counts | SJ.counts", default="gene.counts")
    optabundancegroup.add_argument("-pn", "--path_norms", type=str, help=" alograw | loess | lowess | qspln | quant", default="loess")
    optabundancegroup.add_argument("-pl", "--plot_me", type=str, help="Column headers in meta file", default=None)

    args = parser.parse_args()
    if args.meta:
        generate_meta_file(args.read_dir, args.sample_meta_data_list, args.select_meta_data_list, args.split_hyphen)
    else:
        generate_abundance_script(args.read_dir, args.meta_file, args.code_dir, args.taxID, args.gtfFile,
                                  args.project_title,
                                  args.baseline, args.SampleID, args.mart_dataset,
                                  args.lmBy, args.gtf_feature, args.read_pattern, args.useme_cols,
                                  args.label_from_colname, args.path_type,
                                  args.path_norms, args.covariate, args.plot_me, args.load_table, args.data_file_path, args.deseq)
        if args.launch_job:
            generate_slurm_submit_script(args.project_title, args.read_dir)
            launch_slurm_submit_script(args.project_title, args.read_dir)

if __name__ == '__main__':
    main()

