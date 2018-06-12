__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

configfile: "omic_config.yaml"

import sys
import yaml

def format_plot_columns({config[meta_columns_to_plot]}):
    factors = {config[meta_columns_to_plot]}.split(',')
    reformat_factors = '"' + '","'.join(factors) + '"'
    return 'c({})'.format(reformat_factors)

rule qc_qa:
 input:
    "{config[omic_counts_data]}/{config[project_id]}_counts.txt"
 output:
    "{config[omic_qc_results]/{config[project_id]}/tables/{project_id}_loess_Normed_with_Ratio_and_Abundance.txt"
 log:
    "{config[omic_qc_results]}/{config[project_id]}/logs"
 shell:
    python GenerateAbundanceFile.py -d {config[results}]/{config[project_id]} -mf {omic_meta_data} -p {config[project_id]} -b {config[baseline]} -lm {config[linear_model]} -id {config[StudyID]} -pl format_plot_columns({config[meta_columns_to_plot]}) -df -da {config[omic_counts_data]}

