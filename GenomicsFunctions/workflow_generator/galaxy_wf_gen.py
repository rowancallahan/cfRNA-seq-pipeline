#!/usr/bin/env python

import subprocess
from bioblend import galaxy
import os
import uuid
import time
import getpass
import shutil
import argparse
import json
import samplesheet_class

import requests.packages.urllib3
requests.packages.urllib3.disable_warnings()

##### Need to send parameters to the cgd_client wrapper.  Can set some of the stuff ahead of time, but barcodeid, especially, needs to be set on the fly.

### Functions interacting with the shell.
def find_env_var(env_var):
    """
    Return environmental variable, if it exists.
    """

    try:
        return os.environ[env_var]
    except:
        raise KeyError("Environemental variable " + env_var + "not found.")


### Concatenate the read 1's and read 2's per sample id, passed from the shell script.
def catFastq(proj_dir, sample, analysis_type, dataset_import_dir, r1_file,
             r2_file):
    """
    Input: proj_dir (system location of proj_dir)
    """

    output_prefix = proj_dir + sample
    analysis_type = analysis_type.replace(" ", "\ ")
#    r1_file = dataset_import_dir + sample + "_R1.fastq.gz"
#    r2_file = dataset_import_dir + sample + "_R2.fastq.gz"

    cmd = "cat " + proj_dir + analysis_type + "/*" + sample + "*R1* >> " + r1_file
    print(cmd)
    p = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
    os.chmod(r1_file, 0444)

    cmd = "cat " + proj_dir + analysis_type + "/*" + sample + "*R2* >> " + r2_file
    print(cmd)
    p = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
    os.chmod(r2_file, 0444)


### Functions specifically creating structure to feed in to the BioBlend API.
def buildParams(gi, read_group_dict, workflow_id, sample_id, history_id,
                runid, barcodeid, tissue, param_set, bam_path, cgd_client):

    """
    Create parameter dictionary.  This will supply tools that need input at runtime with the
    necessary values.

    TODO: Refactor, too many input parameters.
    """

    description = "Automated variant calling workflow."

    for key in gi.workflows.show_workflow(workflow_id)['steps']:
            if gi.workflows.show_workflow(workflow_id)['steps'][key][
                'tool_id'] == 'bwa_mem':
                if tissue in \
                        gi.workflows.show_workflow(workflow_id)[u'steps'][key][
                            u'tool_inputs'][u'params'] and tissue != 'NONE':
                    param_set[key] = {'params|readGroup|rgid': sample_id,
                                      'params|readGroup|rgpu':
                                          read_group_dict['rgpu'],
                                      'params|readGroup|rglb':
                                          read_group_dict['rglb'],
                                      'params|readGroup|rgds': description}
                elif tissue == 'NONE':
                    param_set[key] = {'params|readGroup|rgid': sample_id,
                                      'params|readGroup|rgpu':
                                          read_group_dict['rgpu'],
                                      'params|readGroup|rglb':
                                          read_group_dict['rglb'],
                                      'params|readGroup|rgds': description,
                                      'params|readGroup|rgsm': sample_id}

                    # else:
                #     param_set[key] = {'params|readGroup|rgid': sample_id,
                #                       'params|readGroup|rgpu': read_group_dict[
                #                           'rgpu'],
                #                       'params|readGroup|rglb': read_group_dict[
                #                           'rglb'],
                #                       'params|readGroup|rgds': description,
                #                       'params|readGroup|rgsm': tissue}

            elif gi.workflows.show_workflow(workflow_id)['steps'][key][
                'tool_id'] == 'send_to_cgd':
                if tissue in gi.workflows.show_workflow(workflow_id)[
                    u'steps'][key][u'tool_inputs'][u'tissue']:
                        param_set[key] = {'endpoint_choice|runid': runid,
                        'endpoint_choice|barcodeid': barcodeid,
                                  'cgd_client': cgd_client}

            # This needs to be adapted to differentiate between the TUMOR
            # and NORMAL samples in a matched BAM workflow.
            elif gi.workflows.show_workflow(workflow_id)['steps'][key][
                'tool_id'] == 'move_bam':
                if tissue in gi.workflows.show_workflow(workflow_id)[
                    u'steps'][key][u'tool_inputs'][u'sample_type']:
                    param_set[key] = {'run_id': runid, 'sample_id': sample_id,
                                  'history_id': history_id, 'bam_path': bam_path}

            elif gi.workflows.show_workflow(workflow_id)['steps'][key][
                'tool_id'] == 'append_tumor_to_cnv':
                param_set[key] = {'sample_id': sample_id}

    return param_set


def buildDatamap(file_label_dict, gi, workflow_id, library_id, fileid_r1,
                 fileid_r2, fileid_t1_r1, fileid_t2_r2, fileid_n1_r1,
                 fileid_n2_r2, fileid_tsv, fileid_tbam, fileid_nbam):
    """
    Build a datamap structure.  This will supply the necessary connection between 
    workflow inputs and workflow step id values.

    TODO: Refactor all fileid inputs to arrays.  Add option denoting which
    type we are inputting.

    Datamaps have format:
    datamap['<workflow_step_id>'] = {'src':'ldda', 'id':'<library filename id>'}
    """

    datamap = {}
    # Iterate over all workflow inputs
    for key in gi.workflows.show_workflow(workflow_id)['inputs']:
        # Set workflow label, this is the dict value in file_label_dict
        label = gi.workflows.show_workflow(workflow_id)['inputs'][key]['label']
        print(label)
        # If label is R1, just set to the library_id we grabbed after importing the file.
        # This is the same for R2.
        if label == "R1":
            datamap[key] = {'src': 'ldda', 'id': fileid_r1}
        elif label == "R2":
            datamap[key] = {'src': 'ldda', 'id': fileid_r2}
        elif label == "T1":
            datamap[key] = {'src': 'ldda', 'id': fileid_t1_r1}
        elif label == "T2":
            datamap[key] = {'src': 'ldda', 'id': fileid_t2_r2}
        elif label == "N1":
            datamap[key] = {'src': 'ldda', 'id': fileid_n1_r1}
        elif label == "N2":
            datamap[key] = {'src': 'ldda', 'id': fileid_n2_r2}
        elif label == "TSV":
            datamap[key] = {'src': 'ldda', 'id': fileid_tsv}
        elif label == "TBAM":
            datamap[key] = {'src': 'ldda', 'id': fileid_tbam}
        elif label == "NBAM":
            datamap[key] = {'src': 'ldda', 'id': fileid_nbam}
        else:
            # Iterate through the file and label relationship.
            for filename in file_label_dict:
                # If the label corresponding to a particular file is the same as the current label.
                if file_label_dict[filename] == label:
                    # Iterate through all files in data library.
                    for lib in gi.libraries.show_library(library_id, contents=True):
                        # If the file in the data library is the same as the current file from file_label_dict.
                        if lib['name'].split('/')[-1] == filename:
                            # Set the proper datamap format.
                            datamap[key] = {'src': 'ldda', 'id': lib['id']}

    return datamap


def uploadGalaxyDataset(gi, my_library, library_id, filename,
                            file_ext, sleep_time=20):
    """
    Check for the existence of a particular dataset in a Galaxy library before
    uploading.
    """
    if checkDatasetExist(my_library, filename) == False:
        print("Uploading " + filename + " to library " + library_id + ".")
        if file_ext == 'fastq':
            gi.libraries.upload_from_galaxy_filesystem(library_id, filename,
                                                       file_type=file_ext,
                                                       link_data_only='link_to_files')
        elif file_ext == 'bam':
            gi.libraries.upload_from_galaxy_filesystem(library_id, filename,
                                                       file_type=file_ext)
        else:
            raise Exception("Bad file extension passed to "
                            "uploadGalaxyDataset.")

        print("Waiting " + str(sleep_time) +
              " seconds to allow for import to complete.")
        time.sleep(sleep_time)

    else:
        print("Dataset was found in data library.  If filename is not unique, "
              "please delete the file and rerun.")


# Functions that are working with FASTQ or BAM files.
def get_seq_id(filename):
    """
    Find the seqeuncer id field from a FASTQ or BAM.
    :param filename:
    :return: read_groups
    """

    read_groups = []

    if filename.lower().endswith('.gz'):
        cmd = "gzip -cd " + filename + " | head -1"
    elif filename.lower().endswith('.fastq') or \
            filename.lower().endswith('.fq'):
        cmd = "head -1 " + filename
    elif filename.lower().endswith('.bam'):
        cmd = "samtools view " + filename + " | head -1 | cut -f 1"
    else:
        raise Exception("File must end in fastq.gz, fastq, fq, or bam.")

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    while True:
        line = p.stdout.readline()
        if line != '':
            read_groups.extend(line.rstrip('\n').split(':'))
        else:
            break

    print(read_groups)
    return read_groups


def parseFastqId(fastq_id, proj_base):
    """
    According to Broad:
    H0164ALXX140820:2:1101:10003:23460
    H0164ALXX140820:2:1101:15118:25288

    H0164____________ #portion of @RG ID and PU fields indicating Illumina flow cell
    _____ALXX140820__ #portion of @RG PU field indicating barcode or index in a multiplexed run
    _______________:2 #portion of @RG ID and PU fields indicating flow cell lane

    TODO: Also add the other non-required read groups.

    RGCN=String
    CN=String                     Read Group sequencing center name  Default value: null.

    RGDS=String
    DS=String                     Read Group description  Default value: null.

    RGDT=Iso8601Date
    DT=Iso8601Date                Read Group run date  Default value: null.

    RGPI=Integer
    PI=Integer                    Read Group predicted insert size  Default value: null.

    :param fastq_id:
    :return: read_group_dict
    """

    read_group_dict = {}

    read_group_dict['rgid'] = '.'.join([fastq_id[0][:5], fastq_id[1]])
    read_group_dict['rgpu'] = '.'.join([fastq_id[0], fastq_id[1]])
    read_group_dict['rglb'] = proj_base

    return read_group_dict


def matchNameId(json_blob, phrase):
    for entry in json_blob:
        if entry['name'] == phrase:
            return entry['id']


def checkDualIndexed(header):
    """
    Check the SampleSheet before parsing to determine whether it is dual indexed or not.
    """
    if "I5_Index_ID" in header:
        return True

    return False


def parseSampleSheet(sample_sheet, new_sheet, proj_base, CGD_TESTS):
    """
    Iterate through each line of the incoming raw SampleSheet and process valid samples in to a data structure.
    """

    check = False
    sample_info = {}

    sheet1 = sample_sheet
    with sheet1 as sheet:
        for oline in sheet:
            line = oline.rstrip('\n').split(',')

            if check == False and line[0] != "Sample_ID":
                if line[0] == "Experiment Name":
                    new_sheet.write(
                        "Experiment Name," + proj_base + ",,,,,,,,,\n")
                else:
                    new_sheet.write(oline.rstrip('\n') + ',,,\n')

            if line[0] == "Sample_ID":
                check = True
                if checkDualIndexed(oline) == True:
                    new_sheet.write(
                        "Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description,workflow_name,galaxy_id,bam_path\n")
                else:
                    new_sheet.write(oline.rstrip(
                        '\n') + ',workflow_name,galaxy_id,bam_path\n')

            if check == True and line[0] != "Sample_ID":

                sample_id = line[0]
                sample_name = line[1]
                sample_plate = line[2]
                sample_well = line[3]
                i7_index_id = line[4]
                barcode = line[5]

                if len(line) == 8:

                    analysis = line[6]
                    description = line[7]

                elif len(line) == 10:
                    i5_index_id = line[6]
                    i5_barcode = line[7]
                    i7_index_id = i7_index_id + '_' + i5_index_id
                    barcode = barcode + '_' + i5_barcode
                    analysis = line[8]
                    description = line[9]

                # uniq_key = sample_id + "_" + barcode

                if analysis not in sample_info and analysis in CGD_TESTS:
                    sample_info[analysis] = [
                        [sample_id, sample_name, sample_plate, sample_well,
                         i7_index_id, barcode, analysis, description]]
                elif analysis in sample_info:
                    sample_info[analysis].append(
                        [sample_id, sample_name, sample_plate, sample_well,
                         i7_index_id, barcode, analysis, description])
                elif analysis not in sample_info and analysis not in CGD_TESTS:
                    print(
                    analysis + " is not in the list of valid CGD tests. Skipping.")
                else:
                    raise Exception(
                        "There should not be multiple analyses for the same sample.")

    return sample_info


def checkProjDir(proj_dir):
    if not proj_dir.rstrip('\n').endswith('/'):
        return proj_dir + '/'

    return proj_dir


# This needs to die, and be put in to a config file.
def fileToLabel(analysis):
    if analysis == "Clinical Exome" or analysis == "Clinical Exome Validation":
        file_label_dict = {"1000G_omni2.5.hg19.sites_nochr.vcf": "OMNI", \
                           "agilent_cre.interval_list": "PROBES", \
                           "mills.vcf": "MILLS", \
                           "grch37_p13_142.vcf": "DBSNP", \
                           "1000g_phase1_broad_variant.vcf": "1000G", \
                           "db.conf.txt": "LOGIN", \
                           "agilent_probes.bed": "INT_BED", \
                           "ref_GRCh37.p13_top_level.gff3": "GFF", \
                           "richards_genes.txt": "GENES", \
                           "sift_bins.txt": "SIFT", \
                           "agilent_exons_richards.tsv": "EXONS", \
                           "agilent_genes_richards.tsv": "QC_GENES"}

    elif analysis == "Cancer Exome" or analysis == "Broad Cancer Exome":
        file_label_dict = {"agilent_cre.interval_list": "PROBES", \
                           "mills.vcf": "MILLS", \
                           "grch37_p13_146.vcf": "DBSNP", \
                           "1000g_phase1_broad_variant.vcf": "1000G", \
                           "db.conf.txt": "LOGIN", \
                           "ref_GRCh37.p13_top_level.gff3": "GFF", \
                           "all_refseq_genes.txt": "GENES", \
                           "sift_bins.txt": "SIFT", \
                           "tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt": "TX_OR", \
                           "agilent_probes.bed": "INT_BED", \
                           "cosmic_75_coding_hg19_sorted.vcf": "COSMIC", \
                           "centromeres.bed": "CENTROMERE", \
                           "broad_cancer_exome_refseq.txt": "REFSEQ_TX"}

#                           "1000G_omni2.5.hg19.sites_nochr.vcf": "OMNI", \

    elif analysis == "SolidTumor":
        file_label_dict = {"solid_tumor_iad32788.interval_list": "PROBES", \
                           "mills.vcf": "MILLS", \
                           "grch37_p13_146.vcf": "DBSNP", \
                           "1000g_phase1_broad_variant.vcf": "1000G", \
                           "db.conf.txt": "LOGIN", \
                           "solid_tumor_iad32788.genes": "GENES", \
                           "sift_bins.txt": "SIFT", \
                           "tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt": "TX_OR", \
                           "solid_tumor_iad32788.bed": "INT_BED", \
                           "cosmic_75_coding_hg19_sorted.vcf": "COSMIC", \
                           "ref_GRCh37.p13_top_level.gff3": "GFF", \
                           "solidtumorv1_amplicons.txt": "CNV_AMPLICONS", \
                           "solidtumorv1_counts.csv": "CNV_COUNTS", \
                           "solidtumorv1_samples.txt": "CNV_SAMPLES", \
                           "amplihack_ptrl_pon.vcf": "PANEL_OF_NORMALS"}

    elif analysis == "SolidTumorV2":
        file_label_dict = {"solid_tumor_v02_iad86314.interval_list": "PROBES",
                           "mills.vcf": "MILLS",
                           "grch37_p13_146.vcf": "DBSNP",
                           "db.conf.txt": "LOGIN",
                           "solid_tumor_v02_iad86314.genes": "GENES",
                           "sift_bins.txt": "SIFT",
                           "tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt": "TX_OR",
                           "solid_tumor_v02_iad86314.bed": "INT_BED",
                           "cosmic_75_coding_hg19_sorted.vcf": "COSMIC",
                           "ref_GRCh37.p13_top_level.gff3": "GFF",
                           "solidtumorv2_amplicons.txt": "CNV_AMPLICONS",
                           "solidtumorv2_counts.csv": "CNV_COUNTS",
                           "solidtumorv2_samples.txt": "CNV_SAMPLES",
                           "solidtumor_v2_pon.vcf": "PANEL_OF_NORMALS"}

    elif analysis == "SolidTumorV2 TEST":
        file_label_dict = {"solidtumorv2_no_overlap.interval_list": "PROBES",
                           "mills.vcf": "MILLS",
                           "grch37_p13_146.vcf": "DBSNP",
                           "db.conf.txt": "LOGIN",
                           "solid_tumor_v02_iad86314.genes": "GENES",
                           "sift_bins.txt": "SIFT",
                           "tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt": "TX_OR",
                           "solidtumorv2_no_overlap.bed": "INT_BED",
                           "cosmic_75_coding_hg19_sorted.vcf": "COSMIC",
                           "ref_GRCh37.p13_top_level.gff3": "GFF",
                           "solidtumorv2_refseq_transcripts.txt": "REFSEQ_TX",
                           "solidtumorv2_no_overlap_amplicons.txt": "CNV_AMPLICONS",
                           "solidtumorv2_no_overlap_counts.csv": "CNV_COUNTS",
                           "solidtumorv2_no_overlap_samples.txt": "CNV_SAMPLES",
                           "solidtumor_v2_pon.vcf": "PANEL_OF_NORMALS"}

    elif analysis == "SolidTumor TEST":
        file_label_dict = {"solidtumorv1.interval_list": "PROBES", \
                           "mills.vcf": "MILLS", \
                           "grch37_p13_146.vcf": "DBSNP", \
                            "db.conf.txt": "LOGIN", \
                           "solid_tumor_iad32788.genes": "GENES", \
                           "sift_bins.txt": "SIFT", \
                           "tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt": "TX_OR", \
                           "solidtumorv1.bed": "INT_BED", \
                           "cosmic_75_coding_hg19_sorted.vcf": "COSMIC", \
                           "ref_GRCh37.p13_top_level.gff3": "GFF", \
                           "solidtumorv1_refseq.txt": "REFSEQ_TX",
                           "solidtumorv1_no_overlap_amplicons.txt": "CNV_AMPLICONS", \
                           "solidtumorv1_no_overlap_counts.csv": "CNV_COUNTS", \
                           "solidtumorv1_no_overlap_samples.txt": "CNV_SAMPLES", \
                           "amplihack_ptrl_pon.vcf": "PANEL_OF_NORMALS"}

    elif analysis == "CSER":
        file_label_dict = {"db.conf.txt": "LOGIN", \
                           "jhl_genes.txt": "GENES", \
                           "sift_bins.txt": "SIFT"}

    elif analysis == "gene_panels":
        file_label_dict = {"tso.interval_list": "PROBES", \
                           "grch37_p13_146.vcf": "DBSNP", \
                           "db.conf.txt": "LOGIN", \
                           "tso_manifest_genes_richards.txt": "GENES", \
                           "sift_bins.txt": "SIFT", \
                           "tso_exons_richards.tsv": "EXONS", \
                           "tso.bed": "INT_BED", \
                           "ref_GRCh37.p13_top_level.gff3": "GFF", \
                           "mills.vcf": "MILLS", \
                           "tso_manifest_richards.txt": "MANIFEST"}

    elif analysis == "GATK Preprocess":
        file_label_dict = {"grch37_p13_146.vcf": "DBSNP", \
                           "1000g_phase1_broad_variant.vcf": "1000G", \
                           "mills.vcf": "MILLS"}


    elif analysis == "Medical Exome Validation":
        file_label_dict = {"tso.bed": "PROBES",
                           "tso.interval_list": "INTERVALS",
                           "grch37_p13_146.vcf": "DBSNP",
                           "db.conf.txt": "LOGIN",
                           "tso_manifest_genes_richards.txt": "GENES",
                           "sift_bins.txt": "SIFT",
                           "mills.vcf": "MILLS",
                           "ref_GRCh37.p13_top_level.gff3": "GFF"}

    elif analysis == "IDT Validation":
        file_label_dict = {"xgen_idt_probes.interval_list": "PROBES", \
                           "grch37_p13_142.vcf": "DBSNP", \
                           "db.conf.txt": "LOGIN", \
                           "tso_manifest_genes_richards.txt": "GENES", \
                           "sift_bins.txt": "SIFT", \
                           "Homo_sapiens_assembly19.fasta": "REF", \
                           "tso_exons_richards.tsv": "EXONS", \
                           "tso_manifest_richards.txt": "MANIFEST"}

    return file_label_dict


def absolute_to_samba(abs_path, group_name, server="exanas01.ohsu.edu"):
    """
    Take an absolute file path and convert it to an appropriate samba share path.
    """

    # Absolute path looks like:
    # /home/exacloud/clinical/RichardsLab/BAM/d200_bam/DNA-16-00001-1.bam
    # Samba Share Format should be:
    # \\exanas01.ohsu.edu\RichardsLabBAM\d200_bam\DNA-16-00001-1.bam

    abs_path = abs_path.rstrip('\n').split("/")
    print(abs_path)
    my_index = abs_path.index(group_name)
    samba_path = '\\\\' + server + '\\'
    samba_path += abs_path[4] + abs_path[5] + '\\'
    samba_path += '\\'.join(abs_path[6:])
    return samba_path


def updateSampleSheet(curr_hist, samplesheet_old, new_sheet, dest_path,
                      dont_send_bam):
    """
    Include additional information with samplesheet, such as workflow_name and galaxy_id.

    curr_hist = curr_hist_id[hist_name] = [sample, gi.histories.get_histories(name=hist_name)[0]['id']]
    hist_name = proj_base + '_' + sample_proj + '_' + uuid.uuid4().hex[:8]
    samplesheet_old = sample_info[sample_proj] = [[sample_id, sample_name, sample_plate, sample_well, i7_index_id, barcode, analysis, description], [], [], ...]
    uniq_key = sample_id + "_" + barcode
    new_sheet = file handle to write new samplesheet

    Also, update the base path to include a sample BAM.
    """

    for hist in curr_hist:
        ### TODO: Figure out why this needs to be curr_hist[hist][0][0] for CSER workflow.

        # Changed this for CSER workflow, if problems arise, look here.
        sample = curr_hist[hist][0]

        # This is a Samba Share formatted path.
        new_path = dest_path + sample + ".bam"

        for key in samplesheet_old:
            if sample == samplesheet_old[key][0] and "Validation" not in \
                    samplesheet_old[key][6]:
                if dont_send_bam == True:
                    new_sheet.write(','.join(
                        [samplesheet_old[key][0], samplesheet_old[key][1],
                         samplesheet_old[key][2], samplesheet_old[key][3],
                         samplesheet_old[key][4], samplesheet_old[key][5],
                         samplesheet_old[key][6], samplesheet_old[key][7],
                         hist, curr_hist[hist][1], "NA"]))
                else:
                    new_sheet.write(','.join(
                        [samplesheet_old[key][0], samplesheet_old[key][1],
                         samplesheet_old[key][2], samplesheet_old[key][3],
                         samplesheet_old[key][4], samplesheet_old[key][5],
                         samplesheet_old[key][6], samplesheet_old[key][7],
                         hist, curr_hist[hist][1], new_path]))
                new_sheet.write('\n')


# new_sheet.close()


def sendToCgd(infile, runid, cgd_client, endpoint="samplesheet"):
    ### Create command that will send the SampleSheet to the CGD.

    # Make Java version an argument, if we decide we will ever need <java.
    cmd = "/opt/installed/java8/bin/java -jar " + cgd_client + " -f " + infile + " -r " + runid + " -n " + endpoint + " -d"
    p = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

    print(cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except:
        raise Exception("CGD endpoint did not accept " + infile)

    time.sleep(10)  ### This may help prevent error that is occurring.


def checkDatasetExist(ldda_obj, filename):
    ### Check for existence of a dataset in the "automated" data library.

    if filename.split('.')[-1] == "gz":
        for entry in ldda_obj:
            if entry['name'].split('/')[1] == filename.split('/')[-1][:-3]:
                return True
    else:
        for entry in ldda_obj:
            if entry['name'].split('/')[1] == filename.split('/')[-1]:
                return True

    return False


def mergeCancerExomeSampleSheet(sample_sheet):
    """
    INPUT: {Cancer Exome: [[0,1,2,3,4,5,6,7], [], ...]}
           send list of lists
    Take all SampleSheet entries with Cancer Exome Sample Project and merge in to one entry to send to CGD.
    """

    new_sample = []
    if len(sample_sheet) % 3 != 0:
        raise Exception(
            "NextSeq workflows require 3 sample loadings, therefore we would expect 2 tumor and 1 normal sample in the SampleSheet.")
    else:
        for sample in sample_sheet:
            sample_id = sample[0]
            descr = sample[7]
            if descr == "tumor" and new_sample == []:
                new_sample = sample
            elif descr == "tumor":
                new_sample[4] = '_'.join([new_sample[4], sample[4]])
                new_sample[5] = '_'.join([new_sample[5], sample[5]])
            elif descr == "normal":
                new_sample[4] = '_'.join([new_sample[4], sample[4]])
                new_sample[5] = '_'.join([new_sample[5], sample[5]])
                new_sample[
                    7] = 'normal ' + sample_id  ## The Data line should contain both the phrase "normal" and the actual ID of the normal sample
    return new_sample


def checkImportStatus(gi, ldda_obj, filename):
    if filename.split('.')[-1] == "gz":
        for entry in ldda_obj:
            if entry['name'].split('/')[1] == filename.split('/')[-1][:-3]:
                this_obj = gi.datasets.show_dataset(entry['id'], deleted=False,
                                                    hda_ldda='ldda')
                for k in this_obj:
                    if k == "state":
                        if this_obj[k] != "ok":
                            time.sleep(5)
                            return False
    else:
        for entry in ldda_obj:
            if entry['name'].split('/')[1] == filename.split('/')[-1]:
                this_obj = gi.datasets.show_dataset(entry['id'], deleted=False,
                                                    hda_ldda='ldda')
                for k in this_obj:
                    if k == "state":
                        if this_obj[k] != "ok":
                            time.sleep(5)
                            return False


def checkWorkflowNone(workflow_id):
    """
    Check to ensure workflow_id has been filled.  If it hasn't, the most
     likely culprit is that the Galaxy workflow that is being requested has
     not been shared with the person trying to run it.
    """

    if workflow_id == None:
        raise Exception(
            "Please make sure the Galaxy workflow you are running has been "
            "shared with you by the owner.")
    else:
        print("Running workflow: " + workflow_id)


def workflow_parser():

    """
    :return: args
    """

    parser = argparse.ArgumentParser(
        description='Galaxy BioBlend workflow generation script.')

    parser.add_argument('proj_dir', help='Project directory to process.')
    parser.add_argument('cgd_client', help='Path to the cgd_client.')
    parser.add_argument('access_group',
                        help='Access group to which the user whom invoked the '
                             'workflow belongs.  Should be either RichardsLab '
                             'or CorlessLab right now.')
    parser.add_argument('api_key', default='.galaxy_api_key.json',
                        help="Name of API key config file.")
    parser.add_argument('config',
                        default='/home/users/letaw/clinical/installedTest/'
                                'cgd_panels.json',
                        help="Absolute PATH to the config file.")
    parser.add_argument('--bam_dest', help="Location of the directory where "
                                           "BAM files will be served.")
    parser.add_argument('--dont_send_bam_path', action='store_true',
                        help="Don't create a SampleSheet with the bam_path "
                             "heading.")
    parser.add_argument('--send_to_cgd', action='store_true',
                        help="Should a SampleSheet be sent to the CGD?")
    parser.add_argument('--bam_input', action='store_true',
                        help="Input files are BAM, not FASTQ.")

    args = parser.parse_args()

    return args


def set_bam_dest(args, proj_base, *filenames):
    """
    If argument exists, return the path to the samba share.
    If it doesn't, return None.
    :return:
    """

    samba_path = ""
    if args.bam_dest and not args.dont_send_bam_path:
        bam_dest = checkProjDir(args.bam_dest)
        for f in filenames:
            print(f)
            dest_path = bam_dest + proj_base + '/' + f

            if samba_path == "":
                samba_path = absolute_to_samba(dest_path, args.access_group)
            else:
                next_samba = absolute_to_samba(dest_path, args.access_group)
                samba_path = ';'.join([samba_path, next_samba])
        return samba_path
    else:
        print("Not sending BAM to share directory.")
        return None


def retrieve_api_key(filename, username):
    """
    Retrieve an API key associated with a user and return it.
    :return: my_api
    """

    with open(filename, 'r') as api:
        api_key = json.load(api)
        my_api = api_key[username]

    return my_api


def check_file_exist(file_r1, file_r2, proj_dir, sample_id, sample_proj,
                     dataset_import_dir):
    """
    Look for both files, if they don't exist, concatenate them.
    TODO: MD5SUM
    :return:
    """
    if not os.path.isfile(file_r1) and not os.path.isfile(file_r2):
        catFastq(proj_dir, sample_id, sample_proj, dataset_import_dir,
                 file_r1, file_r2)
    else:
        print("File(s) exists.")


def assign_fileid(gi, library_id, file_ext, filename):
    """
Need to rerun the library command so that we can find the uploaded dataset.
    my_library = gi.libraries.show_library(library_id, contents=True)

    :return:
    """
    if file_ext == 'bam':
        while True:
            my_library = gi.libraries.show_library(library_id, contents=True)
            if matchFileId(my_library, filename) != False:
                fileid = matchFileId(my_library, filename)
                return fileid
    elif file_ext == 'fastq':
        while True:
            my_library = gi.libraries.show_library(library_id, contents=True)
            if matchFileId(my_library, os.path.splitext(filename)[0]) != False:
                fileid = matchFileId(my_library, os.path.splitext(filename)[0])
                return fileid

    return None


def matchFileId(my_library, filename):
    for entry in my_library:
#        print(entry['name'].split('/')[1])
#        print(filename.split('/')[-1])
        if entry['name'].split('/')[1] == filename.split('/')[-1]:
            return entry['id']
    return False


def add_samplesheet_fields(samplesheet, specimen, hist_name, hist_id,
                           samba_path):
    """
    These are the fields that will be added to this particular samplesheet
    record.
    :param samplesheet:
    :return:
    """

    samplesheet.data[specimen]['workflow_name'] = hist_name
    samplesheet.sep_count += 1
    samplesheet.data[specimen]['galaxy_id'] = hist_id
    samplesheet.sep_count += 1
    if samba_path != None:
        samplesheet.data[specimen]['bam_path'] = samba_path
        samplesheet.sep_count += 1

    return samplesheet


def check_validation_sample(this_proj, samplesheet, specimen):
    """

    :return:
    """
    if "Validation" in this_proj:
        samplesheet.pop(specimen)

    return samplesheet


def main():

    curr_user = getpass.getuser()  # Find current user.

    # Set arguments for workflow.
    args = workflow_parser()

    # Operations to set the project directory and project basename.
    proj_dir = args.proj_dir  # Set proj_dir to argument.
    print("Project directory is specified on the command line as: " + proj_dir)
    print("Check to make sure trailing slash exists on project directory.")
    proj_dir = checkProjDir(proj_dir)
    proj_base = os.path.basename(os.path.normpath(proj_dir))
    print("Project directory basename is: " + proj_base)

    # Try to set the bam_dest variable based on input, otherwise set to
    # None and don't copy files to BAM folder.  This will cause the tool in
    # Galaxy that does this to fail.

    # Upon installation, a file names 'args.api_key' would need to be
    # created in the user's $HOME directory.
    # Here we will create the PATH to the API and then retrieve the user's
    # API.
    my_home = find_env_var('HOME')
    full_api_key = '/'.join([my_home, args.api_key])
    api_key = retrieve_api_key(full_api_key, curr_user)

    # Open config file.
    with open(args.config, 'r') as panels:
        panel_config = json.load(panels)

    PANELS = panel_config['panels']
    CGD_TESTS = panel_config['cgd_tests'] # Right now CGD_TESTS includes
    # PANELS, which it should not.

    galaxy_url = panel_config['galaxy_url']  ### Galaxy Instance URL
    comp_email_suffix = '@' + panel_config['email_suffix'] + '/'

    dataset_import_dir = (panel_config['dataset_import_dir'] + curr_user +
                          comp_email_suffix)  # All files must go through
                         # this directory.

    gi = galaxy.GalaxyInstance(url=galaxy_url, key=api_key)

    wf_dict = {}

    # Pertaining to the new SampleSheet.  Need to do this with
    # SampleSheetWriter.
    sheet_name = proj_dir + "SampleSheet_cgd_" + proj_base + ".csv"
    #samplesheet_new = open(sheet_name, 'w')

    orig_sheet = proj_dir + "SampleSheet.csv"
    #samplesheet_old = open(orig_sheet, 'rU')
    # proj_base, CGD_TESTS)

    print("Parsing SampleSheet.csv")
    samplesheet_data = samplesheet_class.SampleSheetReader(orig_sheet)
    samplesheet = samplesheet_data.data

    library_id = matchNameId(gi.libraries.get_libraries(),
                             panel_config['galaxy_sequence_data'])
    print("Files will be uploaded to library id: " + library_id)

    richards_lib = matchNameId(gi.libraries.get_libraries(),
                               panel_config['galaxy_resource_files'])

    curr_hist_id = {}  # Collect history id's with history names here,
    # and write them all to new sample sheet.

    my_library = gi.libraries.show_library(library_id, contents=True)

    for specimen in samplesheet:

        try:
            this_proj = samplesheet[specimen]['Sample_Project']
        except:
            raise KeyError("We are looking for a Sample_Project field but "
                           "can't find one.")

        if (this_proj in panel_config['paired_fastq'] or this_proj in
            PANELS):

            # Set Experiment Name to proj_base, which will be the runid in CGD.
            samplesheet_data.header['Experiment Name'][0] = proj_base

            # Match the workflow id's to the particular test.
            if this_proj in PANELS:
                workflow_id = matchNameId(gi.workflows.get_workflows(),
                                          "gene_panels")
            else:
                workflow_id = matchNameId(gi.workflows.get_workflows(),
                                          this_proj)

            print("Running " + this_proj + " analyses.")

            sample_id = samplesheet[specimen]['Sample_ID']
            barcode = samplesheet[specimen]['I7_Index_ID']
            print("Working on sample: " + sample_id)
            checkWorkflowNone(workflow_id)
            print("Concatenating project FASTQ's.")

            # As per EP-270.
            file_r1 = dataset_import_dir + proj_base + '_' + barcode \
                      + '_' + sample_id + '_R1.fastq.gz'
            file_r2 = dataset_import_dir + proj_base + '_' + barcode \
                      + '_' + sample_id + '_R2.fastq.gz'
            print("File 1: " + file_r1)
            print("File 2: " + file_r2)
            check_file_exist(file_r1, file_r2, proj_dir, sample_id,
                             this_proj, dataset_import_dir)

            print("Grabbing id from first entry in input FASTQ.")
            fastq_id = get_seq_id(file_r1)
            print("Your FASTQ ID is: " + ':'.join(fastq_id))

            file_base = sample_id + '.bam'

            samba_path = set_bam_dest(args, proj_base, file_base)

            print("Parsing the FASTQ ID.")
            rg_dict = parseFastqId(fastq_id, proj_base)
            print("The following read groups will be assigned: ")
            print(rg_dict)

            print("Adding " + file_r1 + " and " + file_r2 + " to library " +
            library_id + " if not already present.")
            uploadGalaxyDataset(gi, my_library, library_id, file_r1,
                                'fastq', 1)
            uploadGalaxyDataset(gi, my_library, library_id, file_r2,
                                'fastq', 1)


            print("Finding dataset ID's to match uploaded file names.")
            fileid_r1 = assign_fileid(gi, library_id, 'fastq', file_r1)
            fileid_r2 = assign_fileid(gi, library_id, 'fastq', file_r2)
            print("File ID for R1 is: " + fileid_r1)
            print("File ID for R2 is: " + fileid_r2)

            if this_proj in PANELS:
                file_label_dict = fileToLabel("gene_panels")
            else:
                file_label_dict = fileToLabel(this_proj)

            print("Creating datamap from library to workflow.")
            datamap = buildDatamap(file_label_dict, gi, workflow_id,
                                   richards_lib, fileid_r1, fileid_r2,
                                   None, None, None, None, None, None, None)
            print(datamap)

            print("Invoking Galaxy workflow id: " + workflow_id)

            # Set name to be used in history.
            hist_name = (proj_base + '_' + sample_id + '_' + this_proj +
                         '_' + uuid.uuid4().hex[:8])
            print(hist_name)

            # Create a history that will store the output.
            gi.histories.create_history(name=hist_name)

            print("Finding history id of newly created history.")

            # Find history id associated with current history.
            hist_id = gi.histories.get_histories(name=hist_name)[0]['id']
            curr_hist_id[hist_name] = [sample_id, hist_id]

            print("Creating parameter dictionary.")
            param_set = {}
            # The buildParams function probably needs to be refactored.
            param_set = buildParams(gi, rg_dict, workflow_id, sample_id,
                                    hist_id, proj_base, barcode, "NONE",
                                    param_set, args.bam_dest, args.cgd_client)
            print(param_set)

            wf_dict[sample_id] = [workflow_id, datamap, param_set, hist_id]

            print(samplesheet[specimen])

            samplesheet_data = add_samplesheet_fields(samplesheet_data,
                                                      specimen,
                                                 hist_name, hist_id,
                                                 samba_path)

            print(samplesheet[specimen])

            samplesheet_data.data = check_validation_sample(this_proj,
                                                        samplesheet, specimen)


        # Broad Cancer Exome.
        # Should this also be Cancer Exome?  Otherwise, config files should
        # contain a "fastq_input_matched".
        elif this_proj in panel_config['bam_input_matched']:

            print(samplesheet[specimen])

            sample_id = samplesheet[specimen]['Sample_ID']

            tumor_id = samplesheet[specimen]['Sample_ID']
            normal_id = samplesheet[specimen]['Description'].split(' ')[1]


            # Barcode is hard coded to this string for Broad samples.
            # Set this to an argument instead, as we want to have the bam
            # input generalized for the workflow generator.
            barcode = "BROAD01_BROAD02_BROAD03"
            tumor_barcode = "BROAD01_BROAD02"
            normal_barcode = "BROAD03"

            samplesheet[specimen]['I7_Index_ID'] = barcode
            samplesheet[specimen]['index'] = barcode

            workflow_id = matchNameId(gi.workflows.get_workflows(),
                                      this_proj)

            print("Working on sample: " + tumor_id)
            checkWorkflowNone(workflow_id)
            print("Analysis type: " + this_proj)

            print("Copying project BAM's.")
            file_tbase = tumor_id + '.bam'
            file_nbase = normal_id + '.bam'
            file_tbam = dataset_import_dir + proj_base + '_' + barcode \
                        + '_' + file_tbase
            file_nbam = dataset_import_dir + proj_base + '_' + barcode \
                        + '_' + file_nbase

            samba_path = set_bam_dest(args, proj_base, file_tbase, file_nbase)

            print("Tumor BAM destination: " + file_tbam)
            print("Normal BAM destination: " + file_nbam)


            if not os.path.isfile(file_tbam) and not os.path.isfile(
                    file_nbam):
                curr_tbam = proj_dir + this_proj + '/' + tumor_id + ".bam"
                curr_nbam = proj_dir + this_proj + '/' + normal_id + ".bam"
                shutil.copyfile(curr_tbam, file_tbam)
                shutil.copyfile(curr_nbam, file_nbam)
            else:
                print("File(s) exists.")

            print("Grabbing id from first entry in input FASTQ.")
            tumor_fastq_id = get_seq_id(file_tbam)
            normal_fastq_id = get_seq_id(file_nbam)
            print("Your FASTQ ID's are: " + ':'.join(tumor_fastq_id))
            print("Your FASTQ ID's are: " + ':'.join(normal_fastq_id))

            print("Parsing the FASTQ ID.")
            tumor_rg_dict = parseFastqId(tumor_fastq_id, proj_base)
            normal_rg_dict = parseFastqId(normal_fastq_id, proj_base)
            print("The following read groups will be assigned: ")
            print(tumor_rg_dict)
            print(normal_rg_dict)

            print("Adding " + file_tbam + " and " + file_nbam + " to library "
                  + library_id + " if not already present.")
            uploadGalaxyDataset(gi, my_library, library_id, file_tbam,
                                'bam', 1)
            uploadGalaxyDataset(gi, my_library, library_id, file_nbam,
                                'bam', 1)

            print("Finding dataset ID's to match uploaded file names.")
            fileid_tbam = assign_fileid(gi, library_id, 'bam', file_tbam)
            fileid_nbam = assign_fileid(gi, library_id, 'bam', file_nbam)
            print("File ID for the tumor sample is: " + fileid_tbam)
            print("File ID for the normal sample is: " + fileid_nbam)

            file_label_dict = fileToLabel(this_proj)

            print("Creating datamap from library to workflow.")
            datamap = buildDatamap(file_label_dict, gi, workflow_id,
                                   richards_lib, None, None,
                                   None, None, None, None, None,
                                   fileid_tbam, fileid_nbam)
            print(datamap)

            print("Invoking Galaxy workflow id: " + workflow_id)

            # Set name to be used in history.
            hist_name = (proj_base + '_' + sample_id + '_' + this_proj +
                         '_' + uuid.uuid4().hex[:8])
            print(hist_name)

            # Create a history that will store the output.
            gi.histories.create_history(name=hist_name)

            print("Finding history id of newly created history.")

            # Find history id associated with current history.
            hist_id = gi.histories.get_histories(name=hist_name)[0]['id']
            print(hist_id)

            curr_hist_id[hist_name] = [sample_id, hist_id]


            # Fill parameter set for Galaxy workflow variables.
            print("Creating parameter dictionary.")
            param_set = {}
            # The buildParams function probably needs to be refactored.
            param_set = buildParams(gi, tumor_rg_dict, workflow_id, tumor_id,
                                    hist_id, proj_base, tumor_barcode, "TUMOR",
                                    param_set, args.bam_dest, args.cgd_client)
            param_set = buildParams(gi, normal_rg_dict, workflow_id, normal_id,
                                    hist_id, proj_base, normal_barcode,
                                    "NORMAL", param_set, args.bam_dest,
                                    args.cgd_client)

            wf_dict[sample_id] = [workflow_id, datamap, param_set, hist_id]

            samplesheet_data = add_samplesheet_fields(samplesheet_data,
                                                      specimen, hist_name,
                                                      hist_id, samba_path)


        elif this_proj == "CSER":

            workflow_id = matchNameId(gi.workflows.get_workflows(),
                                      sample_proj)

            for sample in samplesheet[sample_proj]:
                sample_id = sample[0]
                barcode = sample[4]

                print("Working on sample: " + sample_id)
                checkWorkflowNone(workflow_id)
                print("Analysis type: " + sample_proj)

                file_tsv = dataset_import_dir + sample_id + ".tsv"
                src_dir = proj_dir + sample_proj + "/" + sample_id + ".tsv"
                shutil.copy(src_dir, file_tsv)

                os.chmod(file_tsv, 0644)

                gi.libraries.upload_from_galaxy_filesystem(library_id,
                                                           file_tsv,
                                                           file_type='tabular',
                                                           link_data_only='link_to_files')
                ### Add processing of CSER coverage tsv here.
                print("Waiting 30 seconds for files to be available.")
                time.sleep(10)

                my_library = gi.libraries.show_library(library_id,
                                                       contents=True)

                print("Finding dataset ID's to match uploaded file names.")
                fileid_tsv = matchFileId(my_library, library_id, file_tsv)
                print("File ID for TSV is: " + fileid_tsv)

                file_label_dict = fileToLabel(sample_proj)

                print("Creating datamap from library to workflow.")
                datamap = buildDatamap(file_label_dict, gi, workflow_id,
                                       richards_lib, None, None, None, None,
                                       None, None, fileid_tsv)

                hist_name = proj_base + '_' + sample_proj + '_' + uuid.uuid4().hex[
                                                                  :8]

                ### Create a history that will store the output.
                gi.histories.create_history(name=hist_name)

                print("Finding history id of newly created history.")

                ### Find history id associated with current history.
                hist_id = gi.histories.get_histories(name=hist_name)[0]['id']
                curr_hist_id[hist_name] = [sample, hist_id]

                print("Creating parameter dictionary.")
                rg_dict = {}
                param_set = {}
                param_set = buildParams(gi, rg_dict, workflow_id, sample_id,
                                        hist_id, proj_base, barcode, "none",
                                        param_set, bam_dest, args.cgd_client)

                wf_dict[sample_id] = [workflow_id, datamap, param_set, hist_id]
                samplesheet_data[sample_id] = sample



        elif this_proj == "Cancer Exome":

            samplesheet = mergeCancerExomeSampleSheet(
                samplesheet[sample_proj])
            print(samplesheet)

            tumor_id = samplesheet[0]
            normal_id = samplesheet[7].split(' ')[1]
            barcode = samplesheet[4]
            workflow_id = matchNameId(gi.workflows.get_workflows(),
                                      sample_proj)

            print("Working on sample: " + tumor_id)
            checkWorkflowNone(workflow_id)
            print("Analysis type: " + sample_proj)

            print("Concatenating project FASTQ's.")
            file_t1_r1 = dataset_import_dir + tumor_id + "_R1.fastq.gz"
            file_t2_r2 = dataset_import_dir + tumor_id + "_R2.fastq.gz"
            file_n1_r1 = dataset_import_dir + normal_id + "_R1.fastq.gz"
            file_n2_r2 = dataset_import_dir + normal_id + "_R2.fastq.gz"

            print("Tumor File 1: " + file_t1_r1)
            print("Tumor File 2: " + file_t2_r2)
            print("Normal File 1: " + file_n1_r1)
            print("Normal File 2: " + file_n2_r2)

            if not os.path.isfile(file_t1_r1) and not os.path.isfile(
                    file_t2_r2) and not os.path.isfile(
                    file_n1_r1) and not os.path.isfile(file_n2_r2):
                catFastq(proj_dir, tumor_id, sample_proj)
                catFastq(proj_dir, normal_id, sample_proj)
            else:
                print("File(s) exists.")

            print("Grabbing id from first entry in input FASTQ.")
            tumor_fastq_id = getFastqId(file_t1_r1)
            normal_fastq_id = getFastqId(file_n1_r1)
            print("Your FASTQ ID's are: " + ':'.join(tumor_fastq_id))
            print("Your FASTQ ID's are: " + ':'.join(normal_fastq_id))

            print("Parsing the FASTQ ID.")
            tumor_rg_dict = parseFastqId(tumor_fastq_id)
            normal_rg_dict = parseFastqId(normal_fastq_id)
            print("The following read groups will be assigned: ")
            print(tumor_rg_dict)
            print(normal_rg_dict)

            my_library = gi.libraries.show_library(library_id, contents=True)

            uploadGalaxyDataset(gi, my_library, library_id, file_t1_r1)
            uploadGalaxyDataset(gi, my_library, library_id, file_t2_r2)
            uploadGalaxyDataset(gi, my_library, library_id, file_n1_r1)
            uploadGalaxyDataset(gi, my_library, library_id, file_n2_r2)

            print("Finding dataset ID's to match uploaded file names.")
            fileid_t1_r1 = matchFileId(my_library, library_id,
                                       os.path.splitext(file_t1_r1)[0])
            print("File ID for tumor R1 is: " + fileid_t1_r1)
            fileid_t2_r2 = matchFileId(my_library, library_id,
                                       os.path.splitext(file_t2_r2)[0])
            print("File ID for tumor R2 is: " + fileid_t2_r2)
            fileid_n1_r1 = matchFileId(my_library, library_id,
                                       os.path.splitext(file_n1_r1)[0])
            print("File ID for normal R1 is: " + fileid_n1_r1)
            fileid_n2_r2 = matchFileId(my_library, library_id,
                                       os.path.splitext(file_n2_r2)[0])
            print("File ID for normal R2 is: " + fileid_n2_r2)

            file_label_dict = fileToLabel(sample_proj)

            print("Creating datamap from library to workflow.")
            datamap = buildDatamap(file_label_dict, gi, workflow_id,
                                   richards_lib, None, None, fileid_t1_r1,
                                   fileid_t2_r2, fileid_n1_r1, fileid_n2_r2,
                                   None)

            print("Invoking Galaxy workflow id: " + workflow_id)

            ### Set name to be used in history.
            hist_name = proj_base + '_' + sample_proj + '_' + uuid.uuid4().hex[
                                                              :8]

            ### Create a history that will store the output.
            gi.histories.create_history(name=hist_name)

            print("Finding history id of newly created history.")

            ### Find history id associated with current history.
            hist_id = gi.histories.get_histories(name=hist_name)[0]['id']
            curr_hist_id[hist_name] = [samplesheet[0], hist_id]

            print("Creating parameter dictionary.")
            param_set = {}
            param_set = buildParams(gi, tumor_rg_dict, workflow_id, sample_id,
                                    hist_id, proj_base, barcode, "tumor",
                                    param_set, bam_dest, args.cgd_client)
            param_set = buildParams(gi, normal_rg_dict, workflow_id, sample_id,
                                    hist_id, proj_base, barcode, "normal",
                                    param_set, bam_dest, args.cgd_client)
            print(param_set)

            wf_dict[samplesheet[0]] = [workflow_id, datamap, param_set,
                                        hist_id]

            samplesheet_data[samplesheet[0]] = samplesheet

            #    samplesheet_new = open('SampleSheet_new.csv', 'w')

#    if args.dont_send_bam_path:
#        updateSampleSheet(curr_hist_id, samplesheet_data, samplesheet_new,
#                          samba_path, True)
#    else:
#        updateSampleSheet(curr_hist_id, samplesheet_data, samplesheet_new,
#                          samba_path, False)

    # Change permissions on outgoing samplesheet.

#    samplesheet_class.SampleSheetWriter(samplesheet_data, sheet_name)

    try:
        print(wf_dict)
        samplesheet_class.SampleSheetWriter(samplesheet_data, sheet_name)
        os.chmod(sheet_name, 0660)
        for wf in wf_dict.itervalues():
            gi.workflows.run_workflow(wf[0], wf[1], params=wf[2], history_id=wf[3],
                                     import_inputs_to_history=True)
    except:
        os.remove(sheet_name)
        raise

    # Right after updating, send the sheet to CGD.
    if args.send_to_cgd:
        sendToCgd(sheet_name, proj_base, args.cgd_client)
    else:
        print("We will not be sending a SampleSheet to the CGD.")


if __name__ == "__main__":
    main()
