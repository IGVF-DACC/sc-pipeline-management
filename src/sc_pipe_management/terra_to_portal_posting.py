import pandas as pd
import igvf_and_terra_api_tools as api_tools
import logging
import multiprocessing
from functools import partial
import datetime
import os
import dataclasses
import json
import re
from itertools import chain
import subprocess
import requests


# TODO:\
# fapi._session() can time out, how to handle that?
# Figure out how to resume/fix bad posts


# NOTE:
# 1) File aliases are lab:gcp-bucket_submissionID_workflowID_subworkflowID_OutputFileColumn
# 2) RNA doesn't have an assembly (matrix file)


# Hard code lab and award for sandbox and production
POSTING_LAB = '/labs/igvf-dacc-processing-pipeline/'
POSTING_AWARD = '/awards/HG012012/'

# Hard code for Terra output table column types
TERRA_OUTPUT_TABLE_COLUMN_TYPES = {
    # ATACseq data
    'alignment_file': {
        'atac_bam': {
            'file_format': 'bam',
            'description': 'Raw Aligned bam file from Chromap',
            'content_type': 'alignments'
        }
    },
    'index_file': {
        'atac_fragments_index': {
            'file_format': 'tbi',
            'description': 'Raw Fragment file index from Chromap',
            'content_type': 'index'
        },
        'atac_bam_index': {
            'file_format': 'bai',
            'description': 'Raw Aligned bam file index from Chromap',
            'content_type': 'index'
        },
    },
    'tabular_file': {
        'atac_fragments': {
            'file_format': 'bed',
            'description': 'Raw Fragment file from Chromap',
            'content_type': 'fragments',
            'file_format_specifications': ['buenrostro-bernstein:igvf-single-cell-pipeline-fragment-file-specification']
        }
    },
    # RNAseq data
    'matrix_file': {
        'rna_kb_output_folder_tar_gz': {
            'file_format': 'tar',
            'description': 'Raw Tarball containing all the matrices, logs, and bus files generated from kb',
            'content_type': 'kallisto single cell RNAseq output',
            'file_format_specifications': ['buenrostro-bernstein:igvf-sc-pipeline-matrix-tar-specification',
                                           'igvf:igvf-sc-pipeline-rna-tar-mtx-per-file-specification']
        },
        'rna_kb_h5ad': {
            'file_format': 'h5ad',
            'description': 'Raw h5ad containing four separated count matrices: Spliced, Unspliced, Ambiguous, and Total',
            'content_type': 'sparse gene count matrix',
            'file_format_specifications': ['buenrostro-bernstein:igvf-sc-pipeline-matrix-h5-specification']
        },
    },
}

# Hard code for columns with IGVF accessions
ACCESSION_HEADERS_BY_ASSAY_TYPES = {'atac': ['atac_read1_accessions', 'atac_read2_accessions', 'atac_barcode_accessions'],
                                    'rna': ['rna_read1_accessions', 'rna_read2_accessions', 'rna_barcode_accessions']
                                    }

# NOTE: These may end up not being needed if genome tsv reference and assembly column is available in input
# Convert Terra table's info to portal enum
GENOME_ASSEMBLY_INFO = {'Homo sapiens': 'GRCh38', 'Mus musculus': 'GRCm39'}

# Analysis step versions
ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES = {'atac': 'igvf:single-cell-uniform-pipeline-chromap-atacseq-step-v1.1.0',
                                         'rna': 'igvf:single-cell-uniform-pipeline-kalliso-bustools-rnaseq-step-v1.1.0'
                                         }

# QC objects
# Metadata: Table column names
# Attachment: keys are portal properties and values are column IDs
# Metadata: keys are JSON keys, values are portal properties
TERRA_QC_OUTPUTS = {'rnaseq': {'gene_count': {'metadata': ['rna_kb_library_qc_metrics_json', 'rna_kb_run_info_json'],
                                              'description': 'RNAseq Kallisto Bustools QC metric',
                                              'attachment': {'rnaseq_kb_info': 'rna_kb_parameters_json'},
                                              'object_type': 'single_cell_rna_seq_quality_metric',
                                              'metadata_map': {'numRecords': 'n_records',   # starting here is inspect.json
                                                               'numReads': 'n_reads',
                                                               'numBarcodes': 'n_barcodes',
                                                               'medianReadsPerBarcode': 'median_reads_per_barcode',
                                                               'meanReadsPerBarcode': 'mean_reads_per_barcode',
                                                               'numUMIs': 'total_umis',
                                                               'numBarcodeUMIs': 'n_barcode_umis',
                                                               'medianUMIsPerBarcode': 'median_umis_per_barcode',
                                                               'meanUMIsPerBarcode': 'mean_umis_per_barcode',
                                                               'gtRecords': 'gt_records',
                                                               'numBarcodesOnOnlist': 'num_barcodes_on_onlist',
                                                               'percentageBarcodesOnOnlist': 'percentage_barcodes_on_onlist',
                                                               'numReadsOnOnlist': 'num_reads_on_onlist',
                                                               'percentageReadsOnOnlist': 'percentage_reads_on_onlist',
                                                               'n_targets': 'n_targets',    # starting here is run_info.json
                                                               'n_bootstraps': 'n_bootstraps',
                                                               'n_processed': 'n_processed',
                                                               'n_pseudoaligned': 'n_pseudoaligned',
                                                               'n_unique': 'n_unique',
                                                               'p_pseudoaligned': 'p_pseudoaligned',
                                                               'p_unique': 'p_unique',
                                                               'index_version': 'index_version',
                                                               'k-mer length': 'kmer_length'
                                                               }
                                              }
                               },
                    'atacseq': {'alignment': {'metadata': None,
                                              'description': 'ATACseq chromap alignment QC metric',
                                              'attachment': {'atac_bam_summary_stats': 'atac_bam_summary_stats'},
                                              'object_type': 'single_cell_atac_seq_quality_metric'
                                              },
                                'fragment': {'metadata': ['atac_fragments_metrics'],
                                             'description': 'ATACseq chromap fragments QC metric',
                                             'attachment': {'atac_fragment_summary_stats': 'atac_fragments_barcode_summary',
                                                            'atac_fragments_alignment_stats': 'atac_fragments_alignment_stats',
                                                            },
                                             'object_type': 'single_cell_atac_seq_quality_metric',
                                             'metadata_map': {'Number_of_reads': 'n_reads',
                                                              'Number_of_mapped_reads': 'n_mapped_reads',
                                                              'Number_of_uniquely_mapped_reads': 'n_uniquely_mapped_reads',
                                                              'Number_of_reads_have_multi-mappings': 'n_reads_with_multi_mappings',
                                                              'Number_of_candidates': 'n_candidates',
                                                              'Number_of_mappings': 'n_mappings',
                                                              'Number_of_uni-mappings': 'n_uni_mappings',
                                                              'Number_of_multi-mappings': 'n_multi_mappings',
                                                              'Number_of_barcodes_in_whitelist': 'n_barcodes_on_onlist',
                                                              'Number_of_corrected_barcodes': 'n_corrected_barcodes',
                                                              'Number_of_output_mappings_(passed_filters)': 'n_output_mappings',
                                                              'uni-mappings': 'uni_mappings',
                                                              'multi-mappings': 'multi_mappings',
                                                              'total': 'total',
                                                              'percentage_duplicates': 'pct_duplicates',
                                                              }
                                             }
                                }
                    }

# IGVF file url parsing regex for accession
IGVF_URL_PATH_REGEX = re.compile(
    r'\/.*-files\/(IGVF[A-Z0-9]+|TST[A-Z0-9]+)\/@@download')

# GS file path regex
GS_FILE_PATH_REGEX = re.compile(
    r'gs://([a-z0-9\-]+)/submissions/final-outputs/([a-z0-9\-]+)/single_cell_pipeline/([a-z0-9\-]+)/[a-z\-]+/[a-z\_]+/([a-z0-9\-]+)/.*')


# Define how to return errors
@dataclasses.dataclass(frozen=True)
class PostResult:
    col_header: str
    accession: str | None
    error: Exception | None = None

    @classmethod
    def Failure(cls, col_header, error):
        return cls(col_header, None, error)

    @classmethod
    def Success(cls, col_header, accession):
        return cls(col_header, accession, None)

    def Description(self) -> str:
        if self.error:
            return 'POST Failed: ' + str(self.error)
        return self.accession


@dataclasses.dataclass
class RunResult:
    analysis_set_acc: str
    post_results: list[PostResult]
    finish_time: str


@dataclasses.dataclass(frozen=True)
class PipelineOutputIds:
    gcloud_bucket: str
    submission_id: str
    workflow_id: str
    subworkflow_id: str

    def aliases(self) -> str:
        """Generate a string of the form submissionID_workflowID_subworkflowID as file."""
        return f'{self.gcloud_bucket}_{self.submission_id}_{self.workflow_id}_{self.subworkflow_id}'

    def config_aliases(self) -> str:
        """Generate a string of the form submissionID_workflowID."""
        return f'{self.submission_id}_{self.workflow_id}'


def parse_workflow_uuids_from_gs_path(gs_path: str, gs_path_regex: re.Pattern = GS_FILE_PATH_REGEX) -> PipelineOutputIds:
    """Parse workflow UUID from a gs path.

    Args:
        gs_path (str): The gs path to the workflow
        gs_path_regex (re.Pattern, optional): The regex pattern to match the gs path. Defaults to GS_FILE_PATH_REGEX.

    Returns:
        PipelineOutputIds: A dataclass containing the parsed UUIDs.
    """
    matches = re.search(gs_path_regex, gs_path)
    if not matches:
        raise ValueError('Unable to parse workflow UUIDs from file GCP path.')
    uuids = matches.groups()
    return PipelineOutputIds(gcloud_bucket=uuids[0],
                             submission_id=uuids[1],
                             workflow_id=uuids[2],
                             subworkflow_id=uuids[3])


def parse_terra_str_list(terra_str_lists: list[str]) -> list:
    """Parse a string list from Terra to a Python list.

    Args:
        terra_str_list (list[str]): The string list from Terra

    Returns:
        list: A Python list
    """
    parse_list_of_strs = []
    for value in terra_str_lists:
        if str(value) == 'nan':
            continue
        if str(value) == '[]':
            continue
        if str(value) == '[""]':
            continue
        if str(value).startswith("[") and str(value).endswith("]"):
            parse_list_of_strs.extend([ele.strip("'")
                                      for ele in value[1:-1].split(', ')])
        else:
            parse_list_of_strs.append(value)
    return parse_list_of_strs


def parse_igvf_accessions_from_urls(igvf_file_urls: list[str]) -> list[str]:
    """Parse sequence specification URLs to get accessions.

    Args:
        seqspec_download_urls (list): List of seqspec download URLs. Possible to have NaN values.

    Returns:
        list: A list of accessions extracted from the URLs
    """
    cleaned_igvf_file_urls = parse_terra_str_list(
        terra_str_lists=igvf_file_urls)
    return [IGVF_URL_PATH_REGEX.search(str(url)).group(1) for url in cleaned_igvf_file_urls if IGVF_URL_PATH_REGEX.search(str(url)) is not None]


def get_genome_assembly(taxa: str) -> str:
    """Get genome assembly based on taxa.

    Args:
        taxa (str): human or mouse

    Returns:
        str: genome assembly, e.g., GRCh38 or GRCm39
    """
    return GENOME_ASSEMBLY_INFO.get(taxa)  # Default to 'Unknown assembly' if taxa not found


# Helper util functions
def check_single_or_multiome(terra_data_record: pd.Series) -> list:
    """Check based on input measurement set IDs to see if the pipeline is RNAseq-only, ATACseq-only, or multiome-seq

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)

    Returns:
        list: content of the pipeline run
    """
    multiome_map = {'atac_MeaSetIDs': 'ATACseq', 'rna_MeaSetIDs': 'RNAseq'}
    pipeline_output_content = []
    for key, value in multiome_map.items():
        if len(parse_terra_str_list(terra_str_lists=[terra_data_record[key]])) > 0:
            pipeline_output_content.append(value)
    return pipeline_output_content


def get_seqfile_accs_from_table(terra_data_record: pd.Series, seqfile_acc_cols: list) -> list:
    """Get a list of sequence file accessions (for derived_from) from Terra table.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        seqfile_acc_cols (list): columns with sequence file accessions of interest

    Returns:
        list: A list of sequence file accessions
    """
    seqfile_accessions = parse_terra_str_list(
        terra_str_lists=[terra_data_record[seqfile_col] for seqfile_col in seqfile_acc_cols])
    return [seqfile_accession for seqfile_accession in seqfile_accessions if seqfile_accession.startswith(('IGVF', 'TSTF'))]


def get_access_status(sequence_file_accessions: list, igvf_api) -> bool:
    """Get file controlled_access status. If any seq data is controlled access, the output data inherit that.

    Args:
        sequence_file_accessions (list): Input sequence file accessions
        igvf_api (_type_): _description_

    Returns:
        bool: True or False
    """
    access_status = 0
    for seq_file_acc in sequence_file_accessions:
        seqfile_obj = igvf_api.get_by_id(
            f'/sequence-files/{seq_file_acc}').actual_instance
        access_status += seqfile_obj.controlled_access
    if access_status > 0:
        return True
    else:
        return False


def get_seqfile_access_lvl_and_accessions(terra_data_record: pd.Series, assay_type: str, igvf_api) -> tuple:
    """Get sequence files' controlled access info and accession IDs

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        assay_type (str): rna or atac
        igvf_api (_type_): _description_

    Returns:
        tuple: (controlled_access_boolean, sequence file accessions)
    """
    seqfile_accessions = get_seqfile_accs_from_table(
        terra_data_record=terra_data_record, seqfile_acc_cols=ACCESSION_HEADERS_BY_ASSAY_TYPES[assay_type])
    output_data_controlled_access = get_access_status(
        sequence_file_accessions=seqfile_accessions, igvf_api=igvf_api)
    return (output_data_controlled_access, seqfile_accessions)


def get_file_alias(col_header: str, lab: str, terra_data_record: pd.Series, terra_uuids: str) -> str:
    """ Generate a file alias based on the output data name, lab, analysis set accession, and Terra UUIDs.

    Args:
        col_header (str): A data table column header, which is the output file name
        lab (str): The data submitter lab
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        terra_uuids (str): A string of the form GCPbucket_submissionID_workflowID_subworkflowID

    Returns:
        str: A file alias in the form of lab:analysis-set-acc_terra-uuids_col-header_uniform-pipeline
    """
    return f'{lab.split("/")[-2]}:{terra_data_record["analysis_set_acc"]}_{terra_uuids}_{col_header}_uniform-pipeline'


# Helper function to post data object via IGVF utils
def single_post_to_portal(igvf_data_payload: dict, igvf_utils_api, upload_file: bool = False) -> str:
    """POST a single data object to the portal

    Args:
        igvf_data_payload (dict): submission POST payload
        igvf_utils_api (_type_): IGVF utils API client
        upload_file (bool, optional): whether to upload data file. Defaults to False.

    Returns:
        str: accession generated by POST
    """
    _schema_property = igvf_utils_api.get_profile_from_payload(
        igvf_data_payload).properties
    # Post without including long strings
    stdout = igvf_utils_api.post(
        igvf_data_payload, upload_file=upload_file, return_original_status_code=True, truncate_long_strings_in_payload_log=True)
    # NOTE: Some POSTs will not have accessions, so UUID is probably the best to use.
    if 'accession' in stdout[0]:
        return stdout[0]['accession']
    else:
        return stdout[0]['uuid']


def download_qc_file_from_gcp(gs_file_path: str, downloaded_dir: str) -> str:
    """Download Google cloud file to a local directory.

    Args:
        gs_file_path (str): GCP file path, e.g., gs://bucket/submissions/final-outputs/...
        downloaded_dir (str): The local directory to download the file to

    Raises:
        Exception: If the download fails, raise an exception with the error message

    Returns:
        str: Downloaded file path
    """
    if not gs_file_path.startswith('gs://'):
        raise Exception(
            f'File path {gs_file_path} is not a valid GCP file path.')
    if not os.path.exists(downloaded_dir):
        os.makedirs(downloaded_dir)
    file_name = os.path.basename(gs_file_path)
    downloaded_file_path = os.path.join(downloaded_dir, file_name)
    # Download
    gcp_download = subprocess.run(
        ['gcloud', 'storage', 'cp', gs_file_path, downloaded_file_path], capture_output=True)
    if gcp_download.returncode != 0:
        err_msg = ','.join([entry for entry in gcp_download.stderr.decode().split(
            '\n') if entry.startswith('ERROR:')])
        raise Exception(f'GCP download error: {err_msg}.')
    return downloaded_file_path


# Helper function to read JSON file
def read_json_file(json_file_paths: str) -> dict:
    """Read a JSON file and return its content as a dictionary.
    Args:
        json_file_path (str): The path to the JSON file.
    Returns:
        dict: The content of the JSON file as a dictionary.
    """
    json_dict = {}
    for json_file_path in json_file_paths:
        with open(json_file_path, 'r') as file:
            json_dict.update(json.load(file))
    return json_dict


def dump_json(input_json: dict, analysis_set_acc: str, output_root_dir: str = './run_config') -> str:
    """Save a JSON object to a file in the specified output directory.

    Args:
        input_json (dict): The JSON object to save (workflow config json)
        analysis_set_acc (str): Analysis set accession, used to name the output file
        output_root_dir (str, optional): Where the file will be saved to. Defaults to './run_config'.

    Returns:
        str: path to the saved JSON file
    """
    file_dir = os.path.join(
        output_root_dir, datetime.datetime.now().strftime("%m%d%Y"))
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
    output_file_path = os.path.join(
        file_dir, f'{analysis_set_acc}_run_config.json')
    with open(output_file_path, 'w') as file:
        json.dump(input_json, file)
    return output_file_path


def mk_doc_aliases(curr_workflow_config: PipelineOutputIds) -> list:
    """Create a list of document aliases for the workflow configuration.

    Args:
        curr_workflow_config (PipelineOutputIds): The current workflow configuration

    Returns:
        list: A list of document aliases
    """
    return [f'igvf-dacc-processing-pipeline:{curr_workflow_config.config_aliases()}_pipeline_config']


def mk_doc_payload(lab: str, award: str, doc_aliases: list, local_file_path: str) -> dict:
    """Create a document payload for posting to the portal.

    Args:
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        doc_aliases (list): Document aliases
        local_file_path (str): Path to the document attachment

    Returns:
        dict: Document payload ready for posting
    """
    return dict(
        aliases=doc_aliases,
        lab=lab,
        award=award,
        document_type='pipeline parameters',
        description='Terra workflow configuration for the single-cell pipeline run',
        attachment={'path': local_file_path},
        _profile='document'
    )


def mk_anaset_docs_patching_payload(doc_aliases: list, analysis_set_acc: str, igvf_utils_api) -> dict:
    """Create a payload for patching the analysis set with the document accession.

    Args:
        doc_aliases (str): Document aliases
        analysis_set_acc (str): Analysis set accession
        igvf_utils_api (_type_): igvf utils API client

    Returns:
        dict: Payload for patching the analysis set
    """
    return {
        'documents': doc_aliases,
        igvf_utils_api.IGVFID_KEY: f"/analysis-sets/{analysis_set_acc}/",
        '_profile': 'analysis_set'
    }


def post_single_document_to_portal(terra_data_record: pd.Series, lab: str, award: str, terra_namespace: str, terra_workspace: str, igvf_utils_api, output_root_dir: str = '/igvf/data/') -> list[PostResult]:
    """Post workflow configuration document to the portal and link it to the analysis set.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        terra_namespace (str): Terra project namespace, e.g., DACC_ANVIL
        terra_workspace (str): Terra workspace name, e.g., IGVF Single-Cell Data Processing
        igvf_utils_api (_type_): igvf utils API client
        output_root_dir (str, optional): Where the file will be saved to. Defaults to './run_config'.

    Returns:
        list[PostResult]: list(PostResult(col_header='workflow_config_document', accession=IGVF UUID, error=error if any))
    """
    workflow_config_post_results = []
    # All about posting the workflow config document
    try:
        # NOTE: A bit of hack here to get submission ID (RNA data will always be available)
        terra_ids = parse_workflow_uuids_from_gs_path(
            gs_path=terra_data_record['rna_kb_h5ad'])
        curr_workflow_config = api_tools.get_workflow_input_config(
            terra_namespace=terra_namespace, terra_workspace=terra_workspace, submission_id=terra_ids.submission_id, workflow_id=terra_ids.workflow_id)
        # Remove credential file path
        curr_workflow_config['igvf_credentials'] = 'Google cloud path to a txt file with credentials to access the IGVF data portal.'
        # Save config to a local file
        json_output_dir = os.path.join(
            output_root_dir, 'workflow_configs', terra_data_record['analysis_set_acc'])
        local_file_path = dump_json(
            input_json=curr_workflow_config, analysis_set_acc=terra_data_record['analysis_set_acc'], output_root_dir=json_output_dir)
        # Make a document object payload
        doc_aliases = mk_doc_aliases(curr_workflow_config=terra_ids)
        doc_payload = mk_doc_payload(
            lab=lab,
            award=award,
            doc_aliases=doc_aliases,
            local_file_path=local_file_path
        )
        # Post the document to the portal
        curr_post_res = single_post_to_portal(igvf_data_payload=doc_payload,
                                              igvf_utils_api=igvf_utils_api,
                                              upload_file=False)
        workflow_config_post_results.append(PostResult.Success(
            col_header='workflow_config_document', accession=curr_post_res))
        # Patch the analysis set with the document accession
        try:
            igvf_utils_api.patch(mk_anaset_docs_patching_payload(
                doc_aliases=doc_aliases,
                analysis_set_acc=terra_data_record['analysis_set_acc'],
                igvf_utils_api=igvf_utils_api
            ))
            workflow_config_post_results.append(PostResult.Success(
                col_header='analysis_set_workflow_config_patch', accession=terra_data_record["analysis_set_acc"]))
        # Log error of patching analysis set with workflow config document
        except requests.exceptions.HTTPError as e:
            workflow_config_post_results.append(PostResult.Failure(
                col_header='analysis_set_workflow_config_patch', error=e))
    # Log error of posting workflow config document
    except requests.exceptions.HTTPError as e:
        workflow_config_post_results.append(PostResult.Failure(
            col_header='workflow_config_document', error=e))
    # Return the post results
    return workflow_config_post_results


def post_single_qc_metric_to_portal(terra_data_record: pd.Series, qc_data_info: dict, qc_of_file_accs: list, lab: str, award: str, analysis_step_version: str, qc_prefix: str, igvf_utils_api, output_root_dir: str = '/igvf/data/') -> list[PostResult]:
    """Post a single QC metric to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        qc_data_info (dict): The QC data collection
        qc_of_file_accs (list): The QC file accessions
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        analysis_step_version (str): The analysis step version
        qc_prefix (str): The prefix for the QC object type (fragment, gene count, etc.)
        igvf_utils_api (_type_): IGVF python client api
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.

    Returns:
        PostResult: list(PostResult(col_header=qc_obj_type, accession=..., error=...))
    """
    try:
        # make some base QC obj with necessary info
        qc_aliases = [
            f'{lab.split("/")[-2]}:{terra_data_record["analysis_set_acc"]}_{qc_data_info["object_type"]}_{qc_prefix}_uniform-pipeline']
        # NOTE: GCP VM persistent disk is /igvf/data
        curr_qc_file_dir = os.path.join(
            output_root_dir, 'qc_metrics', terra_data_record['analysis_set_acc'])
        if not os.path.exists(curr_qc_file_dir):
            os.makedirs(curr_qc_file_dir)
        # Base QC object payload
        qc_payload = dict(lab=lab,
                          award=award,
                          aliases=qc_aliases,
                          analysis_step_version=analysis_step_version,
                          description=qc_data_info['description'],
                          quality_metric_of=qc_of_file_accs,
                          _profile=qc_data_info['object_type']
                          )
        # Add attachments (key is portal property name, value is table colname)
        # Do this first because some QC obj only has attachments
        if qc_data_info['attachment'] is not None:
            for key, value in qc_data_info['attachment'].items():
                # Download the file first
                curr_attachment_file = download_qc_file_from_gcp(
                    gs_file_path=terra_data_record[value], downloaded_dir=curr_qc_file_dir)
                qc_payload.update({key: {'path': curr_attachment_file}})
        # Add on other metadata
        if qc_data_info['metadata'] is not None:
            # Downloadec JSON file paths
            curr_metadata_files = [download_qc_file_from_gcp(
                gs_file_path=terra_data_record[metadata_col], downloaded_dir=curr_qc_file_dir) for metadata_col in qc_data_info['metadata']]
            # Parse Json files
            qc_json_parsed = read_json_file(
                json_file_paths=curr_metadata_files)
            # Update base qc payload
            for key, value in qc_data_info['metadata_map'].items():
                qc_payload.update({value: qc_json_parsed[key]})
        curr_post_res = single_post_to_portal(igvf_data_payload=qc_payload,
                                              igvf_utils_api=igvf_utils_api,
                                              upload_file=False
                                              )
        # Post result
        return PostResult.Success(col_header=f'{qc_prefix}:{qc_data_info["object_type"]}', accession=curr_post_res)
    except requests.exceptions.HTTPError as e:
        return PostResult.Failure(col_header=f'{qc_prefix}:{qc_data_info["object_type"]}', error=e)


# Posting RNA data to the portal
def post_single_matrix_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_derived_from: list, curr_description: str, lab: str, award: str, curr_content_type: str, curr_file_specification: list, igvf_utils_api, upload_file: bool) -> PostResult:
    """Post a single tabular file (for RNA output) to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_assembly (str): The genome assembly
        curr_content_type (str): The content type
        curr_file_specification (list): The file format specification
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        PostResult: PostResult(col_header=col_header, accession=..., error=...)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_alias = get_file_alias(col_header=col_header,
                                         lab=lab,
                                         terra_data_record=terra_data_record,
                                         terra_uuids=parse_workflow_uuids_from_gs_path(
                                             gs_path=curr_gs_cloud_link).aliases()
                                         )
        # Calculate md5sum
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        # Only one reference file
        reference_files = parse_igvf_accessions_from_urls(
            igvf_file_urls=[terra_data_record['kb_index']])
        curr_mtx_file_payload = dict(award=award,
                                     lab=lab,
                                     analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES['rna'],
                                     aliases=[curr_file_alias],
                                     content_type=curr_content_type,
                                     md5sum=curr_md5sum,
                                     file_format=curr_file_format,
                                     derived_from=curr_derived_from,
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     principal_dimension='cell',
                                     secondary_dimensions=['gene'],
                                     filtered=False,
                                     reference_files=reference_files,
                                     description=curr_description,
                                     file_format_specifications=curr_file_specification,
                                     _profile='matrix_file')
        # Post the metadata
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_mtx_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        # Upload files
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except requests.exceptions.HTTPError as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_all_rna_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/') -> list[PostResult]:
    """POST all RNA data as matrix files to the portal (H5AD, tarball, and QC metrics).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        igvf_utils_api (_type_): _description_
        upload_file (bool): Whether to upload the file to the portal
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.

    Returns:
        list[PostResult]: list[PostResult], e.g., [PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    # Get sequence file accessions
    curr_seqfile_accs = get_seqfile_accs_from_table(
        terra_data_record=terra_data_record, seqfile_acc_cols=ACCESSION_HEADERS_BY_ASSAY_TYPES['rna'])

    # Get seqspec and barcode replacement file accessions (barcode file can be NaN on Terra)
    curr_seqspec_and_barcode_replacement = parse_igvf_accessions_from_urls(
        igvf_file_urls=[terra_data_record['rna_seqspec_urls'], terra_data_record['barcode_replacement_file']])

    # Combine seqfile and seqspec into one list as derived_from
    all_derived_from = list(set(
        chain(curr_seqfile_accs, curr_seqspec_and_barcode_replacement)))

    # Prepare arguments for multiprocessing
    tasks = [
        (
            terra_data_record,
            col_header,
            file_info['file_format'],
            all_derived_from,
            file_info['description'],
            lab,
            award,
            file_info['content_type'],
            file_info['file_format_specifications'],
            igvf_utils_api,
            upload_file
        )
        for col_header, file_info in TERRA_OUTPUT_TABLE_COLUMN_TYPES['matrix_file'].items()
    ]

    # Use multiprocessing to post RNA matrix files
    with multiprocessing.Pool() as pool:
        results = pool.starmap(post_single_matrix_file, tasks)

    # Post RNA QC metrics
    posted_rna_accessions = [
        res.accession for res in results if res.error is None]
    # NOTE: Has some hard coded values for now
    qc_metric_post_result = post_single_qc_metric_to_portal(terra_data_record=terra_data_record,
                                                            qc_data_info=TERRA_QC_OUTPUTS['rnaseq']['gene_count'],
                                                            qc_of_file_accs=posted_rna_accessions,
                                                            lab=lab,
                                                            award=award,
                                                            analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                                                'rna'],
                                                            qc_prefix='gene_count',
                                                            igvf_utils_api=igvf_utils_api,
                                                            output_root_dir=output_root_dir
                                                            )
    # Add QC post report
    results.append(qc_metric_post_result)
    return results


# Post ATAC data to the portal
def post_single_alignment_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_assembly: str, curr_ctrl_access: bool, curr_derived_from: list, igvf_utils_api, upload_file: bool) -> PostResult:
    """Post one alignment file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_assembly (str): The genome assembly
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        PostResult: e.g., PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_alias = get_file_alias(col_header=col_header,
                                         lab=lab,
                                         terra_data_record=terra_data_record,
                                         terra_uuids=parse_workflow_uuids_from_gs_path(
                                             gs_path=curr_gs_cloud_link).aliases()
                                         )
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        # Get reference files
        reference_files = parse_igvf_accessions_from_urls(
            igvf_file_urls=[terra_data_record['chromap_index'], terra_data_record['genome_fasta']])
        # Make payload
        curr_alignment_payload = dict(award=award,
                                      lab=lab,
                                      analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES['atac'],
                                      aliases=[curr_file_alias],
                                      assembly=curr_assembly,
                                      controlled_access=curr_ctrl_access,
                                      file_format=curr_file_format,
                                      content_type='alignments',
                                      filtered=False,
                                      derived_from=curr_derived_from,
                                      file_set=terra_data_record['analysis_set_acc'],
                                      md5sum=curr_md5sum,
                                      submitted_file_name=curr_gs_cloud_link,
                                      reference_files=reference_files,
                                      redacted=False,
                                      description=curr_description,
                                      _profile='alignment_file'
                                      )
        # Post
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_alignment_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except requests.exceptions.HTTPError as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_single_tabular_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_assembly: str, curr_ctrl_access: bool, curr_derived_from: list, curr_file_specification: list, igvf_utils_api, upload_file: bool) -> PostResult:
    """Post one fragment tabular file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_assembly (str): The genome assembly
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The sequence file accessions (for derived_from)
        curr_file_specification (list): The file format specification
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        PostResult: e.g., PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_alias = get_file_alias(col_header=col_header,
                                         lab=lab,
                                         terra_data_record=terra_data_record,
                                         terra_uuids=parse_workflow_uuids_from_gs_path(
                                             gs_path=curr_gs_cloud_link).aliases()
                                         )
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        curr_tab_file_payload = dict(award=award,
                                     lab=lab,
                                     analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES['atac'],
                                     aliases=[curr_file_alias],
                                     content_type='fragments',
                                     controlled_access=curr_ctrl_access,
                                     file_format=curr_file_format,
                                     md5sum=curr_md5sum,
                                     derived_from=curr_derived_from,
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     description=curr_description,
                                     filtered=False,
                                     assembly=curr_assembly,
                                     file_format_specifications=curr_file_specification,
                                     _profile='tabular_file'
                                     )
        # Work on schema dependency
        if curr_file_format == 'bed':
            curr_tab_file_payload['file_format_type'] = 'bed3+'
        # Post
        logging.info(f'Posting {curr_file_alias} to the portal.')
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_tab_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except requests.exceptions.HTTPError as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_single_index_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_ctrl_access: bool, curr_derived_from: list, igvf_utils_api, upload_file: bool) -> PostResult:
    """Post one index file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The file that this index file is index for (for derived_from)
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        PostResult: e.g., PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_alias = get_file_alias(col_header=col_header,
                                         lab=lab,
                                         terra_data_record=terra_data_record,
                                         terra_uuids=parse_workflow_uuids_from_gs_path(
                                             gs_path=curr_gs_cloud_link).aliases()
                                         )
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        curr_index_file_payload = dict(award=award,
                                       lab=lab,
                                       analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                           'atac'],
                                       aliases=[curr_file_alias],
                                       content_type='index',
                                       controlled_access=curr_ctrl_access,
                                       file_format=curr_file_format,
                                       md5sum=curr_md5sum,
                                       derived_from=list(
                                           set(curr_derived_from)),
                                       submitted_file_name=curr_gs_cloud_link,
                                       file_set=terra_data_record['analysis_set_acc'],
                                       description=curr_description,
                                       _profile='index_file'
                                       )
        logging.info(f'Posting {curr_file_alias} to the portal.')
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_index_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except requests.exceptions.HTTPError as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_all_atac_alignment_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, curr_assembly: str, curr_ctrl_access: bool, curr_derived_from: list, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/') -> list[PostResult]:
    """Post all ATAC alignment data to the portal (alignment, bam index files, and QC).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        curr_assembly (str): The genome assembly
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        award (str): The data submitter lab's award
        igvf_utils_api (_type_): IGVF utils api
        upload_file (bool): Whether to upload the file to the portal
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.

    Returns:
        list[PostResult]: [PostResult(
            col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    curr_post_summary = []
    # Post alignment bam
    alignment_file_col = 'atac_bam'
    alignment_result = post_single_alignment_file(terra_data_record=terra_data_record,
                                                  col_header=alignment_file_col,
                                                  curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                      'alignment_file'][alignment_file_col]['file_format'],
                                                  curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                      'alignment_file'][alignment_file_col]['description'],
                                                  lab=lab,
                                                  award=award,
                                                  curr_assembly=curr_assembly,
                                                  curr_ctrl_access=curr_ctrl_access,
                                                  curr_derived_from=curr_derived_from,
                                                  igvf_utils_api=igvf_utils_api,
                                                  upload_file=upload_file)
    curr_post_summary.append(alignment_result)

    # Post bam index .bam.bai files
    bam_index_col = 'atac_bam_index'
    bam_index_result = post_single_index_file(terra_data_record=terra_data_record,
                                              col_header=bam_index_col,
                                              curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                  'index_file'][bam_index_col]['file_format'],
                                              curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                  'index_file'][bam_index_col]['description'],
                                              lab=lab,
                                              award=award,
                                              curr_ctrl_access=curr_ctrl_access,
                                              curr_derived_from=[
                                                  res.Description() for res in curr_post_summary],
                                              igvf_utils_api=igvf_utils_api,
                                              upload_file=upload_file)
    curr_post_summary.append(bam_index_result)

    # Submit Alignment QC
    # NOTE: Has some hard coded values for now
    posted_alignment_accs = [
        res.accession for res in curr_post_summary if res.error is None]
    alignment_metric_post_result = post_single_qc_metric_to_portal(terra_data_record=terra_data_record,
                                                                   qc_data_info=TERRA_QC_OUTPUTS['atacseq']['alignment'],
                                                                   qc_of_file_accs=posted_alignment_accs,
                                                                   lab=lab,
                                                                   award=award,
                                                                   analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                                                       'atac'],
                                                                   qc_prefix='alignment',
                                                                   igvf_utils_api=igvf_utils_api,
                                                                   output_root_dir=output_root_dir
                                                                   )
    # Add QC post report
    curr_post_summary.append(alignment_metric_post_result)
    # Return all post res
    return curr_post_summary


def post_all_atac_fragment_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, curr_assembly: str, curr_derived_from: list, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/') -> list[PostResult]:
    """Post all ATAC fragment data to the portal (fragment files, index files, and QC).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The controlled access status (computed based on SeqFile status)
        curr_assembly (str): The genome assembly
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        curr_file_specification (list): The file format specifications
        igvf_api (_type_): IGVF python client api
        igvf_utils_api (_type_): IGVF utils api
        upload_file (bool): Whether to upload a file or not
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.

    Returns:
        list[PostResult]: [PostResult(
            col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    curr_post_summary = []
    # Post fragment file
    fragment_file_col = 'atac_fragments'
    fragment_file_result = post_single_tabular_file(terra_data_record=terra_data_record,
                                                    col_header=fragment_file_col,
                                                    curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                        'tabular_file'][fragment_file_col]['file_format'],
                                                    curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                        'tabular_file'][fragment_file_col]['description'],
                                                    lab=lab,
                                                    award=award,
                                                    curr_assembly=curr_assembly,
                                                    curr_ctrl_access=False,
                                                    curr_derived_from=curr_derived_from,
                                                    curr_file_specification=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                        'tabular_file'][fragment_file_col]['file_format_specifications'],
                                                    igvf_utils_api=igvf_utils_api,
                                                    upload_file=upload_file)
    curr_post_summary.append(fragment_file_result)

    # Post fragment index file
    fragment_index_file_col = 'atac_fragments_index'
    fragment_index_file_result = post_single_index_file(terra_data_record=terra_data_record,
                                                        col_header=fragment_index_file_col,
                                                        curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                            'index_file'][fragment_index_file_col]['file_format'],
                                                        curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                            'index_file'][fragment_index_file_col]['description'],
                                                        lab=lab,
                                                        award=award,
                                                        curr_ctrl_access=False,
                                                        curr_derived_from=[
                                                            res.Description() for res in curr_post_summary],
                                                        igvf_utils_api=igvf_utils_api,
                                                        upload_file=upload_file
                                                        )
    curr_post_summary.append(fragment_index_file_result)
    # Post Fragment QC
    posted_fragment_accs = [
        res.accession for res in curr_post_summary if res.error is None]
    fragment_metric_post_result = post_single_qc_metric_to_portal(terra_data_record=terra_data_record,
                                                                  qc_data_info=TERRA_QC_OUTPUTS['atacseq']['fragment'],
                                                                  qc_of_file_accs=posted_fragment_accs,
                                                                  lab=lab,
                                                                  award=award,
                                                                  analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                                                      'atac'],
                                                                  qc_prefix='fragment',
                                                                  igvf_utils_api=igvf_utils_api,
                                                                  output_root_dir=output_root_dir
                                                                  )
    # Add QC post report
    curr_post_summary.append(fragment_metric_post_result)
    # Return results
    return curr_post_summary


def execute_partial_function(func):
    """Helper function to execute a partial function."""
    return func()


def post_all_atac_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_api, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/') -> list[PostResult]:
    """Post all ATAC data to the portal (alignment, fragment, and index files).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.

    Returns:
        list[PostResult]: list[PostResult], e.g., [PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    # Current assembly
    curr_assembly = get_genome_assembly(taxa=terra_data_record['taxa'])

    # Get file controlled access level and sequence files
    curr_ctrl_access, curr_seqfile_accs = get_seqfile_access_lvl_and_accessions(
        terra_data_record=terra_data_record, assay_type='atac', igvf_api=igvf_api)

    # Get seqspec file accessions
    curr_seqspec_accs = parse_igvf_accessions_from_urls(
        igvf_file_urls=[terra_data_record['atac_seqspec_urls']])

    # Derived from will be a combination of seqfile, seqspec accessions, and barcode replacement file
    # Beware that some assays will have ATAC barcodes and Read2 concatenated into one fastq
    all_derived_from = list(set(chain(curr_seqfile_accs, curr_seqspec_accs)))

    # Multiprocessing tasks for posting ATAC data
    tasks = [
        partial(post_all_atac_alignment_data_to_portal,
                terra_data_record=terra_data_record,
                lab=lab,
                award=award,
                curr_assembly=curr_assembly,
                curr_ctrl_access=curr_ctrl_access,
                curr_derived_from=all_derived_from,
                igvf_utils_api=igvf_utils_api,
                upload_file=upload_file,
                output_root_dir=output_root_dir),
        partial(post_all_atac_fragment_data_to_portal,
                terra_data_record=terra_data_record,
                lab=lab,
                award=award,
                curr_assembly=curr_assembly,
                curr_derived_from=all_derived_from,
                igvf_utils_api=igvf_utils_api,
                upload_file=upload_file,
                output_root_dir=output_root_dir)
    ]

    # Use multiprocessing to execute tasks
    atac_post_results = []
    with multiprocessing.Pool() as pool:
        results = pool.map(execute_partial_function, tasks)
        for result in results:
            atac_post_results.extend(result)
    return atac_post_results


# Full data posting from a single pipeline run
def post_all_data_from_one_run(terra_data_record: pd.Series, igvf_api, igvf_utils_api, upload_file: bool, terra_namespace: str, terra_workspace: str, pipeline_data_lab: str = POSTING_LAB, pipeline_data_award: str = POSTING_AWARD, output_root_dir: str = '/igvf/data/') -> RunResult:
    """Post all single cell uniform pipeline output to the portal from Terra.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        pipeline_data_lab (str): The data submitter lab, default IGVF pipeline processing lab
        pipeline_data_award (str): The data submitter lab's award, default IGVF DACC award
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.

    Returns:
        RunResult: RunResult(analysis_set_acc=IGVFxxx, post_results=[PostResult(
            col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)], finish_time=timestamp)
    """
    logging.info(
        f'Posting data from {terra_data_record["analysis_set_acc"]} to the portal.')
    # Check if the pipeline is RNAseq-only, ATACseq-only, or multiome-seq
    pipeline_run_output = check_single_or_multiome(
        terra_data_record=terra_data_record)
    # Posting
    post_results = []
    # If RNAseq results are present, post them
    if 'RNAseq' in pipeline_run_output:
        post_results += post_all_rna_data_to_portal(terra_data_record=terra_data_record,
                                                    lab=pipeline_data_lab,
                                                    award=pipeline_data_award,
                                                    igvf_utils_api=igvf_utils_api,
                                                    upload_file=upload_file,
                                                    output_root_dir=output_root_dir
                                                    )
    # If ATACseq results are present, post them
    if 'ATACseq' in pipeline_run_output:
        post_results += post_all_atac_data_to_portal(terra_data_record=terra_data_record,
                                                     lab=pipeline_data_lab,
                                                     award=pipeline_data_award,
                                                     igvf_api=igvf_api,
                                                     igvf_utils_api=igvf_utils_api,
                                                     upload_file=upload_file,
                                                     output_root_dir=output_root_dir
                                                     )

    # Add posting run config document
    post_results += post_single_document_to_portal(terra_data_record=terra_data_record,
                                                   lab=pipeline_data_lab,
                                                   award=pipeline_data_award,
                                                   terra_namespace=terra_namespace,
                                                   terra_workspace=terra_workspace,
                                                   igvf_utils_api=igvf_utils_api,
                                                   output_root_dir=output_root_dir)

    # Return full post res
    return RunResult(
        analysis_set_acc=terra_data_record['analysis_set_acc'],
        post_results=post_results,
        finish_time=str(datetime.datetime.now()))


# Post the entire batch of pipeline outputs to the portal
def post_all_successful_runs(full_terra_data_table: pd.DataFrame, igvf_api, igvf_utils_api, upload_file: bool, terra_namespace: str, terra_workspace: str, output_root_dir: str = '/igvf/data/') -> list[PostResult]:
    """Post all successful runs from a Terra data table to the portal.

    Args:
        full_terra_data_table (pd.DataFrame): The Complete Terra data table
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.

    Returns:
        list[PostResult]: PostResult objects
    """
    results = []
    for _, curr_pipeline in full_terra_data_table.iterrows():
        results.append(post_all_data_from_one_run(
            terra_data_record=curr_pipeline, igvf_api=igvf_api, igvf_utils_api=igvf_utils_api, upload_file=upload_file, terra_namespace=terra_namespace, terra_workspace=terra_workspace, output_root_dir=output_root_dir))
    return results


def summarize_post_status(post_results: list) -> pd.DataFrame:
    """Summarize the POST status of all data from a single pipeline run.

    Args:
        post_status_df (pd.DataFrame): The POST status DataFrame

    Returns:
        pd.DataFrame: The POST status summary DataFrame
    """
    post_status_summary = {}
    for result in post_results:
        accession_results = post_status_summary.setdefault(
            result.analysis_set_acc, {'time_stamp': result.finish_time})
        for post_res in result.post_results:
            accession_results[post_res.col_header] = post_res.Description()
    post_status_summary_table = pd.DataFrame(
        post_status_summary).transpose().fillna('')
    post_status_summary_table['post_status_fail'] = post_status_summary_table.apply(
        lambda x: x.str.startswith(('POST Failed')), axis=1).sum(axis=1)
    return post_status_summary_table


def add_post_status_summary_to_output_data_table(full_terra_data_table: pd.DataFrame, post_status_df: pd.DataFrame) -> pd.DataFrame:
    """Add the POST status summary to the output data table.

    Args:
        full_terra_data_table (pd.DataFrame): The Complete Terra data table
        post_status_df (pd.DataFrame): The POST status DataFrame

    Returns:
        pd.DataFrame: The output data table with POST status summary
    """
    return pd.merge(left=full_terra_data_table, right=post_status_df[['time_stamp', 'post_status_fail']], left_on='analysis_set_acc', right_index=True)


def save_pipeline_postres_tables(pipeline_postres_table: pd.DataFrame, updated_full_data_table: pd.DataFrame, output_root_dir: str):
    """Save the pipeline input table to a TSV file.

    Args:
        pipeline_input_table (pd.DataFrame): The pipeline input table
        output_dir (str): The output directory
    """
    curr_datetime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    # make output dir if it does not exist
    if os.path.exists(output_root_dir) is False:
        os.makedirs(output_root_dir)

    # The original table has an extra index column from merging
    updated_full_data_table.to_csv(os.path.join(
        output_root_dir, f'full_terra_data_table_with_post_status_{curr_datetime}.tsv'), sep='\t', index=False)

    # Change post result detail table index to match that of updated full data table for easy Terra upload
    postres_table_index_name = [
        col for col in updated_full_data_table.columns if col.startswith('entity:')][0]
    pipeline_postres_table.index.name = postres_table_index_name.replace(
        '_id', '_postres_id')
    # Output post result detail table in Terra format
    pipeline_postres_table.to_csv(os.path.join(
        output_root_dir, f'single-cell_uniform_pipeline_postres_detail_{curr_datetime}.tsv'), sep='\t')
