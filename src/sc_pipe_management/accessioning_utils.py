"""
Utility functions for Terra to Portal posting operations.
Functions for parsing, string generation, file handling, and other utilities
that don't directly call .post or .patch operations.
"""

import pandas as pd
import igvf_and_terra_api_tools as api_tools
import datetime
import os
import dataclasses
import json
import re
import subprocess
import firecloud.api as fapi

from constants import (
    ACCESSION_HEADERS_BY_ASSAY_TYPES,
    GENOME_ASSEMBLY_INFO,
    GS_FILE_PATH_REGEX,
    IGVF_URL_PATH_REGEX
)

# Define data classes for structured data handling


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


@dataclasses.dataclass(frozen=True)
class WorkflowConfigInfo:
    doc_aliases: list[str]
    download_path: str


# Functions to parse info from Terra data tables
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
    igvf_file_urls = [url for url in igvf_file_urls if url is not None]
    cleaned_igvf_file_urls = parse_terra_str_list(
        terra_str_lists=igvf_file_urls)
    return [IGVF_URL_PATH_REGEX.search(str(url)).group(1) for url in cleaned_igvf_file_urls if IGVF_URL_PATH_REGEX.search(str(url)) is not None]


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


def get_gs_path_for_terra_output_cols(terra_data_record: pd.Series) -> str:
    """Get GS path from terra data record for submission and workflow IDs parsing, trying RNA first, then ATAC.

    Args:
        terra_data_record (pd.Series): Terra pipeline data record

    Returns:
        str: GS path for workflow UUID parsing

    Raises:
        ValueError: If no valid GS path is found
    """
    # Possible columns to check for GS path (RNA first, then ATAC)
    possible_columns = ['rna_kb_h5ad', 'atac_bam']
    for col in possible_columns:
        # If the column is empty, it returns np.float64(nan)
        if (col in list(terra_data_record.index)) and (pd.notna(terra_data_record[col])):
            return terra_data_record[col]
    raise ValueError(
        'No valid GS path found in the Terra data record for workflow UUID parsing.')


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


# Functions to handle post setups (e.g., file aliases, access levels)
def get_genome_assembly(taxa: str) -> str:
    """Get genome assembly based on taxa.

    Args:
        taxa (str): human or mouse

    Returns:
        str: genome assembly, e.g., GRCh38 or GRCm39
    """
    return GENOME_ASSEMBLY_INFO.get(taxa)  # Default to 'Unknown assembly' if taxa not found


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


def get_existing_analysis_set_docs(analysis_set_acc: str, igvf_utils_api) -> list:
    """Get existing document aliases linked to the analysis set.

    Args:
        analysis_set_acc (str): Analysis set accession
        igvf_utils_api (_type_): igvf utils API client

    Returns:
        list: A list of existing document aliases or an empty list if no documents are linked
    """
    analysis_set_obj = igvf_utils_api.get(f'/analysis-sets/{analysis_set_acc}')
    if not analysis_set_obj.get('documents'):
        return []
    existing_doc_ids = analysis_set_obj['documents']
    all_existing_doc_aliases = set()
    for doc_id in existing_doc_ids:
        doc_obj = igvf_utils_api.get(doc_id)
        if 'aliases' in doc_obj:
            all_existing_doc_aliases.update(doc_obj['aliases'])
    return sorted(all_existing_doc_aliases)


def get_file_aliases(col_header: str, lab: str, terra_data_record: pd.Series, terra_uuids: str) -> list:
    """ Generate a file alias based on the output data name, lab, analysis set accession, and Terra UUIDs.

    Args:
        col_header (str): A data table column header, which is the output file name
        lab (str): The data submitter lab
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        terra_uuids (str): A string of the form GCPbucket_submissionID_workflowID_subworkflowID

    Returns:
        list: A file alias in the form of lab:analysis-set-acc_terra-uuids_col-header_uniform-pipeline
    """
    return [f'{lab.split("/")[-2]}:{terra_data_record["analysis_set_acc"]}_{terra_uuids}_{col_header}_uniform-pipeline']


def mk_doc_aliases(curr_workflow_config: PipelineOutputIds) -> list:
    """Create a list of document aliases for the workflow configuration.

    Args:
        curr_workflow_config (PipelineOutputIds): The current workflow configuration

    Returns:
        list: A list of document aliases
    """
    return [f'igvf-dacc-processing-pipeline:{curr_workflow_config.config_aliases()}_pipeline_config']


def mk_qc_obj_aliases(curr_workflow_config: PipelineOutputIds, analysis_set_acc: str, qc_prefix: str, lab: str) -> list:
    """Create a list of QC objects aliases for the workflow configuration.

    Args:
        curr_workflow_config (PipelineOutputIds): The current workflow configuration
        analysis_set_acc (str): The analysis set accession
        qc_prefix (str): The prefix for the QC object (fragment, gene count, etc.)
        lab (str): The lab name

    Returns:
        list: A list of QC objects aliases
    """
    return [f'{lab.split("/")[-2]}:{analysis_set_acc}_{curr_workflow_config.aliases()}_{qc_prefix}_QC-metrics_uniform-pipeline']


# Download QC files and workflow config JSONs from Terra GCP bucket
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


def download_workflow_config_json(terra_data_record: pd.Series, terra_namespace: str, terra_workspace: str, output_root_dir: str = '/igvf/data/') -> str:
    """Download the workflow configuration JSON file for a given Terra data record.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        terra_namespace (str): Terra project namespace, e.g., DACC_ANVIL
        terra_workspace (str): Terra workspace name, e.g., IGVF Single-Cell Data Processing
        output_root_dir (str, optional): Where the file will be saved to. Defaults to '/igvf/data/'.

    Returns:
        str: path to the saved JSON file
    """
    # NOTE: A bit of hack here to get submission ID (Check RNA data first, then ATAC)
    terra_ids = parse_workflow_uuids_from_gs_path(
        gs_path=get_gs_path_for_terra_output_cols(terra_data_record=terra_data_record))
    curr_workflow_config = api_tools.get_workflow_input_config(terra_namespace=terra_namespace,
                                                               terra_workspace=terra_workspace,
                                                               submission_id=terra_ids.submission_id,
                                                               workflow_id=terra_ids.workflow_id)
    curr_workflow_config['igvf_credentials'] = 'Google cloud path to a txt file with credentials to access the IGVF data portal.'
    json_output_dir = os.path.join(
        output_root_dir, terra_data_record['analysis_set_acc'])
    local_file_path = dump_json(input_json=curr_workflow_config,
                                analysis_set_acc=terra_data_record['analysis_set_acc'],
                                output_root_dir=json_output_dir)
    return WorkflowConfigInfo(doc_aliases=mk_doc_aliases(curr_workflow_config=terra_ids),
                              download_path=local_file_path)


def download_all_workflow_config_jsons(terra_data_table: pd.DataFrame, terra_namespace: str, terra_workspace: str, output_root_dir: str = '/igvf/data/workflow_configs') -> list[str]:
    """Download all workflow configuration JSON files for a given Terra data table.

    Args:
        terra_data_table (pd.DataFrame): The Terra data table containing pipeline runs
        terra_namespace (str): Terra project namespace, e.g., DACC_ANVIL
        terra_workspace (str): Terra workspace name, e.g., IGVF Single-Cell Data Processing
        output_root_dir (str, optional): Where the file will be saved to. Defaults to '/igvf/data/workflow_configs'.

    Returns:
        list[str]: A list of paths to the saved JSON files
    """
    workflow_config_download_paths = {}
    for _, row in terra_data_table.iterrows():
        fapi._set_session()  # Refresh the session for each row
        curr_config_file_info = download_workflow_config_json(terra_data_record=row,
                                                              terra_namespace=terra_namespace,
                                                              terra_workspace=terra_workspace,
                                                              output_root_dir=output_root_dir)
        workflow_config_download_paths[row['analysis_set_acc']
                                       ] = curr_config_file_info
    return workflow_config_download_paths


# Functions to read and write JSON files
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


# Functions to summarize POST status and add it to the output data table
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
        lambda x: x.str.startswith(('POST Failed', 'Conflicted with an existing object')), axis=1).sum(axis=1)
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
