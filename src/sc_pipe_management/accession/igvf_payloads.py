import os
import datetime
import subprocess
import pandas as pd
import firecloud.api as fapi

from sc_pipe_management.accession.parse_terra_metadata import (
    WorkflowConfigInfo,
    TerraJobUUIDs,
    TerraOutputMetadata
)

from sc_pipe_management.igvf_and_terra_api_tools import (
    get_workflow_input_config
)

# NOTE: These may end up not being needed if genome tsv reference and assembly column is available in input
# Convert Terra table's info to portal enum
GENOME_ASSEMBLY_INFO = {'Homo sapiens': 'GRCh38', 'Mus musculus': 'GRCm39'}

fapi._set_session()


class TerraJobInputParams:
    def __init__(self, terra_data_table: pd.DataFrame, terra_uuids: TerraJobUUIDs, output_root_dir: str = '/igvf/data/'):
        self.terra_data_table = terra_data_table
        self.terra_uuids = terra_uuids
        self.terra_namespace = 'DACC_ANVIL'
        self.terra_workspace = 'IGVF Single-Cell Data Processing'
        self.output_root_dir = output_root_dir

    def _get_workflow_input_config(terra_namespace: str, terra_workspace: str, submission_id: str, workflow_id: str) -> dict:
        """Get the workflow input configuration for a given submission and workflow ID in a Terra workspace.

        Args:
            terra_namespace (str): Terra namespace (e.g., DACC_ANVIL)
            terra_workspace (str): Terra workspace name
            submission_id (str): submission ID for the full set of pipeline runs in the Terra workspace
            workflow_id (str): workflow ID (i.e., the ID per pipeline run)

        Raises:
            FireCloudServerError: If submission ID is not found or if workflow ID is not found in the submission.
            FireCloudServerError: Workflow input config is not found.

        Returns:
            dict: Workflow input configuration as a dictionary.
        """
        workflow_data_request = fapi.get_workflow_metadata(namespace=terra_namespace,
                                                           workspace=terra_workspace,
                                                           submission_id=submission_id,
                                                           workflow_id=workflow_id
                                                           )
    if workflow_data_request.status_code != 200:
        raise FireCloudServerError(code=workflow_data_request.status_code,
                                   message=f'Error fetching workflow data: {workflow_id} for submission: {submission_id}.')
    workflow_data_res = workflow_data_request.json()
    try:
        workflow_input_config = workflow_data_res.get('inputs')
        return workflow_input_config
    except KeyError:
        raise FireCloudServerError(code=workflow_data_request.status_code,
                                   message=f'Error parsing workflow input config for workflow: {workflow_id} in submission: {submission_id}. No inputs found.')

    def _download_workflow_config_json(self) -> str:
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
        terra_ids = self._parse_workflow_uuids_from_gs_path()
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


def mk_doc_aliases(curr_workflow_config: TerraJobUUIDs) -> list:
    """Create a list of document aliases for the workflow configuration.

    Args:
        curr_workflow_config (TerraJobUUIDs): The current workflow configuration

    Returns:
        list: A list of document aliases
    """
    return [f'igvf-dacc-processing-pipeline:{curr_workflow_config.config_aliases()}_pipeline_config']


def mk_qc_obj_aliases(curr_workflow_config: TerraJobUUIDs, analysis_set_acc: str, qc_prefix: str, lab: str) -> list:
    """Create a list of QC objects aliases for the workflow configuration.

    Args:
        curr_workflow_config (TerraJobUUIDs): The current workflow configuration
        analysis_set_acc (str): The analysis set accession
        qc_prefix (str): The prefix for the QC object (fragment, gene count, etc.)
        lab (str): The lab name

    Returns:
        list: A list of QC objects aliases
    """
    return [f'{lab.split("/")[-2]}:{analysis_set_acc}_{curr_workflow_config.aliases()}_{qc_prefix}_uniform-pipeline']

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
