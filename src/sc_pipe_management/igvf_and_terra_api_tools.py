import firecloud.api as fapi
import subprocess
import pandas as pd
import io
import dumper
from igvf_client import Configuration
from igvf_client import ApiClient
from igvf_client import IgvfApi
from igvf_utils.connection import Connection
import os


from sc_pipe_management.portal_to_terra_input_from_anaset import mod_input_table_for_terra


fapi._set_session()


def set_up_api_keys(igvf_endpoint: str = 'sandbox') -> dict:
    """Set up API keys for IGVF portal using local env variables. Assumes that the environment variables named with each end site in mind. If not available, it will try to find a generic IGVF API key. If the generic IGVF API key is not set, it will raise an exception.

    Args:
        igvf_endpoint (str, optional): The IGVF endpoint to use. Defaults to 'sandbox'.

    Returns:
        dict: A dictionary containing the public and secret API keys.
    """
    api_keys_by_sites = {
        'production': {'public': 'IGVF_API_KEY_PROD', 'secret': 'IGVF_SECRET_KEY_PROD'},
        'sandbox': {'public': 'IGVF_API_KEY_SANDBOX', 'secret': 'IGVF_SECRET_KEY_SANDBOX'}
    }
    backup_api_keys = {'public': 'IGVF_API_KEY', 'secret': 'IGVF_SECRET_KEY'}

    api_keys = {}
    # If any of the 2 keys are missing from the pair
    if any(env_key not in os.environ.keys() for env_key in api_keys_by_sites[igvf_endpoint].values()):
        # If no backup generic IGVF API keys are set, raise an exception
        if any(env_key not in os.environ.keys() for env_key in backup_api_keys.values()):
            raise Exception(
                'Environment variables for IGVF public and secret keys are not set.')
        # If backup generic IGVF API keys are set, use them
        elif all(env_key in os.environ.keys() for env_key in backup_api_keys.values()):
            for purpose, env_key in backup_api_keys.items():
                api_keys[purpose] = os.environ[env_key]
    # If both keys are in the environment variables, use them
    elif all(env_key in os.environ.keys() for env_key in api_keys_by_sites[igvf_endpoint].values()):
        for purpose, env_key in api_keys_by_sites[igvf_endpoint].items():
            api_keys[purpose] = os.environ[env_key]
    return api_keys


def get_igvf_client_auth(igvf_api_keys: dict, igvf_site: str = 'sandbox'):
    """Set up IGVF data portal access and set up IGVF python client api.

    Args:
        igvf_api_keys (dict): A dictionary containing the public and secret API keys.
        igvf_site (str, optional): The IGVF site to use. Defaults to 'sandbox'.

    Returns:
        IgvfApi: An instance of the IgvfApi client.
    """
    site_urls_by_endpoints = {
        'production': 'https://api.data.igvf.org',
        'sandbox': 'https://api.sandbox.igvf.org'
    }
    config = Configuration(
        access_key=igvf_api_keys['public'],
        secret_access_key=igvf_api_keys['secret'],
        host=site_urls_by_endpoints[igvf_site],
    )
    client = ApiClient(config)
    return IgvfApi(client)


def get_igvf_utils_connection(igvf_api_keys: dict, igvf_utils_mode: str = 'sandbox', submission_mode: bool = False):
    """Set up IGVF utils connection.

    Args:
        igvf_api_keys (dict): A dictionary containing the public and secret API keys.
        igvf_utils_mode (str, optional): IGVF utils endpoint. Defaults to 'sandbox'.
        submission_mode (bool, optional): IGVF utils submission mode. Defaults to False, not submitting.

    Returns:
       igvf_utils.Connection: An instance of the IGVF utils connection.
    """
    iu_conn = Connection(igvf_mode=igvf_utils_mode, submission=submission_mode)
    iu_conn._auth = (igvf_api_keys['public'], igvf_api_keys['secret'])
    return iu_conn


def calculate_gsutil_hash(file_path: str):
    """Calculates the hash of a file using gsutil.

    Args:
        file_path (str): The path to the file for which to calculate the hash.

    Returns:
        str: The MD5 hash of the file.
    """
    result = subprocess.run(
        ['gsutil', 'hash', '-h', file_path], capture_output=True, text=True)
    if result.returncode == 0:
        # Parse the output to extract the hash values
        output_lines = result.stdout.splitlines()
        hash_values = {}
        for line in output_lines:
            if ':' in line:
                key, value = line.split(':')
                hash_values[key.strip()] = value.strip()
        return hash_values['Hash (md5)']
    else:
        raise Exception(f'gsutil hash command failed: {result.stderr}')


def get_terra_tsv_data_table(terra_namespace: str, terra_workspace: str, terra_etype: str) -> pd.DataFrame:
    """Get the data table from Terra workspace.

    Args:
        terra_namespace (str): DACC_ANVIL
        terra_workspace (str): workspace name
        terra_etype (str): Terra entity type name

    Returns:
        pd.DataFrame: The Terra data table
    """
    # TSV get just seems to fail whenever firecloud mode is used.
    query_response_for_tsv = fapi.get_entities_tsv(
        namespace=terra_namespace, workspace=terra_workspace, etype=terra_etype, model='flexible')
    return pd.read_csv(io.StringIO(query_response_for_tsv.content.decode('utf-8')), sep='\t')


def upload_portal_input_tsv_to_terra(terra_namespace: str, terra_workspace: str, terra_etype: str, portal_input_table: pd.DataFrame, local_barcode_file_dir: str, gs_barcode_list_bucket: str, verbose: bool = False):
    """Reformat input data table to Terra format. Upload this table to Terra workspace with user specified entity type name.

    Args:
        terra_namespace (str): DACC_ANVIL
        terra_workspace (str): workspace name
        terra_etype (str): Terra entity type name
        portal_input_table (pd.DataFrame): The input table generated from IGVF portal without any modification
        local_barcode_file_dir (str): Final barcode file directory
        gs_barcode_list_bucket (str): Where the barcode list files will be on gcloud
        verbose (bool, optional): Whether to show detailed log. Defaults to False.
    """
    portal_input_table_terra_format = mod_input_table_for_terra(pipeline_input_table=portal_input_table,
                                                                terra_etype=terra_etype,
                                                                local_barcode_file_dir=local_barcode_file_dir, gs_barcode_list_bucke=gs_barcode_list_bucket
                                                                )
    portal_input_table_as_string = portal_input_table_terra_format.to_csv(
        index=True, header=True, sep='\t')
    input_table_upload = fapi.upload_entities(namespace=terra_namespace,
                                              workspace=terra_workspace,
                                              entity_data=portal_input_table_as_string,
                                              model='flexible',
                                              delete_empty=False
                                              )
    if verbose:
        print(dumper.dump(input_table_upload))


def upload_output_post_res_to_terra(terra_namespace: str, terra_workspace: str, terra_etype: str, output_post_results: pd.DataFrame, verbose: bool = False):
    """Upload pipeline output POSTing to IGVF portal results and statues.

    Args:
        terra_namespace (str): DACC_ANVIL
        terra_workspace (str): workspace name
        terra_etype (str): Terra entity type name
        output_post_results (pd.DataFrame): The POSTing to IGVF portal report table
        verbose (bool, optional): Whether to show detailed log. Defaults to False.
    """
    output_post_results.index.name = f'{terra_etype}_id'
    output_results_table_as_string = output_post_results.to_csv(
        index=True, header=True, sep='\t')
    input_table_upload = fapi.upload_entities(namespace=terra_namespace,
                                              workspace=terra_workspace,
                                              entity_data=output_results_table_as_string,
                                              model='flexible'  # Firecloud mode just seems failing
                                              )
    if verbose:
        print(dumper.dump(input_table_upload))
