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


fapi._set_session()


def get_igvf_auth_and_api(igvf_site: str = 'sandbox'):
    """Set up IGVF data portal access and set up IGVF python client api.
    """
    if igvf_site == 'production':
        host = 'https://api.data.igvf.org'
    elif igvf_site == 'sandbox':
        host = 'https://api.sandbox.igvf.org'
    config = Configuration(
        access_key=os.environ['IGVF_API_KEY'],
        secret_access_key=os.environ['IGVF_SECRET_KEY'],
        host=host
    )
    client = ApiClient(config)
    return IgvfApi(client)


def get_igvf_utils_connection(igvf_utils_mode: str = 'sandbox'):
    return Connection(igvf_mode=igvf_utils_mode, submission=True)


def calculate_gsutil_hash(file_path: str):
    """Calculates the hash of a file using gsutil."""
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


def upload_portal_input_tsv_to_terra(terra_namespace: str, terra_workspace: str, terra_etype: str, portal_input_table: pd.DataFrame, verbose: bool = False):
    """Upload an IGVF pipeline input table to Terra workspace with user specified entity type name.

    Args:
        terra_namespace (str): DACC_ANVIL
        terra_workspace (str): workspace name
        terra_etype (str): Terra entity type name
        portal_input_table (pd.DataFrame): The input table generated from IGVF portal
        verbose (bool, optional): Whether to show detailed log. Defaults to False.
    """
    portal_input_table.index = portal_input_table['analysis_set_acc']
    portal_input_table.index.name = f'entity:{terra_etype}_id'
    portal_input_table_as_string = portal_input_table.to_csv(
        index=True, header=True, sep='\t')
    input_table_upload = fapi.upload_entities(namespace=terra_namespace,
                                              workspace=terra_workspace,
                                              entity_data=portal_input_table_as_string,
                                              model='flexible'  # Firecloud mode just seems failing
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
