import firecloud.api as fapi
import subprocess
import pandas as pd
import io
import dumper


terra_output_table_columns = {
    'alignment_file': ['atac_bam'],
    'document': ['atac_bam_log',
                 'atac_chromap_barcode_metadata',
                 'atac_snapatac2_barcode_metadata',
                 'csv_summary',
                 'html_summary',
                 'joint_barcode_metadata',
                 'rna_barcode_metadata',
                 'rna_kb_output',
                 'rna_log'
                 ],
    'tabular_file': ['atac_filter_fragments',
                     'atac_filter_fragments_index'
                     ],
    'matrix_file': ['rna_aggregated_counts_h5ad',
                    'rna_mtx_tar',
                    'rna_mtxs_h5ad'
                    ]
}


fapi._set_session()


def calculate_gsutil_hash(file_path: str):
    """Calculates the hash of a file using gsutil."""
    result = subprocess.run(["gsutil", "hash", '-h', file_path], capture_output=True, text=True)
    if result.returncode == 0:
        # Parse the output to extract the hash values
        output_lines = result.stdout.splitlines()
        hash_values = {}
        for line in output_lines:
            if ":" in line:
                key, value = line.split(":")
                hash_values[key.strip()] = value.strip()
        return hash_values['Hash (md5)']
    else:
        raise Exception(f"gsutil hash command failed: {result.stderr}")


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
    # fapi._set_session()
    query_response_for_tsv = fapi.get_entities_tsv(namespace=terra_namespace, workspace=terra_workspace, etype=terra_etype, model='flexible')
    return pd.read_csv(io.StringIO(query_response_for_tsv.content.decode('utf-8')), sep='\t')


def upload_portal_input_tsv_to_terra(terra_namespace: str, terra_workspace: str, terra_etype: str, porta_input_table: pd.DataFrame, verbose: bool = False):
    """Upload an IGVF pipeline input table to Terra workspace with user specified entity type name.

    Args:
        terra_namespace (str): DACC_ANVIL
        terra_workspace (str): workspace name
        terra_etype (str): Terra entity type name
        porta_input_table (pd.DataFrame): The input table generated from IGVF portal
    """
    # fapi._set_session()
    porta_input_table.index.name = f'entity:{terra_etype}_id'
    porta_input_table_as_string = porta_input_table.to_csv(index=True, header=True, sep='\t')
    input_table_upload = fapi.upload_entities(namespace=terra_namespace,
                                              workspace=terra_workspace,
                                              entity_data=porta_input_table_as_string,
                                              model='flexible'  # Firecloud mode just seems failing
                                              )
    if verbose:
        print(dumper.dump(input_table_upload))
