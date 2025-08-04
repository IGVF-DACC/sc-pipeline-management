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
    return [f'{lab.split("/")[-2]}:{analysis_set_acc}_{curr_workflow_config.aliases()}_{qc_prefix}_uniform-pipeline']
