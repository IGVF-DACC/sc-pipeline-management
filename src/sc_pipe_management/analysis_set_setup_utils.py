import igvf_utils.connection as iu_conn
import re


ANALYSIS_SET_REGEX = re.compile(r'analysis_sets\/(.*)\/')


def get_analysis_set_obj(analysis_set_id: str, igvf_util_conn) -> dict:
    """Get analysis set object from IGVF

    Args:
        analysis_set_id (str): IGVF analysis set accession
        igvf_util_conn (_type_): IGVF utils connection object

    Returns:
        dict: analysis set object
    """
    return igvf_util_conn.get(analysis_set_id)


def check_if_has_files(analysis_set_obj: dict) -> bool:
    """Check if the analysis set has result files

    Args:
        analysis_set_obj (dict): analysis set object

    Returns:
        bool: _description_
    """
    if analysis_set_obj['files']:
        return True
    else:
        return False


def draft_analysis_set_payload(analysis_set_obj: dict) -> dict:
    """Draft analysis set payload for new analysis set POSTing

    Args:
        analysis_set_obj (dict): analysis set object

    Returns:
        dict: new analysis set payload for POSTing
    """
    analysis_set_payload = dict(
        _profile='analysis_set',
        input_file_sets=[file_set['@id']
                         for file_set in analysis_set_obj['input_file_sets']],
        lab=analysis_set_obj['lab']['@id'],
        award=analysis_set_obj['award']['@id'],
        file_set_type='intermediate analysis'
    )
    return analysis_set_payload


# TODO: Probably should have an exception handling for the case where the analysis set has no input file sets
def singlecell_analysis_set_setup(analysis_set_id: str, igvf_util_conn) -> str:
    """Clone or keep starter analysis set for single cell uniform pipeline setup based on files status.

    Args:
        analysis_set_id (str): analysis set id
        igvf_util_conn (_type_): IGVF utils connection object

    Returns:
        str: Either a new analysis set accession or the original analysis set accession (if no files)
    """
    analysis_set_obj = get_analysis_set_obj(analysis_set_id, igvf_util_conn)
    # If has files, clone it
    if check_if_has_files(analysis_set_obj):
        new_payload = draft_analysis_set_payload(analysis_set_obj)
        _schema_property = igvf_util_conn.get_profile_from_payload(
            new_payload).properties
        stdout = igvf_util_conn.post(
            new_payload, upload_file=False, return_original_status_code=True, require_aliases=False)
        return stdout[0]['accession']
    else:
        # If no files, keep it as is
        return analysis_set_id


def batch_singlecell_analysis_set_setup(analysis_set_ids: list, igvf_util_conn) -> list:
    """Batch run clone or keep starter analysis sets for single cell uniform pipeline setup.

    Args:
        analysis_set_id (list): A list of analysis set accessions from single cell progress tracker
        igvf_util_conn (_type_): IGVF utils connection object

    Returns:
        list: A list of new or original analysis set accessions for starters
    """
    single_cell_starter_analysis_sets = []
    for analysis_set_id in analysis_set_ids:
        single_cell_starter_analysis_sets.append(
            singlecell_analysis_set_setup(analysis_set_id, igvf_util_conn))
    return single_cell_starter_analysis_sets
