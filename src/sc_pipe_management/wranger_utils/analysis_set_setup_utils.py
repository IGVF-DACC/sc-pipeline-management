
import logging

import accession.terra_to_portal_posting as t2portal
import wranger_utils.constant as const


# NOTE: The intent of this script is to be run in a Jupyter notebook so that the results can be easily inspected before posting new analysis sets to the portal.


def write_list_to_file(list_to_write: list, file_path: str) -> None:
    """Write a list of strings to a file, one per line.

    Args:
        list_to_write (list): A list of strings to write to the file.
        file_path (str): The path to the file where the list should be written.
    """
    with open(file_path, 'w') as f:
        for item in list_to_write:
            f.write(f"{item}\n")


def remove_duplicate_sublists(list_of_lists: list[list]) -> list[list]:
    """Remove duplicate sublists from a list of lists.

    Args:
        list_of_lists (list[list]): A list containing sublists.

    Returns:
        list[list]: A new list with duplicate sublists removed.
    """
    return [list(t) for t in set(tuple(sublist) for sublist in list_of_lists)]


def generate_filtered_fields(lab_id: str, award_id: str, preferred_assay_titles: list, statues: list, sample_names: list = []) -> dict:
    """Generate a dictionary of filtered fields for querying.

    Args:
        lab_id (str): e.g. /labs/j-michael-cherry/
        award_id (str): e.g. /awards/HG012012/
        preferred_assay_titles (list): e.g., ["10x multiome"]
        statues (list): a list of statuses to include, e.g., ['released']
        sample_names (list): a list of sample term names to include, e.g., ['induced pluripotent stem cell']. Defaults to [].

    Returns:
        dict: A dictionary of filtered fields for querying.
    """
    base_criteria = {'status': statues,
                     'lab.@id': lab_id,
                     'award.@id': award_id,
                     'preferred_assay_titles': preferred_assay_titles,
                     'audit.NOT_COMPLIANT.category!': const.EXCLUDED_NONCOMP_AUDITS,
                     'audit.ERROR.category!': const.EXCLUDED_ERROR_AUDITS
                     }
    if sample_names:
        base_criteria.update({'samples.sample_terms.term_name': sample_names})
    return base_criteria


def query_for_measets(filtered_fields: dict, igvf_client_api, limit: str | int = 'all') -> list | None:
    """Query for MeasurementSet objects based on filtered fields.

    Args:
        filtered_fields (dict): A dictionary of filtered fields for querying.
        igvf_client_api (_type_): IGVF client API instance.
        limit (str | int, optional): Query result limit. Defaults to 'all'.

    Returns:
        list | None: _description_
    """
    query_res = igvf_client_api.search(
        type=['MeasurementSet'], field_filters=filtered_fields, limit=limit)
    if query_res.total == 0:
        return None
    else:
        if limit == 'all':
            logging.info(f'>>>>>>>> {query_res.total} results found.')
        else:
            logging.info(f'>>>>>>>> {len(query_res.graph)} results found.')
        return query_res.graph


def check_is_uniform_workflow(workflow_id: str, igvf_client_api) -> bool:
    """Check if a workflow is a uniform workflow.

    Args:
        workflow_id (str): The ID of the workflow to check.
        igvf_client_api (_type_): The IGVF client API instance.

    Returns:
        bool: Whether the workflow is a uniform workflow.
    """
    workflow_obj = igvf_client_api.get_by_id(workflow_id).actual_instance
    if workflow_obj.uniform_pipeline:
        return workflow_obj.uniform_pipeline


def check_is_scpipeline(file_set_id: str, igvf_client_api) -> bool:
    """Check if a FileSet object is related to single-cell data. The assumption is that the Measurement Sets used to start all these checks are already single cell assays. But currently, there is no way to tell if an Analysis Set is created for the uniform pipeline. The uniform_pipeline property is calculated once files are linked.

    Args:
        file_set_id (str): e.g., /analysis-sets/IGVFAS0000xxxx/
        igvf_client_api (_type_): IGVF client API instance.
    Returns:
        bool: Whether the FileSet is related to single-cell data. True if it is, False otherwise.
    """
    # If not used as input for any analysis set, not a duplicated single cell uniform pipeline
    if not file_set_id.startswith('/analysis-sets/'):
        return False
    # Get the FileSet object by ID
    file_set_obj = igvf_client_api.get_by_id(file_set_id).actual_instance
    # If the FileSet has a uniform pipeline status, it is a single cell uniform pipeline analysis set
    if file_set_obj.uniform_pipeline_status is not None:
        return True
    # If the FileSet has uniform workflows, it is a single cell uniform pipeline analysis set
    if file_set_obj.workflows:
        return any([check_is_uniform_workflow(workflow_id=entry, igvf_client_api=igvf_client_api) for entry in file_set_obj.workflows])
    # If the FileSet has phrases indicating single cell in its aliases or description, it is a single cell uniform pipeline analysis set
    else:
        possible_descriptions = []
        # Check aliases
        if file_set_obj.aliases:
            possible_descriptions.extend(file_set_obj.aliases)
        if file_set_obj.description:
            possible_descriptions.append(file_set_obj.description)
        return any(kwd in entry for entry in possible_descriptions for kwd in const.POSSIBLE_SINGLE_CELL_PIPELINE_PHRASES)


def check_is_duped_for_all(measurement_set_id: str, igvf_client_api) -> bool:
    """Check if a MeasurementSet is already linked to a single cell uniform pipeline analysis set.

    Args:
        measurement_set_id (str): e.g., /measurement-sets/IGVFM0000xxxx/
        igvf_client_api (_type_): IGVF client API instance.

    Returns:
        bool: Whether the MeasurementSet is already linked to a single cell uniform pipeline analysis set.
    """
    measet_obj = igvf_client_api.measurement_sets(
        id=[measurement_set_id]).graph[0]
    if measet_obj.input_for:
        return any([check_is_scpipeline(file_set_id=entry, igvf_client_api=igvf_client_api) for entry in measet_obj.input_for])


def calc_input_file_sets_single(single_query_res, igvf_client_api) -> list | None:
    """ Calculate input file sets used to create a new single cell uniform pipeline analysis set from a measurement set ID.

    Args:
        single_query_res (_type_): a single IGVF client API MeasurementSet object query result
        igvf_client_api (_type_): IGVF client API instance.

    Returns:
        list | None: A sorted list of input file set IDs, or None if the analysis set is a duplicate.
    """
    # check if duplicating existing starter analysis sets
    if check_is_duped_for_all(measurement_set_id=single_query_res.id, igvf_client_api=igvf_client_api):
        return
    # collect a list of inputs if not duped
    curr_output = set()
    curr_output.add(single_query_res.id)
    if single_query_res.related_multiome_datasets:
        curr_output.update(single_query_res.related_multiome_datasets)
    if single_query_res.barcode_replacement_file:
        brf_obj = igvf_client_api.get_by_id(
            single_query_res.barcode_replacement_file).actual_instance
        curr_output.add(brf_obj.file_set)
    return sorted(curr_output)


def get_all_input_file_sets(query_res: list, igvf_client_api) -> list:
    """ Calculate input file sets used to create new single cell uniform pipeline analysis sets from a list of measurement set IDs.

    Args:
        query_res (list): a list of IGVF client API MeasurementSet object query results
        igvf_client_api (_type_): IGVF client API instance.

    Returns:
        list: A list of lists, each containing sorted input file set IDs for a new analysis set.
    """
    input_fileset_collection = []
    for res in query_res:
        input_fileset_collection.append(calc_input_file_sets_single(
            single_query_res=res, igvf_client_api=igvf_client_api))
    return remove_duplicate_sublists([entry for entry in input_fileset_collection if entry is not None])


def get_sample_accessions(input_file_sets: list, igvf_client_api) -> str:
    """ Get a concatenated string of unique sample accessions from a list of measurement set IDs.

    Args:
        input_file_sets (list): _description_
        igvf_client_api (_type_): _description_

    Returns:
        str: _description_
    """
    sample_accessions = set()
    for file_set_id in input_file_sets:
        # Parse will be curated sets
        if file_set_id.startswith('/measurement-sets/'):
            file_set_obj = igvf_client_api.get_by_id(file_set_id)
            curr_samples = [entry.split('/')[-2]
                            for entry in file_set_obj.samples]
            sample_accessions.update(curr_samples)
    return '_'.join(sorted(sample_accessions))


def create_analysis_set_payload(input_file_sets: list, lab: list, award: list, igvf_client_api) -> dict:
    """ Create a payload dictionary for creating a new analysis set.

    Args:
        input_file_sets (list): a list of input file set IDs for the new analysis set
        lab (list): e.g., /labs/j-michael-cherry/
        award (list): e.g., /awards/HG012012/
        igvf_client_api (_type_): IGVF client API instance.

    Returns:
        dict: A dictionary payload for creating a new analysis set.
    """
    alias_lab = lab.split('/')[-2]
    sample_accessions = get_sample_accessions(
        input_file_sets=input_file_sets, igvf_client_api=igvf_client_api)
    payload = {
        '_profile': 'analysis_set',
        'input_file_sets': input_file_sets,
        'lab': lab,
        'award': award,
        'aliases': [f'{alias_lab}:Sample-{sample_accessions}_single-cell-uniform-pipeline'],
        'uniform_pipeline_status': 'preprocessing',
        'file_set_type': 'intermediate analysis'
    }
    return payload


def create_all_analysis_set_payload(all_input_file_sets: list[list], lab: str, award: str, igvf_utils_api, igvf_client_api) -> list:
    """Post new analysis sets to the portal based on a list of input file sets.

    Args:
        all_input_file_sets (list[list]): a list of lists, each containing input file set IDs for a new analysis set
        lab (str): e.g., /labs/j-michael-cherry/
        award (str): e.g., /awards/HG012012/
        igvf_utils_api (_type_): IGVF utils API instance
        igvf_client_api (_type_): IGVF client API instance.

    Returns:
        list: A list of new analysis set accessions created.
    """
    post_res = []
    for input_file_sets in all_input_file_sets:
        curr_payload = create_analysis_set_payload(
            input_file_sets, lab, award, igvf_client_api)
        igvf_post_mthd = t2portal.IGVFPostService(igvf_utils_api=igvf_utils_api,
                                                  data_obj_payload=curr_payload,
                                                  upload_file=False,
                                                  resumed_posting=False
                                                  )
        post_res.append(igvf_post_mthd.single_post_to_portal())
    return post_res
