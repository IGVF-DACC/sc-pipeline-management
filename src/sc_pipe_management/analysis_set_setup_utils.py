import igvf_and_terra_api_tools as api_tools
import terra_to_portal_posting as t2portal
import argparse
from functools import wraps
from datetime import datetime
import requests


def limit_type(value):
    """Custom type function that accepts 'all' or converts to int."""
    if value.lower() == 'all':
        return 'all'
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"limit must be 'all' or an integer, got '{value}'")


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--post_endpoint', type=str, choices=['sandbox', 'prod', 'staging'],
                        help="""The POST endpoint, sandbox, prod, or staging.""")
    parser.add_argument('--lab', type=str,
                        help="""Lab ID for the analysis set.""")
    parser.add_argument('--award', type=str,
                        help="""Award ID for the analysis set.""")
    parser.add_argument('--preferred_assay_title', type=str,
                        help="""Preferred assay title for the measurement set query.""")
    parser.add_argument('--limit', type=limit_type, default='all',
                        help="""Query result limit. Defaults to all.""")
    parser.add_argument('--output_dir', type=str, default='./',
                        help="""Path to the output directory to save the new analysis set accessions.""")
    return parser


# TODO: this will be replaced by AnalysisSet pipeline status once the ticket is in
FILE_TRACKING_PIPELINE_RUNS = 'src/sc_pipe_management/files/running_list_of_finished_analysis_sets.txt'


POSSIBLE_SINGLE_CELL_PIPELINE_PHRASES = ['scpipe',
                                         'single cell',
                                         'single cell uniform pipeline',
                                         'uniform pipeline',
                                         'sc-pipe',
                                         'scpipeline',
                                         'sc-uniform-pipeline',
                                         'sc-pipeline',
                                         'uniform-pipeline']


FINISHED_ANALYSIS_SETS = open(
    FILE_TRACKING_PIPELINE_RUNS, 'r').read().splitlines()

EXCLUDED_NONCOMP_AUDITS = ['missing sequence specification',
                           'missing barcode replacement file',
                           'missing read names',
                           'missing barcode onlist']

EXCLUDED_ERROR_AUDITS = ['upload status not validated',
                         'unexpected barcode onlist',
                         'inconsistent sequence specifications',
                         'inconsistent preferred assay title']

MEASET_STATUSES = ['released', 'preview']


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


def generate_filtered_fields(lab_id: str, award_id: str, preferred_assay_titles: list, excluded_nc_audits: list, excluded_error_audits: list, statues: list, sample_names: list = []) -> dict:
    """Generate a dictionary of filtered fields for querying.

    Args:
        lab_id (str): e.g. /labs/j-michael-cherry/
        award_id (str): e.g. /awards/HG012012/
        preferred_assay_titles (list): e.g., ["10x multiome"]
        excluded_nc_audits (list): a list of NOT_COMPLIANT audit categories to exclude
        excluded_error_audits (list): a list of ERROR audit categories to exclude
        statues (list): a list of statuses to include, e.g., ['released']
        sample_names (list): a list of sample term names to include, e.g., ['induced pluripotent stem cell']. Defaults to [].

    Returns:
        dict: A dictionary of filtered fields for querying.
    """
    base_criteria = {'status': statues,
                     'lab.@id': lab_id,
                     'award.@id': award_id,
                     'preferred_assay_titles': preferred_assay_titles,
                     'audit.NOT_COMPLIANT.category!': excluded_nc_audits,
                     'audit.ERROR.category!': excluded_error_audits
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
            print(f'>>>>>>>> {query_res.total} results found.')
        else:
            print(f'>>>>>>>> {len(query_res.graph)} results found.')
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


def check_is_scpipeline(file_set_obj, igvf_client_api) -> bool:
    """Check if a FileSet object is related to single-cell data. The assumption is that the Measurement Sets used to start all these checks are already single cell assays. But currently, there is no way to tell if an Analysis Set is created for the uniform pipeline. The uniform_pipeline property is calculated once files are linked.

    Args:
        file_set_obj (_type_): a single IGVF client API FileSet object query result

    Returns:
        bool: Whether the FileSet is related to single-cell data.
    """
    if file_set_obj.workflows:
        return any([check_is_uniform_workflow(workflow_id=entry, igvf_client_api=igvf_client_api) for entry in file_set_obj.workflows])
    else:
        possible_descriptions = []
        # Check aliases
        if file_set_obj.aliases:
            possible_descriptions.extend(file_set_obj.aliases)
        if file_set_obj.description:
            possible_descriptions.append(file_set_obj.description)
        return any(kwd in entry for entry in possible_descriptions for kwd in POSSIBLE_SINGLE_CELL_PIPELINE_PHRASES)


def check_is_duped_scpipeline(file_set_id: str, igvf_client_api, finished_anasets: list = []) -> bool:
    """Check if an analysis set is a single cell uniform pipeline related and hasn't been processed.

    Args:
        file_set_id (str): e.g., /analysis-sets/IGVFAS0000xxxx/
        igvf_client_api (_type_): IGVF client API instance.
        finished_anasets (list, optional): a list of analysis sets already processed. Defaults to [].

    Returns:
        bool: Whether the analysis set is a single cell uniform pipeline, have no duplicates, and haven't been processed.
    """
    if not file_set_id.startswith('/analysis-sets/'):
        return False
    if file_set_id.split('/')[-2] in finished_anasets:
        return True
    sub_file_set_obj = igvf_client_api.get_by_id(file_set_id).actual_instance
    if check_is_scpipeline(file_set_obj=sub_file_set_obj, igvf_client_api=igvf_client_api):
        return True
    else:
        return False


def check_is_duped_for_all(measurement_set_id: str, igvf_client_api, finished_anasets: str = []) -> bool:
    """Check if a MeasurementSet is already linked to a single cell uniform pipeline analysis set.

    Args:
        measurement_set_id (str): e.g., /measurement-sets/IGVFM0000xxxx/
        igvf_client_api (_type_): IGVF client API instance.
        finished_anasets (str, optional): a list of analysis sets already processed. Defaults to [].

    Returns:
        bool: Whether the MeasurementSet is already linked to a single cell uniform pipeline analysis set.
    """
    measet_obj = igvf_client_api.measurement_sets(
        id=[measurement_set_id]).graph[0]
    if measet_obj.input_for:
        return any([check_is_duped_scpipeline(file_set_id=entry, igvf_client_api=igvf_client_api, finished_anasets=finished_anasets) for entry in measet_obj.input_for])


def calc_input_file_sets_single(single_query_res, igvf_client_api, finished_anasets: str = []) -> list | None:
    """ Calculate input file sets used to create a new single cell uniform pipeline analysis set from a measurement set ID.

    Args:
        single_query_res (_type_): a single IGVF client API MeasurementSet object query result
        igvf_client_api (_type_): IGVF client API instance.
        finished_anasets (str, optional): a list of analysis sets already processed. Defaults to [].

    Returns:
        list | None: A sorted list of input file set IDs, or None if the analysis set is a duplicate.
    """
    # check if duplicating existing starter analysis sets
    if check_is_duped_for_all(measurement_set_id=single_query_res.id, igvf_client_api=igvf_client_api, finished_anasets=finished_anasets):
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


def get_all_input_file_sets(query_res: list, igvf_client_api, finished_anasets: str = []) -> list:
    """ Calculate input file sets used to create new single cell uniform pipeline analysis sets from a list of measurement set IDs.

    Args:
        query_res (list): a list of IGVF client API MeasurementSet object query results
        igvf_client_api (_type_): IGVF client API instance.
        finished_anasets (str, optional): a list of analysis sets already processed. Defaults to [].

    Returns:
        list: A list of lists, each containing sorted input file set IDs for a new analysis set.
    """
    input_fileset_collection = []
    for res in query_res:
        input_fileset_collection.append(calc_input_file_sets_single(
            single_query_res=res, igvf_client_api=igvf_client_api, finished_anasets=finished_anasets))
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
        post_res.append(t2portal.single_post_to_portal(
            igvf_data_payload=curr_payload, igvf_utils_api=igvf_utils_api, upload_file=False))
    return post_res


## decorator for preventing time out ##
def retry(tries=1, delay=5, backoff=2):
    import time
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    args:
        tries (int): number of times to try (not retry) before giving up
        delay (int): initial delay between retries in seconds
        backoff (int): backoff multiplier e.g. value of 2 will double the delay each retry
    """
    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except Exception as e:
                    print(str(e))
                    msg = "Retrying in %d seconds..." % (mdelay)
                    print(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry


def main():
    parser = get_parser()
    args = parser.parse_args()

    # Set up igvf client connection
    @retry(tries=1, delay=3, backoff=2)
    def do_igvf_client_api(igvf_api_keys, post_endpoint):
        """Set up IGVF client API connection
        """
        return api_tools.get_igvf_client_auth(igvf_api_keys=igvf_api_keys,
                                              igvf_endpoint=post_endpoint)

    @retry(tries=1, delay=3, backoff=2)
    def do_igvf_utils_api(igvf_api_keys, post_endpoint):
        """Set up IGVF utils API connection
        """
        if post_endpoint != 'staging':
            return api_tools.get_igvf_utils_connection(igvf_api_keys=igvf_api_keys,
                                                       igvf_utils_mode=post_endpoint,
                                                       submission_mode=True)
        else:
            return api_tools.get_igvf_utils_connection(igvf_api_keys=igvf_api_keys,
                                                       igvf_utils_mode=api_tools.SITE_URLS_BY_ENDPOINTS[
                                                           post_endpoint],
                                                       submission_mode=True)

    # Get the IGVF API keys
    igvf_api_keys = api_tools.set_up_api_keys(
        igvf_endpoint=args.post_endpoint)

    # Call IGVF APIs
    igvf_client_api = do_igvf_client_api(
        igvf_api_keys=igvf_api_keys, post_endpoint=args.post_endpoint)
    igvf_utils_api = do_igvf_utils_api(
        igvf_api_keys=igvf_api_keys, post_endpoint=args.post_endpoint)

    # Query for measurement sets
    filtered_fields = generate_filtered_fields(
        lab_id=args.lab,
        award_id=args.award,
        # query one assay title at a time
        preferred_assay_titles=[args.preferred_assay_title],
        excluded_nc_audits=EXCLUDED_NONCOMP_AUDITS,
        excluded_error_audits=EXCLUDED_ERROR_AUDITS,
        statues=MEASET_STATUSES
    )

    # Do the query
    query_res = query_for_measets(
        filtered_fields=filtered_fields, igvf_client_api=igvf_client_api, limit=args.limit)
    if query_res is None:
        print('>>>>>>>>>> No measurement sets found. Exiting.')
        return

    # Get all input file sets for new analysis sets
    all_input_file_sets = get_all_input_file_sets(
        query_res=query_res, igvf_client_api=igvf_client_api, finished_anasets=FINISHED_ANALYSIS_SETS)
    if not all_input_file_sets:
        print('>>>>>>>>>> No new analysis sets to create. Exiting.')
        return

    # Post new analysis sets to the portal
    try:
        post_res = create_all_analysis_set_payload(
            all_input_file_sets=all_input_file_sets,
            lab=args.lab,
            award=args.award,
            igvf_utils_api=igvf_utils_api,
            igvf_client_api=igvf_client_api
        )
    except requests.exceptions.HTTPError as e:
        print(f'>>>>>>>>>>> Error posting new analysis sets: {e}')
        return

    # Save the new analysis set accessions to a file
    unique_name = f'{args.lab.split("/")[-2]}_{args.preferred_assay_title.replace(" ", "-")}_{datetime.now().strftime("%m%d%Y")}'
    output_file_path = f"{args.output_dir.rstrip('/')}/new_anasets_to_run_{unique_name}.txt"
    write_list_to_file(list_to_write=post_res, file_path=output_file_path)
    print(
        f'>>>>>>>>>> New analysis set accessions saved to {output_file_path}')


if __name__ == "__main__":
    main()
