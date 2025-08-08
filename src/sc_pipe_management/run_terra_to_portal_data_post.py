import sys
import os

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)


import igvf_and_terra_api_tools as api_tools
import sc_pipe_management.accession.igvf_payloads as igvf_payloads
import sc_pipe_management.accession.parse_terra_metadata as terra_parse
import sc_pipe_management.accession.constants as const
import sc_pipe_management.accession.terra_to_portal_posting as tr2igvf
import argparse
import firecloud.api as fapi
import os
from datetime import datetime
from functools import wraps
import logging


def read_exclusion_file(file_path: str) -> list:
    """Read exclusion file and return list of accessions.

    Args:
        file_path (str): Path to the exclusion file containing accessions to exclude.

    Raises:
        argparse.ArgumentTypeError: If the file cannot be found or read, an error is raised with a message indicating the issue.
        argparse.ArgumentTypeError: If the file is empty or contains only whitespace, an error is raised.

    Returns:
        list: A list of accessions to exclude from the data table or an empty list if no file is provided.
    """
    if file_path is None:
        return []
    try:
        with open(file_path, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        raise argparse.ArgumentTypeError(
            f"Exclusion file not found: {file_path}")
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Error reading exclusion file: {e}")


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--post_endpoint', type=str, choices=['sandbox', 'prod', 'staging'],
                        help="""The POST endpoint, sandbox, prod, or staging.""")
    parser.add_argument('--terra_namespace', type=str, default='DACC_ANVIL',
                        help="""Terra namespace name for the pipeline platform. Defaults to DACC_ANVIL""")
    parser.add_argument('--terra_workspace', type=str, default='IGVF Single-Cell Data Processing',
                        help="""Terra workspace name for the pipeline platform. Defaults to IGVF Single-Cell Data Processing""")
    parser.add_argument('--terra_etype', type=str,
                        help="""Terra eType name for the pipeline data table.""")
    parser.add_argument('--excluded_accs', type=read_exclusion_file, default=None,
                        help="""A file with a list of accessions to exclude from the data table. If not provided, no accessions will be excluded.""")
    parser.add_argument('--upload_file', action='store_true',
                        help="""If True, upload the file to the portal.""")
    parser.add_argument('--output_dir', type=str, default=None,
                        help="""Path to the output directory. Defaults to $(pwd)/terra_datatables/output/$(date +%m%d%Y)""")
    parser.add_argument('--resumed_posting', action='store_true',
                        help="""Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.""")
    parser.add_argument('--tries', type=int, default=3,
                        help="""Number of tries to attempt the API calls. Default is 3.""")
    parser.add_argument('--delay', type=int, default=5,
                        help="""Initial delay between retries in seconds. Default is 5.""")
    parser.add_argument('--backoff', type=int, default=2,
                        help="""Backoff multiplier e.g. value of 2 will double the delay each retry. Default is 2.""")
    return parser


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
                    logging.debug(str(e))
                    msg = "Retrying in %d seconds..." % (mdelay)
                    logging.info(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry


def main():
    parser = get_parser()
    args = parser.parse_args()

    tries = args.tries
    delay = args.delay
    backoff = args.backoff

    @retry(tries, delay, backoff)
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

    # Set up igvf client connection
    @retry(tries, delay, backoff)
    def do_igvf_client_api(igvf_api_keys, post_endpoint):
        """Set up IGVF client API connection
        """
        return api_tools.get_igvf_client_auth(igvf_api_keys=igvf_api_keys,
                                              igvf_endpoint=post_endpoint)

    # Refresh firecloud API
    @retry(tries, delay, backoff)
    def do_firecloud_api():
        """Set up Fire Cloud API connection
        """
        fapi._set_session()

    # Set up local output directory
    today = datetime.now().strftime("%m%d%Y")
    if args.output_dir is None:
        args.output_dir = os.path.join(
            os.getcwd(), "terra_datatables", "output", today, args.terra_etype)

    # Call FireCloud API
    do_firecloud_api()

    # Get the IGVF API keys
    igvf_api_keys = api_tools.set_up_api_keys(
        igvf_endpoint=args.post_endpoint)

    # Call IGVF APIs
    igvf_client_api = do_igvf_client_api(
        igvf_api_keys=igvf_api_keys, post_endpoint=args.post_endpoint)
    igvf_utils_api = do_igvf_utils_api(
        igvf_api_keys=igvf_api_keys, post_endpoint=args.post_endpoint)

    # Download data table from Terra
    terra_table = api_tools.get_terra_tsv_data_table(terra_namespace=args.terra_namespace,
                                                     terra_workspace=args.terra_workspace,
                                                     terra_etype=args.terra_etype,
                                                     excluded_accessions=args.excluded_accs)
    if args.excluded_accs is None:
        num_excluded = 0
    else:
        num_excluded = len(args.excluded_accs)
    logging.info(
        f'>>>>>>>>>>>>>> A total of {terra_table.shape[0]} pipeline runs found with {num_excluded} of which excluded.')

    # Download all workflow configs first (firecloud times out)
    pipeline_params_info = igvf_payloads.PipelineParamsInfo(
        terra_namespace=args.terra_namespace,
        terra_workspace=args.terra_workspace,
        terra_datatable=terra_table,
        igvf_client_api=igvf_client_api,
        output_root_dir=os.path.join(args.output_dir, 'workflow_configs'))
    anaset_input_params_file_paths = pipeline_params_info.get_all_input_params()
    logging.info(
        f'>>>>>>>>>>>>>> {len(anaset_input_params_file_paths)} configs downloaded')

    # Will return a list of accessioning results
    portal_post_results = tr2igvf.post_all_pipeline_runs_from_one_submission(terra_data_table=terra_table,
                                                                             igvf_client_api=igvf_client_api,
                                                                             anaset_input_params_file_paths=anaset_input_params_file_paths,
                                                                             igvf_utils_api=igvf_utils_api,
                                                                             output_root_dir=args.output_dir,
                                                                             upload_file=args.upload_file,
                                                                             resumed_posting=args.resumed_posting)

    logging.info(
        f'>>>>>>>>>>>>>> A total of {len(portal_post_results)} pipeline runs posted to the portal.')

    # Summarize into a table
    portal_post_summary = tr2igvf.summarize_post_status(
        run_results=portal_post_results)
    updated_terra_table = tr2igvf.add_post_status_summary_to_output_data_table(
        full_terra_data_table=terra_table,
        post_status_df=portal_post_summary)

    # Save the updated terra table
    tr2igvf.save_pipeline_postres_tables(
        pipeline_postres_table=portal_post_summary,
        updated_full_data_table=updated_terra_table,
        output_root_dir=args.output_dir)

    logging.info(
        f'>>>>>>>>>>>>>> Pipeline post results saved to {args.output_dir}')


if __name__ == '__main__':
    main()
