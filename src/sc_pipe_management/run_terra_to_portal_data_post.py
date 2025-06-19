import igvf_and_terra_api_tools as api_tools
import terra_to_portal_posting as terra2portal_transfer
import argparse
import firecloud.api as fapi
import os
from datetime import datetime
from functools import wraps


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
    parser.add_argument('--upload_file', action='store_true',
                        help="""If True, upload the file to the portal.""")
    parser.add_argument('--output_dir', type=str, default=None,
                        help="""Path to the output directory. Defaults to $(pwd)/terra_datatables/output/$(date +%m%d%Y)""")
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
            os.getcwd(), "terra_datatables", "output", today)

    # Call FireCloud API
    do_firecloud_api()

    # Download data table from Terra
    terra_table = api_tools.get_terra_tsv_data_table(terra_namespace=args.terra_namespace,
                                                     terra_workspace=args.terra_workspace,
                                                     terra_etype=args.terra_etype)

    # Download all workflow configs first (firecloud times out)
    config_file_collection = terra2portal_transfer.download_all_workflow_config_jsons(
        terra_namespace=args.terra_namespace,
        terra_workspace=args.terra_workspace,
        terra_data_table=terra_table,
        output_root_dir=os.path.join(args.output_dir, 'workflow_configs'))
    print(f'>>>>>>>>>>>>>> {len(config_file_collection)} configs downloaded')

    # Get the IGVF API keys
    igvf_api_keys = api_tools.set_up_api_keys(
        igvf_endpoint=args.post_endpoint)

    # Call IGVF APIs
    igvf_client_api = do_igvf_client_api(
        igvf_api_keys=igvf_api_keys, post_endpoint=args.post_endpoint)
    igvf_utils_api = do_igvf_utils_api(
        igvf_api_keys=igvf_api_keys, post_endpoint=args.post_endpoint)

    # Will return a list of Postres
    portal_post_results = terra2portal_transfer.post_all_successful_runs(full_terra_data_table=terra_table,
                                                                         igvf_api=igvf_client_api,
                                                                         igvf_utils_api=igvf_utils_api,
                                                                         upload_file=args.upload_file,
                                                                         config_file_collection=config_file_collection,
                                                                         output_root_dir=args.output_dir)

    # Summarize into a table
    portal_post_summary = terra2portal_transfer.summarize_post_status(
        post_results=portal_post_results)
    updated_terra_table = terra2portal_transfer.add_post_status_summary_to_output_data_table(
        full_terra_data_table=terra_table, post_status_df=portal_post_summary)

    # Save the updated terra table
    terra2portal_transfer.save_pipeline_postres_tables(
        pipeline_postres_table=portal_post_summary, updated_full_data_table=updated_terra_table, output_root_dir=args.output_dir)


if __name__ == '__main__':
    main()
