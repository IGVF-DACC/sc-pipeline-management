import igvf_and_terra_api_tools as api_tools
import terra_to_portal_posting as terra2portal_transfer
import argparse
import firecloud.api as fapi
import os
from datetime import datetime


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
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    terra_table = api_tools.get_terra_tsv_data_table(terra_namespace=args.terra_namespace,
                                                     terra_workspace=args.terra_workspace,
                                                     terra_etype=args.terra_etype)

    # Get the API keys
    igvf_api_keys = api_tools.set_up_api_keys(
        igvf_endpoint=args.post_endpoint)

    # Set up igvf utils connection (staging requires using the site URL)
    if args.post_endpoint != 'staging':
        igvf_utils_api = api_tools.get_igvf_utils_connection(igvf_api_keys=igvf_api_keys,
                                                             igvf_utils_mode=args.post_endpoint,
                                                             submission_mode=True)
    else:
        igvf_utils_api = api_tools.get_igvf_utils_connection(igvf_api_keys=igvf_api_keys,
                                                             igvf_utils_mode=api_tools.SITE_URLS_BY_ENDPOINTS[
                                                                 args.post_endpoint],
                                                             submission_mode=False)

    # Set up igvf client connection
    igvf_client_api = api_tools.get_igvf_client_auth(igvf_api_keys=igvf_api_keys,
                                                     igvf_endpoint=args.post_endpoint)

    # Refresh firecloud API
    fapi._set_session()

    # Set up local output directory
    today = datetime.now().strftime("%m%d%Y")
    if args.output_dir is None:
        args.output_dir = os.path.join(
            os.getcwd(), "terra_datatables", "output", today)

    # Will return a list of Postres
    portal_post_results = terra2portal_transfer.post_all_successful_runs(full_terra_data_table=terra_table,
                                                                         terra_namespace=args.terra_namespace,
                                                                         terra_workspace=args.terra_workspace,
                                                                         igvf_api=igvf_client_api,
                                                                         igvf_utils_api=igvf_utils_api,
                                                                         upload_file=args.upload_file,
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
