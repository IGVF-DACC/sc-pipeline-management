import igvf_and_terra_api_tools as api_tools
import terra_to_portal_posting as terra2portal_transfer
import argparse
import logging


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--post_endpoint', type=str,
                        help="""The POST endpoint, sandbox or production.""")
    parser.add_argument('--terra_namespace', type=str,
                        help="""Terra namespace name for the pipeline platform.""")
    parser.add_argument('--terra_workspace', type=str,
                        help="""Terra workspace name for the pipeline platform.""")
    parser.add_argument('--terra_etype', type=str,
                        help="""Terra eType name for the pipeline data table.""")
    parser.add_argument('--upload_file', action='store_true',
                        help="""If True, upload the file to the portal.""")
    parser.add_argument('--output_dir', type=str,
                        help="""Path to the output directory.""")
    return parser


def setup_logging(log_file_path):
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler(log_file_path),
                            logging.StreamHandler(),
                        ])


def main():
    parser = get_parser()
    args = parser.parse_args()
    setup_logging(log_file_path='./Run_Logs/terra_to_portal_data_post.log')
    terra_table = api_tools.get_terra_tsv_data_table(terra_namespace=args.terra_namespace,
                                                     terra_workspace=args.terra_workspace,
                                                     terra_etype=args.terra_etype)
    logging.info(f'Got the Terra data table with shape: {terra_table.shape}')
    igvf_utils_api = api_tools.get_igvf_utils_connection(
        igvf_utils_mode=args.post_endpoint)
    igvf_portal_api = api_tools.get_igvf_auth_and_api(
        igvf_site=args.post_endpoint)
    # Will return a list of Postres
    portal_post_results = terra2portal_transfer.post_all_successful_runs(full_terra_data_table=terra_table,
                                                                         igvf_api=igvf_portal_api,
                                                                         igvf_utils_api=igvf_utils_api,
                                                                         upload_file=args.upload_file)
    # Summarize into a table
    portal_post_summary = terra2portal_transfer.summarize_post_status(
        post_results=portal_post_results)
    updated_terra_table = terra2portal_transfer.add_post_status_summary_to_output_data_table(
        full_terra_data_table=terra_table, post_status_df=portal_post_summary)
    terra2portal_transfer.save_pipeline_postres_tables(
        pipeline_postres_table=portal_post_summary, updated_full_data_table=updated_terra_table, output_dir=args.output_dir)


if __name__ == '__main__':
    main()
