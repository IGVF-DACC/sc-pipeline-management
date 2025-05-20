import igvf_and_terra_api_tools as api_tools
import portal_to_terra_input_from_anaset as portal2terra_transfer
import argparse
import firecloud.api as fapi
from datetime import datetime
import os
import subprocess
import pandas
import logging


def setup_logging(log_file: str):
    """Generate a log file for the script.

    Args:
        log_file (str): Log file path.
    """
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler(log_file),
                            logging.StreamHandler(),
                        ])


def get_parser():
    """
    Get command line arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    # IGVF query endoint
    parser.add_argument('--igvf_endpoint', type=str, choices=['sandbox', 'prod', 'staging'],
                        help="""The IGVF endpoint, sandbox, prod, or staging.""")
    # Either input an analysis set or a file containing a list of analysis sets
    group.add_argument('--input_analysis_set', type=str,
                       help="""A list of analysis set accessions to be processed.""")
    group.add_argument('--input_analysis_set_file', type=str,
                       help="""A file containing a list of analysis set accessions to be processed.""")
    # Project name and workspace name for the pipeline platform
    parser.add_argument('--terra_namespace', type=str, default='DACC_ANVIL',
                        help="""Terra namespace name for the pipeline platform.""")
    parser.add_argument('--terra_workspace', type=str, default='IGVF Single-Cell Data Processing',
                        help="""Terra workspace name for the pipeline platform.""")
    # Terra eType name to modify the index column
    parser.add_argument('--terra_etype', type=str,
                        help="""Terra eType name for the pipeline data table.""")
    # Where the local final barcodes are stored
    parser.add_argument('--local_barcode_file_dir', type=str, default=None,
                        help=""""Local directory for the pipeline final barcode files. Defaults to $(pwd)/final_barcode_list/$(date +%m%d%Y)")""")
    # Where on GCP the barcode files are stored
    parser.add_argument('--gs_barcode_list_bucket', type=str, default=None,
                        help="""Google Storage bucket for the pipeline final barcode files. Defaults to gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/submissions/final_barcode_onlist/$(date +%m%d%Y)""")
    # Output directory for the pipeline input table
    parser.add_argument('--output_dir', type=str, default=None,
                        help="""Path to the output directory.""")
    return parser


def main():
    # Today's date and logging setup
    today = datetime.now().strftime("%m%d%Y")
    log_dir = './Run_Logs'
    os.makedirs(log_dir, exist_ok=True)
    setup_logging(log_file=os.path.join(
        os.getcwd(), log_dir, f'portal_to_terra_input_generation_{today}.log'))
    logging.info('Run Date: %s', today)

    # Set up arg parser
    parser = get_parser()
    args = parser.parse_args()

    # Get the IGVF client API, production or sandbox
    igvf_api_keys = api_tools.set_up_api_keys(
        igvf_endpoint=args.igvf_endpoint)

    # Set up IGVF client API
    igvf_portal_api = api_tools.get_igvf_client_auth(igvf_api_keys=igvf_api_keys,
                                                     igvf_site=args.igvf_endpoint)
    # Refresh Firecloud API
    fapi._set_session()

    # Get a list of analysis set accessions from the arg or file, but not both
    if args.input_analysis_set_file:
        with open(args.input_analysis_set_file, 'r') as f:
            final_input_analysis_sets = ','.join(f.read().splitlines())
    elif args.input_analysis_set:
        final_input_analysis_sets = args.input_analysis_set
    else:
        analysis_input_error_msg = 'Please provide either --input_analysis_set or --input_analysis_set_file.'
        logging.error(analysis_input_error_msg)
        raise ValueError(analysis_input_error_msg)
    logging.info(final_input_analysis_sets)

    # Get default barcode file directory if not provided
    if args.local_barcode_file_dir is None:
        args.local_barcode_file_dir = os.path.join(
            os.getcwd(), "final_barcode_list", today)

    # Get default GCP bucket if not provided
    if args.gs_barcode_list_bucket is None:
        args.gs_barcode_list_bucket = f"gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/submissions/final_barcode_onlist/{today}"

    # Get default input parameter table dir
    if args.output_dir is None:
        args.output_dir = os.path.join(
            os.getcwd(), "terra_datatables", "input", today)

    # Log all input args
    logging.info("Command-line arguments: %s", vars(args))

    # Generate the pipeline input table
    table = portal2terra_transfer.generate_pipeline_input_table(
        query_analysis_set_accs=final_input_analysis_sets.split(','),
        terra_etype=args.terra_etype,
        local_barcode_file_dir=args.local_barcode_file_dir,
        gs_barcode_list_bucket=args.gs_barcode_list_bucket,
        igvf_api=igvf_portal_api)
    logging.info("Pipeline input table generated.")

    # Output the table to local folder as a Copy
    portal2terra_transfer.save_pipeline_input_table(
        pipeline_input_table=table, output_dir=os.path.join(args.output_dir, args.terra_etype))
    logging.info("Pipeline input table saved to local folder.")

    try:
        # Upload the pipeline input table to Terra
        api_tools.upload_portal_input_tsv_to_terra(
            terra_namespace=args.terra_namespace,
            terra_workspace=args.terra_workspace,
            portal_input_table=table,
            verbose=True
        )
        logging.info("Pipeline input table uploaded to Terra.")

        # Copy barcode files to GCP bucket
        subprocess.run(
            ['gcloud', 'storage', 'cp', '-r',
                args.local_barcode_file_dir, args.gs_barcode_list_bucket],
            check=True
        )
        logging.info("Pipeline barcode files copied to GCP bucket.")
    except Exception as e:
        logging.error(str(e))
        raise Exception(
            "Error uploading to Terra or copying to GCP bucket: %s" % str(e))


if __name__ == '__main__':
    main()
