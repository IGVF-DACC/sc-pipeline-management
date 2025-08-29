import sys
import os

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import argparse
from datetime import datetime
import os
import subprocess
import logging

import igvf_and_terra_api_tools as api_tools
import sc_pipe_management.input_params_generation.generate_terra_input_table as terra_forming


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
    # TODO: consider just defaulting to 'prod' and removing this argument
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
                                                     igvf_endpoint=args.igvf_endpoint)

    # Get a list of analysis set accessions from the arg or file, but not both
    if args.input_analysis_set_file:
        with open(args.input_analysis_set_file, 'r') as f:
            final_input_analysis_sets = ','.join(f.read().splitlines())
    elif args.input_analysis_set:
        final_input_analysis_sets = args.input_analysis_set
    logging.info(final_input_analysis_sets)

    # Where all the files generated will be stored locally
    partial_root_dir = os.path.join(
        os.getcwd(), "terra_datatables", "input", today, args.terra_etype)

    # Where barcode onlist files will be stored locally
    local_barcode_file_dir = os.path.join(
        partial_root_dir, "final_barcode_onlist")

    # Where barcode onlist files will be stored on Terra (make sure it doesn't end with a /)
    # Will be used to swap out the local dir in the final barcode file paths
    gs_barcode_list_bucket = f"gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/submissions/final_barcode_onlist/{today}/{args.terra_etype}"

    # Log all input args
    logging.info("Command-line arguments: %s", vars(args))

    # Generate the pipeline input table and save locally
    terra_data_table = terra_forming.CompleteTerraForming(
        analysis_set_accessions=final_input_analysis_sets.split(','),
        igvf_api=igvf_portal_api,
        partial_root_dir=partial_root_dir,
        terra_etype=args.terra_etype,
        local_barcode_file_dir=local_barcode_file_dir,
        gs_barcode_list_bucket=gs_barcode_list_bucket
    ).generate_complete_terra_input_table()

    logging.info("Pipeline input table generated and save locally.")

    # Save the pipeline input table locally
    terra_forming.save_pipeline_input_table(
        terra_data_table=terra_data_table,
        partial_root_dir=partial_root_dir)

    logging.info("Pipeline input table saved locally.")

    # If input generation has errors, log a warning and do not upload to Terra
    if any(x != "None" for x in terra_data_table['possible_errors'].values):
        logging.warning(
            'The pipeline input table has possible errors in parameter generation. Will not upload to Terra. Check the locally saved datable for details.')
    else:
        try:
            # Upload the pipeline input table to Terra
            api_tools.upload_portal_input_tsv_to_terra(
                terra_namespace=args.terra_namespace,
                terra_workspace=args.terra_workspace,
                portal_input_table=terra_data_table,
                verbose=True
            )
            logging.info("Pipeline input table uploaded to Terra.")

            # Copy barcode files to GCP bucket
            subprocess.run(
                ['gcloud', 'storage', 'rsync',
                    local_barcode_file_dir, gs_barcode_list_bucket],
                check=True
            )
            logging.info("Pipeline barcode files copied to GCP bucket.")
        except Exception as e:
            logging.error(str(e))
            raise Exception(
                "Error uploading to Terra or copying to GCP bucket: %s" % str(e))


if __name__ == '__main__':
    main()
