import sc_pipe_management.igvf_and_terra_api_tools as api_tools
import sc_pipe_management.portal_to_terra_input_from_anaset as portal2terra_transfer
import argparse
import logging
import firecloud.api as fapi


def get_parser():
    """
    Get command line arguments for the script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    # IGVF query endoint
    parser.add_argument('--igvf_endpoint', type=str,
                        help="""The IGVF endpoint, sandbox or production.""")
    # Either input an analysis set or a file containing a list of analysis sets
    group.add_argument('--input_analysis_set', type=str,
                       help="""A list of analysis set accessions to be processed.""")
    group.add_argument('--input_analysis_set_file', type=str,
                       help="""A file containing a list of analysis set accessions to be processed.""")
    # Terra eType name to modify the index column
    parser.add_argument('--terra_etype', type=str,
                        help="""Terra eType name for the pipeline data table.""")
    # Where the local final barcodes are stored
    parser.add_argument('--local_barcode_file_dir', type=str,
                        help="""Local directory for the pipeline final barcode files.""")
    # Where on GCP the barcode files are stored
    parser.add_argument('--gs_barcode_list_bucket', type=str,
                        help="""Google Storage bucket for the pipeline final barcode files.""")
    # Output directory for the pipeline input table
    parser.add_argument('--output_dir', type=str,
                        help="""Path to the output directory.""")
    return parser


def setup_logging(log_file):
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler(log_file),
                            logging.StreamHandler(),
                        ],
                        filename='./Run_Logs/porta_to_terra_input_generation.log')


def main():
    parser = get_parser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

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
        raise ValueError(
            'Please provide either --input_analysis_set or --input_analysis_set_file.')
    print(final_input_analysis_sets)

    # Generate the pipeline input table
    table = portal2terra_transfer.generate_pipeline_input_table(
        query_analysis_set_accs=final_input_analysis_sets.split(','),
        terra_etype=args.terra_etype,
        local_barcode_file_dir=args.local_barcode_file_dir,
        gs_barcode_list_bucket=args.gs_barcode_list_bucket,
        igvf_api=igvf_portal_api)

    # Output the table to local
    portal2terra_transfer.save_pipeline_input_table(
        pipeline_input_table=table, output_dir=args.output_dir)


if __name__ == '__main__':
    main()
