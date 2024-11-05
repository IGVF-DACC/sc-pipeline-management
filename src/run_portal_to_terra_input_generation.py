import igvf_and_terra_api_tools as api_tools
import portal_to_terra_input_from_anaset as portal2terra_transfer
import argparse
import logging


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--igvf_endpoint', type=str,
                        help="""The IGVF endpoint, sandbox or production.""")
    parser.add_argument('--input_analysis_set', type=str,
                        help="""A list of analysis set accessions to be processed.""")
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
    igvf_portal_api = api_tools.get_igvf_auth_and_api(
        igvf_site=args.igvf_endpoint)
    print(args.input_analysis_set)
    table = portal2terra_transfer.generate_pipeline_input_table(query_analysis_set_accs=args.input_analysis_set.split(','),
                                                                igvf_api=igvf_portal_api)
    portal2terra_transfer.save_pipeline_input_table(
        pipeline_input_table=table, output_dir=args.output_dir)


if __name__ == '__main__':
    main()
