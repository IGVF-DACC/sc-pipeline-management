import sys
import os

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)


import hashlib
import unittest
import logging


def compute_md5sum(file_path: str) -> str:
    """Compute the MD5 checksum of a file.

    Args:
        file_path (str): File of interest to compute the md5sum for.

    Returns:
        str: MD5 checksum of the file in hexadecimal format.
    """
    md5_hash = hashlib.md5()
    with open(file_path, 'rb') as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b''):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()


# API tool
IGVF_ENDPOINT = 'prod'
IGVF_PROD_API_KEYS = igvf_tools.set_up_api_keys(igvf_endpoint=IGVF_ENDPOINT)
IGVF_PROD_CLIENT_API = igvf_tools.get_igvf_client_auth(igvf_api_keys=IGVF_PROD_API_KEYS,
                                                       igvf_endpoint=IGVF_ENDPOINT)


class TestSeqspecInfoGeneration:
    def test_generate_read_id_input(self):
        """Test if read id input -i arg is generated correctly.
        """
        for assay_title, seqspec_paths in SEQSPEC_FILES_BY_ASSAY_TITLES.items():
            for assay_type, seqspec_path in seqspec_paths.items():
                if seqspec_path is None:
                    # Skip if the seqspec path is None (e.g., SPLiT-seq for atac)
                    logging.info(
                        f'Skipping {assay_title} for {assay_type} as it has no seqspec.')
                    continue
                for read_id_methd in ['onlist', 'index']:
                    logging.info(
                        f'Testing {assay_title} {assay_type} with read_id_methd: {read_id_methd}.')
                    logging.info(seqspec_path)
                    curr_read_input = generate_ordered_read_ids(seqspec_file_path=seqspec_path,
                                                                assay_type=assay_type,
                                                                usage_purpose=read_id_methd,
                                                                igvf_api=IGVF_PROD_CLIENT_API
                                                                )
                    curr_error_detail = (f'Generated index read ids do not match expected for {assay_title} {assay_type} with method {read_id_methd}. '
                                         f'Expected: {SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type][f"{read_id_methd}_input"]}, '
                                         f' Got: {curr_read_input}.'
                                         )
                    self.assertEqual(
                        first=curr_read_input, second=SEQSPEC_OUTPUT_BY_ASSAY_TITLES[
                            assay_title][assay_type][f'{read_id_methd}_input'],
                        msg=curr_error_detail)
        logging.info('All read id input tests done.\n')

    def test_read_index_generation(self):
        """Check if rna and atac read index are correct.
        """
        for assay_title, seqspec_paths in SEQSPEC_FILES_BY_ASSAY_TITLES.items():
            for assay_type, seqspec_path in seqspec_paths.items():
                if seqspec_path is None:
                    # Skip if the seqspec path is None (e.g., SPLiT-seq for atac)
                    logging.info(
                        f'Skipping {assay_title} for {assay_type} as it has no seqspec.')
                    continue
                logging.info(
                    f'Testing {assay_title} {assay_type} for read index generation.')
                # Get the current read index from seqspec
                logging.info(seqspec_path)
                curr_index_input = seqspec_index_get(
                    seqspec_file_path=seqspec_path, assay_type=assay_type, igvf_api=IGVF_PROD_CLIENT_API)
                curr_error_detail = (f'Generated {assay_type} index read index does not match expected for {assay_title} {assay_type}. '
                                     f'Expected: {SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type]["read_index"]}, '
                                     f' Got: {curr_index_input}.'
                                     )
                self.assertEqual(
                    first=curr_index_input,
                    second=SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type]['read_index'],
                    msg=curr_error_detail
                )
        logging.info('All read index tests done.\n')

    def test_finalinclusion_list_generation(self):
        """Check if the final inclusion list file has correct md5sum for each assay type.
        """
        for assay_title, seqspec_paths in SEQSPEC_FILES_BY_ASSAY_TITLES.items():
            for assay_type, seqspec_path in seqspec_paths.items():
                if seqspec_path is None:
                    # Skip if the seqspec path is None (e.g., SPLiT-seq for atac)
                    logging.info(
                        f'Skipping {assay_title} for {assay_type} as it has no seqspec.')
                    continue
                logging.info(
                    f'Testing final inclusion list generation for {assay_title} {assay_type}.')
                logging.info(seqspec_path)
                output_file_path = generate_finalinclusion_list(
                    seqspec_file_path=seqspec_path,
                    assay_type=assay_type,
                    onlist_method=SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type]['onlist_method'],
                    final_inclusion_list_path=os.path.join(
                        'src/tests/test_files/', f'{assay_title}_{assay_type}_final_barcode_list.txt'),
                    igvf_api=IGVF_PROD_CLIENT_API
                )
                curr_md5sum = compute_md5sum(
                    file_path=output_file_path
                )
                # Compare the md5sum of the generated file with the expected md5sum
                expected_md5sum = SEQSPEC_OUTPUT_BY_ASSAY_TITLES[
                    assay_title][assay_type]['onlist_final_list']
                curr_error_detail = (f'Generated final inclusion list md5sum does not match expected for {assay_title} {assay_type}. '
                                     f'Expected: {expected_md5sum}, '
                                     f' Got: {curr_md5sum}.'
                                     )
                self.assertEqual(
                    first=curr_md5sum,
                    second=expected_md5sum,
                    msg=curr_error_detail
                )
        logging.info('All final inclusion list tests done.\n')
