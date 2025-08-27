import sys
import os
import hashlib
import logging

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
project_root = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', '..', '..'))
src_path = os.path.join(project_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)


import sc_pipe_management.igvf_and_terra_api_tools as igvf_tools
import sc_pipe_management.input_params_generation.seqspec_parsing as seqspec_parsing
import constant as test_const


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
    """Test class for seqspec parsing functions.
    """

    def test_seqspec_metadata_generation_idv_steps(self):
        """Test if the modality is correctly parsed from seqspec file.
        """
        for test_file_type, test_seqspec_info in test_const.SEQSPEC_PARSE_TEST_RESULTS.items():
            for assay_type, per_assay_seqspec_info in test_seqspec_info.items():

                # Get the reference seqspec info
                test_ref_seqspec_info = test_const.SEQSPEC_PARSE_TEST_RESULTS[
                    test_file_type][assay_type]

                # Generate the seqspec metadata for testing
                curr_seqspec_metadata_mthd = seqspec_parsing.GetSeqSpecMetadata(
                    seqspec_file_path=per_assay_seqspec_info.seqspec_file_path,
                    igvf_api=IGVF_PROD_CLIENT_API
                )

                # Check if modality match
                curr_seqspec_modality = curr_seqspec_metadata_mthd._get_seqspec_modality()
                assert curr_seqspec_modality == test_ref_seqspec_info.modality, \
                    f"Modality mismatch for {test_file_type} {assay_type}: expected {test_ref_seqspec_info.modality}, got {curr_seqspec_modality}"

                # Check if read IDs match
                curr_ordered_read_ids = curr_seqspec_metadata_mthd._generate_ordered_read_ids(
                    seqfiles_metadata=test_const.TEST_SEQFILES_METADATA)
                assert curr_ordered_read_ids == test_ref_seqspec_info.ordered_read_ids, \
                    f"Ordered read IDs mismatch for {test_file_type} {assay_type}: expected {test_ref_seqspec_info.ordered_read_ids}, got {curr_ordered_read_ids}"

                # Check if onlist files match
                curr_onlist_files = curr_seqspec_metadata_mthd._get_onlist_files()
                assert curr_onlist_files == sorted(test_ref_seqspec_info.onlist_files), \
                    f"Onlist files mismatch for {test_file_type} {assay_type}: expected {test_ref_seqspec_info.onlist_files}, got {curr_onlist_files}"

    def test_generate_seqspec_metadata_integration(self):
        """Integration test for generate_seqspec_metadata."""
        for test_file_type, test_seqspec_info in test_const.SEQSPEC_PARSE_TEST_RESULTS.items():
            for assay_type, per_assay_seqspec_info in test_seqspec_info.items():
                test_ref_seqspec_info = test_const.SEQSPEC_PARSE_TEST_RESULTS[
                    test_file_type][assay_type]

                curr_seqspec_metadata_mthd = seqspec_parsing.GetSeqSpecMetadata(
                    seqspec_file_path=per_assay_seqspec_info.seqspec_file_path,
                    igvf_api=IGVF_PROD_CLIENT_API
                )
                curr_seqspec_metadata = curr_seqspec_metadata_mthd.generate_seqspec_metadata(
                    seqfiles_metadata=test_const.TEST_SEQFILES_METADATA
                )

                assert curr_seqspec_metadata.modality == test_ref_seqspec_info.modality, \
                    f"Modality mismatch for {test_file_type} {assay_type}: expected {test_ref_seqspec_info.modality}, got {curr_seqspec_metadata.modality}"

                assert curr_seqspec_metadata.ordered_read_ids == test_ref_seqspec_info.ordered_read_ids, \
                    f"Ordered read IDs mismatch for {test_file_type} {assay_type}: expected {test_ref_seqspec_info.ordered_read_ids}, got {curr_seqspec_metadata.ordered_read_ids}"

                # If onlist_files may be unsorted, sort both before comparing
                assert sorted(curr_seqspec_metadata.onlist_files) == sorted(test_ref_seqspec_info.onlist_files), \
                    f"Onlist files mismatch for {test_file_type} {assay_type}: expected {test_ref_seqspec_info.onlist_files}, got {curr_seqspec_metadata.onlist_files}"
