import sys
import os
import hashlib

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

    def setup_method(self):
        """Setup method to initialize any required variables or state before each test.
        """
        test_seqspec_metadata_mthds = {}
        # preferred assay titles and seqspecs by assay type
        for test_assay_title, test_seqspec_info in test_const.SEQSPEC_PARSE_TEST_RESULTS.items():
            test_seqspec_metadata_mthds[test_assay_title] = {}
            # assay types and their corresponding test and reference seqspec info
            for assay_type, per_assay_seqspec_info in test_seqspec_info.items():

                # Get the reference seqspec info
                test_ref_seqspec_info = test_const.SEQSPEC_PARSE_TEST_RESULTS[
                    test_assay_title][assay_type]

                # Generate the seqspec metadata for testing
                curr_seqspec_metadata_mthd = seqspec_parsing.GetSeqSpecMetadata(
                    seqspec_file_path=per_assay_seqspec_info.seqspec_file_path,
                    igvf_api=IGVF_PROD_CLIENT_API
                )
                test_seqspec_metadata_mthds[test_assay_title][assay_type] = {
                    "test": curr_seqspec_metadata_mthd,
                    "reference": test_ref_seqspec_info
                }
        return test_seqspec_metadata_mthds

    def test_seqspec_metadata_generation_idv_steps(self):
        """Test if the modality is correctly parsed from seqspec file.
        """
        test_seqspec_metadata_mthds = self.setup_method()
        # preferred assay titles and seqspecs by assay type
        for test_assay_title, test_seqspec_info in test_seqspec_metadata_mthds.items():
            # assay types and their corresponding test and reference seqspec info
            for assay_type, test_pairs in test_seqspec_info.items():
                # Check if modality match
                curr_seqspec_modality = test_pairs['test']._get_seqspec_modality(
                )
                assert curr_seqspec_modality == test_pairs['reference'].modality, \
                    f"Modality mismatch for {test_assay_title} {assay_type}: expected {test_pairs['test'].modality}, got {curr_seqspec_modality}"

                # Check if read IDs match
                curr_ordered_read_ids = test_pairs['test']._generate_ordered_read_ids(
                    seqfiles_metadata=test_const.TEST_SEQFILES_METADATA, assay_type=assay_type)
                assert curr_ordered_read_ids == test_pairs['reference'].ordered_read_ids, \
                    f"Ordered read IDs mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].ordered_read_ids}, got {curr_ordered_read_ids}"

                # Check if onlist files match
                curr_onlist_files = test_pairs['test']._get_onlist_files()
                assert curr_onlist_files == sorted(test_pairs['reference'].onlist_files), \
                    f"Onlist files mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].onlist_files}, got {curr_onlist_files}"

    def test_generate_seqspec_metadata_integration(self):
        """Integration test for generate_seqspec_metadata."""
        test_seqspec_metadata_mthds = self.setup_method()
        # preferred assay titles and seqspecs by assay type
        for test_assay_title, test_seqspec_info in test_seqspec_metadata_mthds.items():
            # assay types and their corresponding test and reference seqspec info
            for assay_type, test_pairs in test_seqspec_info.items():
                # Generate the seqspec metadata
                curr_seqspec_metadata = test_pairs['test'].generate_seqspec_metadata(
                    seqfiles_metadata=test_const.TEST_SEQFILES_METADATA, assay_type=assay_type
                )

                # Check if modality match
                assert curr_seqspec_metadata.modality == test_pairs['reference'].modality, \
                    f"Modality mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].modality}, got {curr_seqspec_metadata.modality}"

                # Check if read IDs match
                assert curr_seqspec_metadata.ordered_read_ids == test_pairs['reference'].ordered_read_ids, \
                    f"Ordered read IDs mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].ordered_read_ids}, got {curr_seqspec_metadata.ordered_read_ids}"

                # Check if onlist files match
                assert sorted(curr_seqspec_metadata.onlist_files) == sorted(test_pairs['reference'].onlist_files), \
                    f"Onlist files mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].onlist_files}, got {curr_seqspec_metadata.onlist_files}"


class TestGetSeqSpecToolOutput:
    """Test class for GetSeqSpecToolOutput.
    """

    def setup_method(self):
        """Setup method to initialize any required variables or state before each test.
        """
        test_seqspec_tool_mthds = {}
        for test_assay_title, test_seqspec_info in test_const.SEQSPEC_PARSE_TEST_RESULTS.items():
            test_seqspec_tool_mthds[test_assay_title] = {}
            for assay_type, per_assay_seqspec_info in test_seqspec_info.items():
                # Get the reference seqspec info
                test_ref_seqspec_info = test_const.SEQSPEC_PARSE_TEST_RESULTS[
                    test_assay_title][assay_type]

                converted_test_seqspec_metadata = test_ref_seqspec_info.convert_to_seqspec_metadata()

                # Generate the seqspec metadata for testing
                curr_seqspec_tool_mthd = seqspec_parsing.GetSeqSpecToolOutput(
                    seqspec_metadata=converted_test_seqspec_metadata,
                    onlist_method=test_ref_seqspec_info.onlist_method,
                    output_barcode_list_file=test_ref_seqspec_info.final_barcode_file
                )
                test_seqspec_tool_mthds[test_assay_title][assay_type] = {
                    "test": curr_seqspec_tool_mthd,
                    "reference": test_ref_seqspec_info
                }
        return test_seqspec_tool_mthds

    def test_get_seqspec_tool_output(self):
        """Test if the seqspec tool output is correctly generated.
        """
        test_seqspec_tool_mthds = self.setup_method()
        for test_assay_title, test_seqspec_info in test_seqspec_tool_mthds.items():
            for assay_type, test_pairs in test_seqspec_info.items():
                # Generate test seqspec tool method
                curr_seqspec_tool_mthd = test_pairs['test']

                # Check if read index string match
                curr_seqspec_index = curr_seqspec_tool_mthd._get_seqspec_index()
                assert curr_seqspec_index == test_pairs['reference'].read_index, \
                    f"Read index string mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].read_index}, got {curr_seqspec_index}"

                # Check if onlist file MD5 match
                curr_final_barcode_file = curr_seqspec_tool_mthd._get_seqspec_onlist()
                curr_final_barcode_file_md5 = compute_md5sum(
                    curr_final_barcode_file)
                assert curr_final_barcode_file_md5 == test_pairs['reference'].final_barcode_file_md5, \
                    f"Onlist file path mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].final_barcode_file_md5}, got {curr_final_barcode_file_md5}"

    def test_get_seqspec_tool_output_integration(self):
        """Integration test for GetSeqSpecToolOutput."""
        test_seqspec_tool_mthds = self.setup_method()
        for test_assay_title, test_seqspec_info in test_seqspec_tool_mthds.items():
            for assay_type, test_pairs in test_seqspec_info.items():
                # Generate test seqspec tool output
                curr_seqspec_tool_output = test_pairs['test'].generate_seqspec_tool_output(
                )

                # Check if read index string match
                assert curr_seqspec_tool_output.read_index_string == test_pairs['reference'].read_index, \
                    f"Read index string mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].read_index}, got {curr_seqspec_tool_output.read_index_string}"

                # Check if onlist file MD5 match
                curr_final_barcode_file_md5 = compute_md5sum(
                    curr_seqspec_tool_output.final_barcode_file)
                assert curr_final_barcode_file_md5 == test_pairs['reference'].final_barcode_file_md5, \
                    f"Onlist file path mismatch for {test_assay_title} {assay_type}: expected {test_pairs['reference'].final_barcode_file_md5}, got {curr_final_barcode_file_md5}"

                # Check if errors is None
                assert curr_seqspec_tool_output.errors is None, \
                    f"Errors is not None for {test_assay_title} {assay_type}: got {curr_seqspec_tool_output.errors}"
