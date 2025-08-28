import pytest
import os
import sys
import pandas as pd

# Add the absolute path to the 'src' directory to sys.path
src_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.input_params_generation.constant as const
import sc_pipe_management.input_params_generation.generate_terra_input_table as terra_form
import sc_pipe_management.igvf_and_terra_api_tools as igvf_tools


# API tool
IGVF_ENDPOINT = 'prod'
IGVF_PROD_API_KEYS = igvf_tools.set_up_api_keys(igvf_endpoint=IGVF_ENDPOINT)
IGVF_PROD_CLIENT_API = igvf_tools.get_igvf_client_auth(igvf_api_keys=IGVF_PROD_API_KEYS,
                                                       igvf_endpoint=IGVF_ENDPOINT)


class TestCompleteTerraForming:

    def setup_method(self):
        # Expected Table
        self.expected_terra_table = pd.read_csv(
            'src/tests/test_files/test_input_datatable_generation_reference.tsv', sep='\t', index_col=0)
        # Analysis Sets
        self.test_analysis_set_accessions = [
            'IGVFDS9664YLAD',    # 10X multiome
            'IGVFDS0223KLTB',    # Parse SPLiT-seq
            'IGVFDS0657NHTA'     # ShareSeq
        ]
        # Generate Table
        self.test_terra_table = terra_form.CompleteTerraForming(
            analysis_set_accessions=self.test_analysis_set_accessions,
            igvf_api=IGVF_PROD_CLIENT_API,
            partial_root_dir='src/tests/test_files',
            terra_etype='unittest_pipeline_tester',
            local_barcode_file_dir=os.path.join(
                os.getcwd(), 'final_barcode_list/'),
            gs_barcode_list_bucket='gs://unittest_mock_bucket/submissions/final_barcode_onlist/'
        ).generate_complete_terra_input_table()

    def test_table_shape(self):
        assert self.expected_terra_table.shape == self.test_terra_table.shape

    def test_table_columns(self):
        assert sorted(list(self.expected_terra_table.columns)) == sorted(
            list(self.test_terra_table.columns))

    def test_table_rows(self):
        assert sorted(list(self.expected_terra_table.index)) == sorted(
            list(self.test_terra_table.index))

    def test_table_content(self):
        pd.testing.assert_frame_equal(
            self.expected_terra_table, self.test_terra_table)
