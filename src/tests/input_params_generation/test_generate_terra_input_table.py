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

    @pytest.fixture
    def expected_terra_table(self):
        return pd.read_csv(
            'src/tests/test_files/test_input_datatable_generation_reference.tsv', sep='\t', index_col=0)

    @pytest.fixture
    def test_terra_table_imported(self):
        test_analysis_set_accessions = [
            'IGVFDS9664YLAD',
            'IGVFDS0223KLTB',
            'IGVFDS0657NHTA'
        ]
        test_terra_table = terra_form.CompleteTerraForming(
            analysis_set_accessions=test_analysis_set_accessions,
            igvf_api=IGVF_PROD_CLIENT_API,
            partial_root_dir='src/tests/test_files',
            terra_etype='unittest_pipeline_tester',
            local_barcode_file_dir=os.path.join(
                os.getcwd(), 'final_barcode_list/'),
            gs_barcode_list_bucket='gs://unittest_mock_bucket/submissions/final_barcode_onlist/'
        ).generate_complete_terra_input_table()
        test_terra_table.to_csv(
            'src/tests/test_files/test_input_datatable_generation_output.tsv', sep='\t')
        return pd.read_csv(
            'src/tests/test_files/test_input_datatable_generation_output.tsv', sep='\t', index_col=0)

    def test_table_shape(self, expected_terra_table, test_terra_table_imported):
        assert expected_terra_table.shape == test_terra_table_imported.shape

    def test_table_columns(self, expected_terra_table, test_terra_table_imported):
        assert sorted(list(expected_terra_table.columns)) == sorted(
            list(test_terra_table_imported.columns))

    def test_table_rows(self, expected_terra_table, test_terra_table_imported):
        assert sorted(list(expected_terra_table.index)) == sorted(
            list(test_terra_table_imported.index))

    def test_table_content(self, expected_terra_table, test_terra_table_imported):
        pd.testing.assert_frame_equal(
            expected_terra_table, test_terra_table_imported, check_like=True)
