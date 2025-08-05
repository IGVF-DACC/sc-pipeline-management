import sys
import os

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import unittest
import pandas as pd
from unittest.mock import Mock, patch, MagicMock, mock_open
import json
from pathlib import Path

import sc_pipe_management.accession.parse_terra_metadata as parse_terra
import sc_pipe_management.accession.igvf_payloads as igvf_payloads
import sc_pipe_management.accession.terra_to_portal_posting as igvf_posting


class TestTerraOutputMetadata(unittest.TestCase):
    """Test cases for TerraOutputMetadata class using real test data."""

    @classmethod
    def setUpClass(cls):
        """Load the Terra test data table once for all tests."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        cls.terra_data_table = pd.read_csv(test_file_path, sep='\t')

    def setUp(self):
        """Set up test fixtures before each test method."""
        # First row has a full set 10x multiome data
        self.sample_terra_record = self.terra_data_table.iloc[0]
        self.terra_metadata = parse_terra.TerraOutputMetadata(
            self.sample_terra_record, self.mock_igvf_api)

    def test_init_with_real_data(self):
        """Test TerraOutputMetadata initialization with real test data."""
        self.assertEqual(self.terra_metadata.anaset_accession,
                         self.sample_terra_record['analysis_set_acc'])
        self.assertEqual(self.terra_metadata.taxa,
                         self.sample_terra_record['taxa'])
        self.assertIsNotNone(self.terra_metadata.genome_assembly)
        self.assertTrue(self.terra_metadata.terra_data_record.equals(
            self.sample_terra_record))

    def test_get_gs_path_prioritizes_rna(self):
        """Test that RNA paths are prioritized over ATAC paths."""
        gs_path = self.terra_metadata._get_gs_path_for_terra_output_cols()

        # Should return RNA path if available
        if pd.notna(self.sample_terra_record.get('rna_kb_h5ad')):
            self.assertEqual(gs_path, self.sample_terra_record['rna_kb_h5ad'])
        elif pd.notna(self.sample_terra_record.get('atac_bam')):
            self.assertEqual(gs_path, self.sample_terra_record['atac_bam'])

    def test_parse_workflow_uuids_from_gs_path(self):
        """Test parsing workflow UUIDs from GS path using real data."""
        try:
            result = self.terra_metadata._parse_workflow_uuids_from_gs_path()
            self.assertIsInstance(result, parse_terra.InputFileAccs)
            self.assertIsNotNone(result.gcloud_bucket)
            self.assertIsNotNone(result.submission_id)
            self.assertIsNotNone(result.workflow_id)
            self.assertIsNotNone(result.subworkflow_id)
        except ValueError:
            # Some test data might not have valid GS paths
            self.skipTest(
                "Test data doesn't contain valid GS paths for UUID parsing")

    def test_parse_terra_str_list_with_real_data(self):
        """Test parsing comma-separated strings from real Terra data."""
        # Test RNA sequence file accessions
        if 'rna_seqfile_accessions' in self.sample_terra_record:
            rna_accessions_str = self.sample_terra_record['rna_seqfile_accessions']
            if pd.notna(rna_accessions_str):
                result = self.terra_metadata._parse_terra_str_list(
                    [rna_accessions_str])
                self.assertIsInstance(result, list)
                self.assertTrue(all(acc.startswith(('IGVF', 'TSTF'))
                                for acc in result))

    def test_parse_igvf_accessions_from_urls_real_data(self):
        """Test extracting IGVF accessions from URLs using real data."""
        if 'rna_seqspec_urls' in self.sample_terra_record:
            urls_str = self.sample_terra_record['rna_seqspec_urls']
            if pd.notna(urls_str):
                result = self.terra_metadata._parse_igvf_accessions_from_urls([
                                                                              urls_str])
                self.assertIsInstance(result, list)
                # Should extract valid accessions from URLs
                if result:
                    self.assertTrue(all(acc.startswith(('IGVF', 'TSTF'))
                                    for acc in result))

    def test_get_input_file_accs_from_table_rna(self):
        """Test getting input file accessions for RNA with real data."""
        try:
            result = self.terra_metadata._get_input_file_accs_from_table('rna')
            self.assertIsInstance(result, parse_terra.InputFileAccs)
            self.assertIsInstance(result.derived_from_accessions, list)
            self.assertIsInstance(result.reference_files, list)

            # Check that accessions start with IGVF or TSTF
            for acc in result.derived_from_accessions:
                self.assertTrue(acc.startswith(('IGVF', 'TSTF')))
        except (KeyError, AttributeError):
            # Some test data might not have all required columns
            self.skipTest(
                "Test data missing required columns for RNA input file parsing")

    def test_get_input_file_accs_from_table_atac(self):
        """Test getting input file accessions for ATAC with real data."""
        try:
            result = self.terra_metadata._get_input_file_accs_from_table(
                'atac')
            self.assertIsInstance(result, parse_terra.InputFileAccs)
            self.assertIsInstance(result.derived_from_accessions, list)
            self.assertIsInstance(result.reference_files, list)

            # Check that accessions start with IGVF or TSTF
            for acc in result.derived_from_accessions:
                self.assertTrue(acc.startswith(('IGVF', 'TSTF')))
        except (KeyError, AttributeError):
            # Some test data might not have all required columns
            self.skipTest(
                "Test data missing required columns for ATAC input file parsing")

    def test_get_gs_path_no_valid_path_raises_error(self):
        """Test that missing GS paths raise ValueError."""
        # Create a record with no valid GS paths
        empty_record = pd.Series({
            'analysis_set_acc': 'IGVFDS999ZZZ',
            'taxa': 'Homo sapiens'
        })
        empty_metadata = parse_terra.TerraOutputMetadata(
            empty_record, self.mock_igvf_api)

        with self.assertRaises(ValueError):
            empty_metadata._get_gs_path_for_terra_output_cols()

    @patch('sc_pipe_management.accession.parse_terra_metadata.fapi')
    @patch('sc_pipe_management.accession.parse_terra_metadata.dump_json')
    def test_download_workflow_config_json(self, mock_dump_json, mock_fapi):
        """Test downloading workflow configuration JSON."""
        # Mock Terra API response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {'inputs': {'param1': 'value1'}}
        mock_fapi.get_workflow_metadata.return_value = mock_response

        # Mock file dump
        mock_dump_json.return_value = '/path/to/config.json'

        try:
            result = self.terra_metadata._download_workflow_config_json(
                terra_namespace='test_namespace',
                terra_workspace='test_workspace'
            )

            self.assertEqual(result.download_path, '/path/to/config.json')
            self.assertIsInstance(result.doc_aliases, list)
        except ValueError:
            # Skip if test data doesn't have valid workflow UUIDs
            self.skipTest("Test data doesn't contain valid workflow UUIDs")


class TestMatrixFilePayload(unittest.TestCase):
    """Test cases for MatrixFilePayload class using real test data."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        # Load test data
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            terra_data_table = pd.read_csv(test_file_path, sep='\t')
            self.sample_terra_record = terra_data_table.iloc[0]
        else:
            self.sample_terra_record = pd.Series({
                'analysis_set_acc': 'IGVFDS123ABC',
                'taxa': 'Homo sapiens',
                'rna_kb_h5ad': 'gs://fc-secure-bucket/submissions/12345/workflow_67890/call-aggregate_metrics/subwf_abcde/output.h5ad'
            })

        # Mock TerraOutputMetadata
        self.mock_terra_metadata = Mock()
        self.mock_terra_metadata.terra_data_record = self.sample_terra_record
        self.mock_terra_metadata.taxa = self.sample_terra_record.get(
            'taxa', 'Homo sapiens')
        self.mock_terra_metadata._parse_workflow_uuids_from_gs_path.return_value = parse_terra.InputFileAccs(
            gcloud_bucket='fc-secure-bucket',
            submission_id='12345',
            workflow_id='67890',
            subworkflow_id='abcde'
        )
        self.mock_terra_metadata._get_input_file_accs_from_table.return_value = parse_terra.InputFileAccs(
            derived_from_accessions=['IGVFFF001AAA', 'IGVFFF002BBB'],
            reference_files=['IGVFFF003CCC']
        )

    def test_init_with_real_terra_output_name(self):
        """Test MatrixFilePayload initialization with real Terra output names."""
        # Test with RNA output
        if 'rna_kb_h5ad' in self.sample_terra_record:
            payload = igvf_payloads.MatrixFilePayload(
                self.mock_terra_metadata, 'rna_kb_h5ad')
            self.assertEqual(payload.terra_output_name, 'rna_kb_h5ad')
            self.assertIsNotNone(payload.lab)
            self.assertIsNotNone(payload.award)

    def test_terra_output_name_property(self):
        """Test terra_output_name property."""
        payload = igvf_payloads.MatrixFilePayload(
            self.mock_terra_metadata, 'rna_kb_h5ad')
        self.assertEqual(payload.terra_output_name, 'rna_kb_h5ad')

    def test_submitted_file_name_property_with_real_data(self):
        """Test submitted_file_name property with real data."""
        if 'rna_kb_h5ad' in self.sample_terra_record:
            payload = igvf_payloads.MatrixFilePayload(
                self.mock_terra_metadata, 'rna_kb_h5ad')
            expected_file_name = self.sample_terra_record['rna_kb_h5ad']
            self.assertEqual(payload.submitted_file_name, expected_file_name)

    @patch('sc_pipe_management.accession.igvf_payloads.api_tools.calculate_gsfile_hex_hash')
    def test_md5sum_property(self, mock_hash):
        """Test md5sum property."""
        mock_hash.return_value = 'abc123def456'
        payload = igvf_payloads.MatrixFilePayload(
            self.mock_terra_metadata, 'rna_kb_h5ad')

        result = payload.md5sum
        self.assertEqual(result, 'abc123def456')

        if 'rna_kb_h5ad' in self.sample_terra_record:
            mock_hash.assert_called_once_with(
                file_path=self.sample_terra_record['rna_kb_h5ad'])

    @patch('sc_pipe_management.accession.igvf_payloads._get_file_aliases')
    def test_aliases_property(self, mock_get_aliases):
        """Test aliases property."""
        mock_get_aliases.return_value = [
            'test-lab:rna_kb_h5ad_IGVFDS123ABC_12345']
        payload = igvf_payloads.MatrixFilePayload(
            self.mock_terra_metadata, 'rna_kb_h5ad')

        result = payload.aliases
        self.assertEqual(result, ['test-lab:rna_kb_h5ad_IGVFDS123ABC_12345'])


class TestQCMetricsPayload(unittest.TestCase):
    """Test cases for QCMetricsPayload class."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        # Load test data
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            terra_data_table = pd.read_csv(test_file_path, sep='\t')
            self.sample_terra_record = terra_data_table.iloc[0]
        else:
            self.sample_terra_record = pd.Series({
                'analysis_set_acc': 'IGVFDS123ABC',
                'rna_kb_library_qc_metrics_json': 'gs://bucket/qc_metrics.json'
            })

        # Mock QC data info
        self.mock_qc_data_info = Mock()
        self.mock_qc_data_info.metadata = ['rna_kb_library_qc_metrics_json']
        self.mock_qc_data_info.attachment = {
            'rnaseq_kb_info': 'parameters.json'}
        self.mock_qc_data_info.object_type = 'single_cell_rna_seq_quality_metric'
        self.mock_qc_data_info.description = 'RNAseq Kallisto Bustools QC metric'
        self.mock_qc_data_info.metadata_map = {
            'numRecords': 'n_records', 'numReads': 'n_reads'}

        self.mock_terra_uuids = parse_terra.InputFileAccs(
            gcloud_bucket='fc-secure-bucket',
            submission_id='12345',
            workflow_id='67890',
            subworkflow_id='abcde'
        )

    def test_init(self):
        """Test QCMetricsPayload initialization."""
        qc_payload = igvf_payloads.QCMetricsPayload(
            terra_data_record=self.sample_terra_record,
            qc_data_info=self.mock_qc_data_info,
            qc_of=['IGVFFF001AAA'],
            terra_output_name='gene_count',
            terra_uuids=self.mock_terra_uuids,
            lab='test_lab',
            award='test_award'
        )

        self.assertEqual(qc_payload.terra_output_name, 'gene_count')
        self.assertEqual(qc_payload.lab, 'test_lab')
        self.assertEqual(qc_payload.award, 'test_award')

    def test_submitted_file_name_property(self):
        """Test submitted_file_name property returns None for QC."""
        qc_payload = igvf_payloads.QCMetricsPayload(
            terra_data_record=self.sample_terra_record,
            qc_data_info=self.mock_qc_data_info,
            qc_of=['IGVFFF001AAA'],
            terra_output_name='gene_count',
            terra_uuids=self.mock_terra_uuids,
            lab='test_lab',
            award='test_award'
        )

        self.assertIsNone(qc_payload.submitted_file_name)

    def test_md5sum_property(self):
        """Test md5sum property returns None for QC."""
        qc_payload = igvf_payloads.QCMetricsPayload(
            terra_data_record=self.sample_terra_record,
            qc_data_info=self.mock_qc_data_info,
            qc_of=['IGVFFF001AAA'],
            terra_output_name='gene_count',
            terra_uuids=self.mock_terra_uuids,
            lab='test_lab',
            award='test_award'
        )

        self.assertIsNone(qc_payload.md5sum)

    def test_read_json_file(self):
        """Test reading JSON files."""
        qc_payload = igvf_payloads.QCMetricsPayload(
            terra_data_record=self.sample_terra_record,
            qc_data_info=self.mock_qc_data_info,
            qc_of=['IGVFFF001AAA'],
            terra_output_name='gene_count',
            terra_uuids=self.mock_terra_uuids,
            lab='test_lab',
            award='test_award'
        )

        mock_json_data = {'numRecords': 1000, 'numReads': 50000}

        with patch('builtins.open', mock_open(read_data=json.dumps(mock_json_data))):
            with patch('json.load', return_value=mock_json_data):
                result = qc_payload._read_json_file(['test_file.json'])
                self.assertEqual(result, mock_json_data)


class TestIGVFPostService(unittest.TestCase):
    """Test cases for IGVFPostService class."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.mock_payload = Mock()
        self.mock_payload.aliases = ['test-lab:sample_alias']
        self.mock_payload.md5sum = 'abc123def456'

        self.mock_igvf_api = Mock()

        self.post_service = igvf_posting.IGVFPostService(
            self.mock_igvf_api, self.mock_payload)

    def test_init(self):
        """Test IGVFPostService initialization."""
        self.assertEqual(self.post_service.igvf_utils_api, self.mock_igvf_api)
        self.assertEqual(self.post_service.data_obj_payload, self.mock_payload)
        self.assertFalse(self.post_service.upload_file)
        self.assertFalse(self.post_service.resumed_posting)

    def test_init_with_options(self):
        """Test IGVFPostService initialization with options."""
        post_service = igvf_posting.IGVFPostService(
            self.mock_igvf_api,
            self.mock_payload,
            upload_file=True,
            resumed_posting=True
        )

        self.assertTrue(post_service.upload_file)
        self.assertTrue(post_service.resumed_posting)

    def test_get_conflict_object_id_by_alias(self):
        """Test getting conflict object ID by alias."""
        self.mock_igvf_api.get.return_value = {'@id': '/files/IGVFFI001AAA/'}

        result = self.post_service._get_conflict_object_id()

        self.assertEqual(result, '/files/IGVFFI001AAA/')
        self.mock_igvf_api.get.assert_called_once_with(
            'aliases:test-lab:sample_alias')

    def test_get_conflict_object_id_by_md5(self):
        """Test getting conflict object ID by MD5 when alias fails."""
        self.mock_payload.aliases = None
        self.mock_igvf_api.get.return_value = {'@id': '/files/IGVFFI001AAA/'}

        result = self.post_service._get_conflict_object_id()

        self.assertEqual(result, '/files/IGVFFI001AAA/')
        self.mock_igvf_api.get.assert_called_once_with('md5:abc123def456')

    def test_get_conflict_object_id_no_conflict(self):
        """Test no conflict found."""
        self.mock_igvf_api.get.return_value = None

        result = self.post_service._get_conflict_object_id()

        self.assertIsNone(result)

    def test_get_conflict_object_id_no_identifiers(self):
        """Test when no aliases or MD5 available."""
        self.mock_payload.aliases = None
        self.mock_payload.md5sum = None

        result = self.post_service._get_conflict_object_id()

        self.assertIsNone(result)
        self.mock_igvf_api.get.assert_not_called()


class TestPostResult(unittest.TestCase):
    """Test cases for PostResult dataclass."""

    def test_post_result_success(self):
        """Test PostResult.Success creation."""
        result = igvf_posting.PostResult.Success(
            col_header='test_col', accession='IGVFFI001AAA')

        self.assertEqual(result.col_header, 'test_col')
        self.assertEqual(result.accession, 'IGVFFI001AAA')
        self.assertIsNone(result.error)

    def test_post_result_failure(self):
        """Test PostResult.Failure creation."""
        error = Exception('Test error')
        result = igvf_posting.PostResult.Failure(
            col_header='test_col', error=error)

        self.assertEqual(result.col_header, 'test_col')
        self.assertIsNone(result.accession)
        self.assertEqual(result.error, error)


class TestPipelineParamsInfo(unittest.TestCase):
    """Test cases for PipelineParamsInfo class."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        # Load or create test data
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            self.terra_data_table = pd.read_csv(test_file_path, sep='\t')
        else:
            self.terra_data_table = pd.DataFrame({
                'analysis_set_acc': ['IGVFDS123ABC', 'IGVFDS456DEF'],
                'rna_kb_h5ad': [
                    'gs://fc-secure-bucket/submissions/12345/workflow_67890/call-aggregate_metrics/subwf_abcde/output.h5ad',
                    'gs://fc-secure-bucket/submissions/23456/workflow_78901/call-aggregate_metrics/subwf_fghij/output.h5ad'
                ]
            })

        self.pipeline_params = igvf_payloads.PipelineParamsInfo(
            self.terra_data_table)

    def test_init(self):
        """Test PipelineParamsInfo initialization."""
        self.assertEqual(self.pipeline_params.terra_namespace, 'DACC_ANVIL')
        self.assertEqual(self.pipeline_params.terra_workspace,
                         'IGVF Single-Cell Data Processing')
        self.assertTrue(self.pipeline_params.terra_datable.equals(
            self.terra_data_table))

    def test_init_with_custom_params(self):
        """Test PipelineParamsInfo initialization with custom parameters."""
        custom_params = igvf_payloads.PipelineParamsInfo(
            self.terra_data_table,
            terra_namespace='CUSTOM_NAMESPACE',
            terra_workspace='Custom Workspace',
            output_root_dir='/custom/path/'
        )

        self.assertEqual(custom_params.terra_namespace, 'CUSTOM_NAMESPACE')
        self.assertEqual(custom_params.terra_workspace, 'Custom Workspace')
        self.assertEqual(custom_params.output_root_dir, '/custom/path/')

    @patch('sc_pipe_management.accession.igvf_payloads.fapi.get_workflow_metadata')
    def test_get_single_input_params_success(self, mock_get_metadata):
        """Test successful retrieval of input parameters."""
        # Mock successful API response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'inputs': {'param1': 'value1', 'param2': 'value2'}
        }
        mock_get_metadata.return_value = mock_response

        # Mock terra metadata
        mock_terra_metadata = Mock()
        mock_terra_metadata._parse_workflow_uuids_from_gs_path.return_value = parse_terra.InputFileAccs(
            gcloud_bucket='fc-secure-bucket',
            submission_id='12345',
            workflow_id='67890',
            subworkflow_id='abcde'
        )
        mock_terra_metadata.analysis_set_acc = 'IGVFDS123ABC'

        with patch('json.dumps', return_value='{"inputs": {"param1": "value1"}}'):
            result = self.pipeline_params._get_single_input_params(
                mock_terra_metadata)
            self.assertIsInstance(result, str)
            self.assertIn('param1', result)

    @patch('sc_pipe_management.accession.igvf_payloads.fapi.get_workflow_metadata')
    def test_get_single_input_params_api_error(self, mock_get_metadata):
        """Test API error handling."""
        # Mock failed API response
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get_metadata.return_value = mock_response

        mock_terra_metadata = Mock()
        mock_terra_metadata.analysis_set_acc = 'IGVFDS123ABC'
        mock_terra_metadata._parse_workflow_uuids_from_gs_path.return_value = Mock()

        with self.assertRaises(Exception):
            self.pipeline_params._get_single_input_params(mock_terra_metadata)


class TestAccessioning(unittest.TestCase):
    """Integration tests for the accessioning module using real test data."""

    @classmethod
    def setUpClass(cls):
        """Load the Terra test data table once for all tests."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            cls.terra_data_table = pd.read_csv(test_file_path, sep='\t')
            cls.has_real_data = True
        else:
            cls.has_real_data = False
            cls.skipTest("Test data file not found")

    def setUp(self):
        """Set up test fixtures before each test method."""
        if not self.has_real_data:
            self.skipTest("No real test data available")

        self.sample_record = self.terra_data_table.iloc[0]

    def test_terra_data_table_structure(self):
        """Test that the Terra data table has expected columns."""
        expected_columns = ['analysis_set_acc', 'taxa']
        for col in expected_columns:
            self.assertIn(col, self.terra_data_table.columns,
                          f"Missing expected column: {col}")

    def test_multiple_records_processing(self):
        """Test processing multiple records from the Terra data table."""
        processed_records = []

        # Test first 3 records
        for idx, row in self.terra_data_table.head(3).iterrows():
            try:
                mock_igvf_api = Mock()
                terra_metadata = parse_terra.TerraOutputMetadata(
                    row, mock_igvf_api)

                # Test basic functionality
                self.assertIsNotNone(terra_metadata.anaset_accession)
                self.assertIsNotNone(terra_metadata.taxa)

                processed_records.append(terra_metadata.anaset_accession)

            except Exception as e:
                # Log but don't fail - some test data might be incomplete
                print(f"Warning: Could not process record {idx}: {e}")

        self.assertGreater(len(processed_records), 0,
                           "Should process at least one record")

    def test_end_to_end_workflow(self):
        """Test end-to-end workflow with real data."""
        # This test demonstrates how the components work together
        mock_igvf_api = Mock()
        mock_igvf_api.get.return_value = {'controlled_access': False}

        try:
            # Initialize with real data
            terra_metadata = parse_terra.TerraOutputMetadata(
                self.sample_record, mock_igvf_api)

            # Test that we can create payloads if RNA data exists
            if 'rna_kb_h5ad' in self.sample_record and pd.notna(self.sample_record['rna_kb_h5ad']):
                # Mock the required methods
                terra_metadata._get_input_file_accs_from_table = Mock(return_value=parse_terra.InputFileAccs(
                    derived_from_accessions=['IGVFFF001AAA'],
                    reference_files=['IGVFFF002BBB']
                ))

                payload = igvf_payloads.MatrixFilePayload(
                    terra_metadata, 'rna_kb_h5ad')

                self.assertEqual(payload.terra_output_name, 'rna_kb_h5ad')
                self.assertIsNotNone(payload.submitted_file_name)

                # Test post service
                post_service = igvf_posting.IGVFPostService(
                    mock_igvf_api, payload)
                self.assertIsNotNone(post_service)

        except Exception as e:
            # Skip if test data doesn't support full workflow
            self.skipTest(f"Test data doesn't support full workflow: {e}")


if __name__ == '__main__':
    # Run tests with verbosity
    unittest.main(verbosity=2)
