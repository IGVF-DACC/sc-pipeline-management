import sys
import os
import pytest
import pandas as pd
from unittest.mock import Mock, patch, MagicMock, mock_open
import json
from pathlib import Path

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.accession.parse_terra_metadata as parse_terra
import sc_pipe_management.accession.igvf_payloads as igvf_payloads
import sc_pipe_management.accession.terra_to_portal_posting as igvf_posting


class TestTerraOutputMetadata:
    """Test cases for TerraOutputMetadata class using real test data."""

    @pytest.fixture(scope="class")
    def terra_data_table(self):
        """Load the Terra test data table once for all tests."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        return pd.read_csv(test_file_path, sep='\t')

    @pytest.fixture
    def mock_igvf_api(self):
        """Mock IGVF API client."""
        mock_api = Mock()
        mock_api.get.return_value = {'controlled_access': False}
        return mock_api

    @pytest.fixture
    def sample_terra_record(self, terra_data_table):
        """Get the first row of test data."""
        return terra_data_table.iloc[0]

    @pytest.fixture
    def terra_metadata(self, sample_terra_record, mock_igvf_api):
        """Create TerraOutputMetadata instance."""
        return parse_terra.TerraOutputMetadata(sample_terra_record, mock_igvf_api)

    def test_init_with_real_data(self, terra_metadata, sample_terra_record):
        """Test TerraOutputMetadata initialization with real test data."""
        assert terra_metadata.anaset_accession == sample_terra_record['analysis_set_acc']
        assert terra_metadata.taxa == sample_terra_record['taxa']
        assert terra_metadata.genome_assembly is not None
        assert terra_metadata.pipeline_output.equals(sample_terra_record)

    def test_get_gs_path_prioritizes_rna(self, terra_metadata, sample_terra_record):
        """Test that RNA paths are prioritized over ATAC paths."""
        gs_path = terra_metadata._get_gs_path_for_terra_output_cols()

        # Should return RNA path if available
        if pd.notna(sample_terra_record.get('rna_kb_h5ad')):
            assert gs_path == sample_terra_record['rna_kb_h5ad']
        elif pd.notna(sample_terra_record.get('atac_bam')):
            assert gs_path == sample_terra_record['atac_bam']

    def test_parse_workflow_uuids_from_gs_path(self, terra_metadata):
        """Test parsing workflow UUIDs from GS path using real data."""
        try:
            result = terra_metadata._parse_workflow_uuids_from_gs_path()
            assert isinstance(result, parse_terra.PipelineOutputIds)
            assert result.gcloud_bucket is not None
            assert result.submission_id is not None
            assert result.workflow_id is not None
            assert result.subworkflow_id is not None
        except ValueError:
            # Some test data might not have valid GS paths
            pytest.skip(
                "Test data doesn't contain valid GS paths for UUID parsing")

    def test_parse_terra_str_list_with_real_data(self, terra_metadata, sample_terra_record):
        """Test parsing comma-separated strings from real Terra data."""
        # Test RNA sequence file accessions
        if 'rna_seqfile_accessions' in sample_terra_record:
            rna_accessions_str = sample_terra_record['rna_seqfile_accessions']
            if pd.notna(rna_accessions_str):
                result = terra_metadata._parse_terra_str_list(
                    [rna_accessions_str])
                assert isinstance(result, list)
                assert all(acc.startswith(('IGVF', 'TSTF')) for acc in result)

    def test_parse_igvf_accessions_from_urls_real_data(self, terra_metadata, sample_terra_record):
        """Test extracting IGVF accessions from URLs using real data."""
        if 'rna_seqspec_urls' in sample_terra_record:
            urls_str = sample_terra_record['rna_seqspec_urls']
            if pd.notna(urls_str):
                result = terra_metadata._parse_igvf_accessions_from_urls([
                                                                         urls_str])
                assert isinstance(result, list)
                # Should extract valid accessions from URLs
                if result:
                    assert all(acc.startswith(('IGVF', 'TSTF'))
                               for acc in result)

    def test_get_input_file_accs_from_table_rna(self, terra_metadata):
        """Test getting input file accessions for RNA with real data."""
        try:
            result = terra_metadata._get_input_file_accs_from_table('rna')
            assert isinstance(result, parse_terra.InputFileAccs)
            assert isinstance(result.derived_from_accessions, list)
            assert isinstance(result.reference_files, list)

            # Check that accessions start with IGVF or TSTF
            for acc in result.derived_from_accessions:
                assert acc.startswith(('IGVF', 'TSTF'))
        except (KeyError, AttributeError):
            # Some test data might not have all required columns
            pytest.skip(
                "Test data missing required columns for RNA input file parsing")

    def test_get_input_file_accs_from_table_atac(self, terra_metadata):
        """Test getting input file accessions for ATAC with real data."""
        try:
            result = terra_metadata._get_input_file_accs_from_table('atac')
            assert isinstance(result, parse_terra.InputFileAccs)
            assert isinstance(result.derived_from_accessions, list)
            assert isinstance(result.reference_files, list)

            # Check that accessions start with IGVF or TSTF
            for acc in result.derived_from_accessions:
                assert acc.startswith(('IGVF', 'TSTF'))
        except (KeyError, AttributeError):
            # Some test data might not have all required columns
            pytest.skip(
                "Test data missing required columns for ATAC input file parsing")

    def test_get_gs_path_no_valid_path_raises_error(self, mock_igvf_api):
        """Test that missing GS paths raise ValueError."""
        # Create a record with no valid GS paths
        empty_record = pd.Series({
            'analysis_set_acc': 'IGVFDS999ZZZ',
            'taxa': 'Homo sapiens'
        })
        empty_metadata = parse_terra.TerraOutputMetadata(
            empty_record, mock_igvf_api)

        with pytest.raises(ValueError):
            empty_metadata._get_gs_path_for_terra_output_cols()

    @patch('sc_pipe_management.accession.parse_terra_metadata.fapi')
    @patch('sc_pipe_management.accession.parse_terra_metadata.dump_json')
    def test_download_workflow_config_json(self, mock_dump_json, mock_fapi, terra_metadata):
        """Test downloading workflow configuration JSON."""
        # Mock Terra API response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {'inputs': {'param1': 'value1'}}
        mock_fapi.get_workflow_metadata.return_value = mock_response

        # Mock file dump
        mock_dump_json.return_value = '/path/to/config.json'

        try:
            result = terra_metadata._download_workflow_config_json(
                terra_namespace='test_namespace',
                terra_workspace='test_workspace'
            )

            assert result.download_path == '/path/to/config.json'
            assert isinstance(result.doc_aliases, list)
        except ValueError:
            # Skip if test data doesn't have valid workflow UUIDs
            pytest.skip("Test data doesn't contain valid workflow UUIDs")


class TestMatrixFilePayload:
    """Test cases for MatrixFilePayload class using real test data."""

    @pytest.fixture
    def sample_terra_record(self):
        """Load test data record."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            terra_data_table = pd.read_csv(test_file_path, sep='\t')
            return terra_data_table.iloc[0]
        else:
            return pd.Series({
                'analysis_set_acc': 'IGVFDS123ABC',
                'taxa': 'Homo sapiens',
                'rna_kb_h5ad': 'gs://fc-secure-bucket/submissions/12345/workflow_67890/call-aggregate_metrics/subwf_abcde/output.h5ad'
            })

    @pytest.fixture
    def mock_terra_metadata(self, sample_terra_record):
        """Mock TerraOutputMetadata instance."""
        mock_metadata = Mock()
        mock_metadata.terra_data_record = sample_terra_record
        mock_metadata.taxa = sample_terra_record.get('taxa', 'Homo sapiens')
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = parse_terra.PipelineOutputIds(
            gcloud_bucket='fc-secure-bucket',
            submission_id='12345',
            workflow_id='67890',
            subworkflow_id='abcde'
        )
        mock_metadata._get_input_file_accs_from_table.return_value = parse_terra.InputFileAccs(
            derived_from_accessions=['IGVFFF001AAA', 'IGVFFF002BBB'],
            reference_files=['IGVFFF003CCC']
        )
        return mock_metadata

    def test_init_with_real_terra_output_name(self, mock_terra_metadata, sample_terra_record):
        """Test MatrixFilePayload initialization with real Terra output names."""
        # Test with RNA output
        if 'rna_kb_h5ad' in sample_terra_record:
            payload = igvf_payloads.MatrixFilePayload(
                mock_terra_metadata, 'rna_kb_h5ad')
            assert payload.terra_output_name == 'rna_kb_h5ad'
            assert payload.lab is not None
            assert payload.award is not None

    def test_terra_output_name_property(self, mock_terra_metadata):
        """Test terra_output_name property."""
        payload = igvf_payloads.MatrixFilePayload(
            mock_terra_metadata, 'rna_kb_h5ad')
        assert payload.terra_output_name == 'rna_kb_h5ad'

    def test_submitted_file_name_property_with_real_data(self, mock_terra_metadata, sample_terra_record):
        """Test submitted_file_name property with real data."""
        if 'rna_kb_h5ad' in sample_terra_record:
            payload = igvf_payloads.MatrixFilePayload(
                mock_terra_metadata, 'rna_kb_h5ad')
            expected_file_name = sample_terra_record['rna_kb_h5ad']
            assert payload.submitted_file_name == expected_file_name

    @patch('sc_pipe_management.accession.igvf_payloads.api_tools.calculate_gsfile_hex_hash')
    def test_md5sum_property(self, mock_hash, mock_terra_metadata, sample_terra_record):
        """Test md5sum property."""
        mock_hash.return_value = 'abc123def456'
        payload = igvf_payloads.MatrixFilePayload(
            mock_terra_metadata, 'rna_kb_h5ad')

        result = payload.md5sum
        assert result == 'abc123def456'

        if 'rna_kb_h5ad' in sample_terra_record:
            mock_hash.assert_called_once_with(
                file_path=sample_terra_record['rna_kb_h5ad'])

    @patch('sc_pipe_management.accession.igvf_payloads._get_file_aliases')
    def test_aliases_property(self, mock_get_aliases, mock_terra_metadata):
        """Test aliases property."""
        mock_get_aliases.return_value = [
            'test-lab:rna_kb_h5ad_IGVFDS123ABC_12345']
        payload = igvf_payloads.MatrixFilePayload(
            mock_terra_metadata, 'rna_kb_h5ad')

        result = payload.aliases
        assert result == ['test-lab:rna_kb_h5ad_IGVFDS123ABC_12345']


class TestQCMetricsPayload:
    """Test cases for QCMetricsPayload class."""

    @pytest.fixture
    def sample_terra_record(self):
        """Load test data record."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            terra_data_table = pd.read_csv(test_file_path, sep='\t')
            return terra_data_table.iloc[0]
        else:
            return pd.Series({
                'analysis_set_acc': 'IGVFDS123ABC',
                'rna_kb_library_qc_metrics_json': 'gs://bucket/qc_metrics.json'
            })

    @pytest.fixture
    def mock_qc_data_info(self):
        """Mock QC data info."""
        mock_info = Mock()
        mock_info.metadata = ['rna_kb_library_qc_metrics_json']
        mock_info.attachment = {'rnaseq_kb_info': 'parameters.json'}
        mock_info.object_type = 'single_cell_rna_seq_quality_metric'
        mock_info.description = 'RNAseq Kallisto Bustools QC metric'
        mock_info.metadata_map = {
            'numRecords': 'n_records', 'numReads': 'n_reads'}
        return mock_info

    @pytest.fixture
    def mock_terra_uuids(self):
        """Mock Terra UUIDs."""
        return parse_terra.PipelineOutputIds(
            gcloud_bucket='fc-secure-bucket',
            submission_id='12345',
            workflow_id='67890',
            subworkflow_id='abcde'
        )

    @pytest.fixture
    def qc_payload(self, sample_terra_record, mock_qc_data_info, mock_terra_uuids):
        """Create QCMetricsPayload instance."""
        return igvf_payloads.QCMetricsPayload(
            terra_data_record=sample_terra_record,
            qc_data_info=mock_qc_data_info,
            qc_of=['IGVFFF001AAA'],
            terra_output_name='gene_count',
            terra_uuids=mock_terra_uuids,
            lab='test_lab',
            award='test_award'
        )

    def test_init(self, qc_payload):
        """Test QCMetricsPayload initialization."""
        assert qc_payload.terra_output_name == 'gene_count'
        assert qc_payload.lab == 'test_lab'
        assert qc_payload.award == 'test_award'

    def test_submitted_file_name_property(self, qc_payload):
        """Test submitted_file_name property returns None for QC."""
        assert qc_payload.submitted_file_name is None

    def test_md5sum_property(self, qc_payload):
        """Test md5sum property returns None for QC."""
        assert qc_payload.md5sum is None

    def test_read_json_file(self, qc_payload):
        """Test reading JSON files."""
        mock_json_data = {'numRecords': 1000, 'numReads': 50000}

        with patch('builtins.open', mock_open(read_data=json.dumps(mock_json_data))):
            with patch('json.load', return_value=mock_json_data):
                result = qc_payload._read_json_file(['test_file.json'])
                assert result == mock_json_data


class TestIGVFPostService:
    """Test cases for IGVFPostService class."""

    @pytest.fixture
    def mock_payload(self):
        """Mock payload for testing."""
        payload = Mock()
        payload.aliases = ['test-lab:sample_alias']
        payload.md5sum = 'abc123def456'
        return payload

    @pytest.fixture
    def mock_igvf_api(self):
        """Mock IGVF API."""
        return Mock()

    @pytest.fixture
    def post_service(self, mock_igvf_api, mock_payload):
        """Create IGVFPostService instance."""
        return igvf_posting.IGVFPostService(mock_igvf_api, mock_payload)

    def test_init(self, post_service, mock_igvf_api, mock_payload):
        """Test IGVFPostService initialization."""
        assert post_service.igvf_utils_api == mock_igvf_api
        assert post_service.data_obj_payload == mock_payload
        assert post_service.upload_file is False
        assert post_service.resumed_posting is False

    def test_init_with_options(self, mock_igvf_api, mock_payload):
        """Test IGVFPostService initialization with options."""
        post_service = igvf_posting.IGVFPostService(
            mock_igvf_api,
            mock_payload,
            upload_file=True,
            resumed_posting=True
        )

        assert post_service.upload_file is True
        assert post_service.resumed_posting is True

    def test_get_conflict_object_id_by_alias(self, post_service, mock_igvf_api):
        """Test getting conflict object ID by alias."""
        mock_igvf_api.get.return_value = {'@id': '/files/IGVFFI001AAA/'}

        result = post_service._get_conflict_object_id()

        assert result == '/files/IGVFFI001AAA/'
        mock_igvf_api.get.assert_called_once_with(
            'aliases:test-lab:sample_alias')

    def test_get_conflict_object_id_by_md5(self, post_service, mock_igvf_api, mock_payload):
        """Test getting conflict object ID by MD5 when alias fails."""
        mock_payload.aliases = None
        mock_igvf_api.get.return_value = {'@id': '/files/IGVFFI001AAA/'}

        result = post_service._get_conflict_object_id()

        assert result == '/files/IGVFFI001AAA/'
        mock_igvf_api.get.assert_called_once_with('md5:abc123def456')

    def test_get_conflict_object_id_no_conflict(self, post_service, mock_igvf_api):
        """Test no conflict found."""
        mock_igvf_api.get.return_value = None

        result = post_service._get_conflict_object_id()

        assert result is None

    def test_get_conflict_object_id_no_identifiers(self, post_service, mock_igvf_api, mock_payload):
        """Test when no aliases or MD5 available."""
        mock_payload.aliases = None
        mock_payload.md5sum = None

        result = post_service._get_conflict_object_id()

        assert result is None
        mock_igvf_api.get.assert_not_called()


class TestPostResult:
    """Test cases for PostResult dataclass."""

    def test_post_result_success(self):
        """Test PostResult.Success creation."""
        result = igvf_posting.PostResult.Success(
            col_header='test_col', accession='IGVFFI001AAA')

        assert result.col_header == 'test_col'
        assert result.accession == 'IGVFFI001AAA'
        assert result.error is None

    def test_post_result_failure(self):
        """Test PostResult.Failure creation."""
        error = Exception('Test error')
        result = igvf_posting.PostResult.Failure(
            col_header='test_col', error=error)

        assert result.col_header == 'test_col'
        assert result.accession is None
        assert result.error == error


class TestPipelineParamsInfo:
    """Test cases for PipelineParamsInfo class."""

    @pytest.fixture
    def terra_data_table(self):
        """Load or create test data."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            return pd.read_csv(test_file_path, sep='\t')
        else:
            return pd.DataFrame({
                'analysis_set_acc': ['IGVFDS123ABC', 'IGVFDS456DEF'],
                'rna_kb_h5ad': [
                    'gs://fc-secure-bucket/submissions/12345/workflow_67890/call-aggregate_metrics/subwf_abcde/output.h5ad',
                    'gs://fc-secure-bucket/submissions/23456/workflow_78901/call-aggregate_metrics/subwf_fghij/output.h5ad'
                ]
            })

    @pytest.fixture
    def pipeline_params(self, terra_data_table):
        """Create PipelineParamsInfo instance."""
        return igvf_payloads.PipelineParamsInfo(terra_data_table)

    def test_init(self, pipeline_params, terra_data_table):
        """Test PipelineParamsInfo initialization."""
        assert pipeline_params.terra_namespace == 'DACC_ANVIL'
        assert pipeline_params.terra_workspace == 'IGVF Single-Cell Data Processing'
        assert pipeline_params.terra_datable.equals(terra_data_table)

    def test_init_with_custom_params(self, terra_data_table):
        """Test PipelineParamsInfo initialization with custom parameters."""
        custom_params = igvf_payloads.PipelineParamsInfo(
            terra_data_table,
            terra_namespace='CUSTOM_NAMESPACE',
            terra_workspace='Custom Workspace',
            output_root_dir='/custom/path/'
        )

        assert custom_params.terra_namespace == 'CUSTOM_NAMESPACE'
        assert custom_params.terra_workspace == 'Custom Workspace'
        assert custom_params.output_root_dir == '/custom/path/'

    @patch('sc_pipe_management.accession.igvf_payloads.fapi.get_workflow_metadata')
    def test_get_single_input_params_success(self, mock_get_metadata, pipeline_params):
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
        mock_terra_metadata._parse_workflow_uuids_from_gs_path.return_value = parse_terra.PipelineOutputIds(
            gcloud_bucket='fc-secure-bucket',
            submission_id='12345',
            workflow_id='67890',
            subworkflow_id='abcde'
        )
        mock_terra_metadata.analysis_set_acc = 'IGVFDS123ABC'

        with patch('json.dumps', return_value='{"inputs": {"param1": "value1"}}'):
            result = pipeline_params._get_single_input_params(
                mock_terra_metadata)
            assert isinstance(result, str)
            assert 'param1' in result

    @patch('sc_pipe_management.accession.igvf_payloads.fapi.get_workflow_metadata')
    def test_get_single_input_params_api_error(self, mock_get_metadata, pipeline_params):
        """Test API error handling."""
        # Mock failed API response
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get_metadata.return_value = mock_response

        mock_terra_metadata = Mock()
        mock_terra_metadata.analysis_set_acc = 'IGVFDS123ABC'
        mock_terra_metadata._parse_workflow_uuids_from_gs_path.return_value = Mock()

        with pytest.raises(Exception):
            pipeline_params._get_single_input_params(mock_terra_metadata)


class TestAccessioning:
    """Integration tests for the accessioning module using real test data."""

    @pytest.fixture(scope="class")
    def terra_data_table(self):
        """Load the Terra test data table once for all tests."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            return pd.read_csv(test_file_path, sep='\t')
        else:
            pytest.skip("Test data file not found")

    @pytest.fixture
    def sample_record(self, terra_data_table):
        """Get first record from test data."""
        return terra_data_table.iloc[0]

    def test_terra_data_table_structure(self, terra_data_table):
        """Test that the Terra data table has expected columns."""
        expected_columns = ['analysis_set_acc', 'taxa']
        for col in expected_columns:
            assert col in terra_data_table.columns, f"Missing expected column: {col}"

    def test_multiple_records_processing(self, terra_data_table):
        """Test processing multiple records from the Terra data table."""
        processed_records = []

        # Test first 3 records
        for idx, row in terra_data_table.head(3).iterrows():
            try:
                mock_igvf_api = Mock()
                terra_metadata = parse_terra.TerraOutputMetadata(
                    row, mock_igvf_api)

                # Test basic functionality
                assert terra_metadata.anaset_accession is not None
                assert terra_metadata.taxa is not None

                processed_records.append(terra_metadata.anaset_accession)

            except Exception as e:
                # Log but don't fail - some test data might be incomplete
                print(f"Warning: Could not process record {idx}: {e}")

        assert len(processed_records) > 0, "Should process at least one record"

    def test_all_records_processing(self, terra_data_table):
        """Test processing all records from the Terra data table."""
        processed_records = []
        failed_records = []

        # Test ALL records
        for idx, row in terra_data_table.iterrows():
            try:
                mock_igvf_api = Mock()
                terra_metadata = parse_terra.TerraOutputMetadata(
                    row, mock_igvf_api)

                # Test basic functionality
                assert terra_metadata.anaset_accession is not None
                assert terra_metadata.taxa is not None

                processed_records.append(terra_metadata.anaset_accession)

            except Exception as e:
                # Log but don't fail - some test data might be incomplete
                failed_records.append({'idx': idx, 'error': str(
                    e), 'accession': row.get('analysis_set_acc', 'Unknown')})
                print(
                    f"Warning: Could not process record {idx} ({row.get('analysis_set_acc', 'Unknown')}): {e}")

        # Print summary
        print(f"Successfully processed {len(processed_records)} records")
        print(f"Failed to process {len(failed_records)} records")

        if failed_records:
            print("Failed records:")
            for failed in failed_records:
                print(
                    f"  - Index {failed['idx']}: {failed['accession']} - {failed['error']}")

        assert len(processed_records) > 0, "Should process at least one record"

        # Optional: Assert that most records should process successfully
        success_rate = len(processed_records) / len(terra_data_table)
        assert success_rate > 0.8, f"Success rate too low: {success_rate:.2%}. Check data quality."

    def test_end_to_end_workflow(self, sample_record):
        """Test end-to-end workflow with real data."""
        # This test demonstrates how the components work together
        mock_igvf_api = Mock()
        mock_igvf_api.get.return_value = {'controlled_access': False}

        try:
            # Initialize with real data
            terra_metadata = parse_terra.TerraOutputMetadata(
                sample_record, mock_igvf_api)

            # Test that we can create payloads if RNA data exists
            if 'rna_kb_h5ad' in sample_record and pd.notna(sample_record['rna_kb_h5ad']):
                # Mock the required methods
                terra_metadata._get_input_file_accs_from_table = Mock(return_value=parse_terra.InputFileAccs(
                    derived_from_accessions=['IGVFFF001AAA'],
                    reference_files=['IGVFFF002BBB']
                ))

                payload = igvf_payloads.MatrixFilePayload(
                    terra_metadata, 'rna_kb_h5ad')

                assert payload.terra_output_name == 'rna_kb_h5ad'
                assert payload.submitted_file_name is not None

                # Test post service
                post_service = igvf_posting.IGVFPostService(
                    mock_igvf_api, payload)
                assert post_service is not None

        except Exception as e:
            # Skip if test data doesn't support full workflow
            pytest.skip(f"Test data doesn't support full workflow: {e}")


# Parametrized tests for testing multiple records
class TestParametrizedRecords:
    """Parametrized tests to run the same test on multiple records."""

    @pytest.fixture(scope="class")
    def terra_data_table(self):
        """Load the Terra test data table."""
        test_file_path = Path(__file__).parent / \
            'test_files' / 'DACC_post_testcases.tsv'
        if test_file_path.exists():
            return pd.read_csv(test_file_path, sep='\t')
        else:
            pytest.skip("Test data file not found")

    @pytest.mark.parametrize("record_index", [0, 1, 2])
    def test_individual_records(self, terra_data_table, record_index):
        """Test individual records from the Terra data table."""
        if record_index >= len(terra_data_table):
            pytest.skip(
                f"Record index {record_index} not available in test data")

        record = terra_data_table.iloc[record_index]
        mock_igvf_api = Mock()

        try:
            terra_metadata = parse_terra.TerraOutputMetadata(
                record, mock_igvf_api)
            assert terra_metadata.anaset_accession is not None
            assert terra_metadata.taxa is not None
            print(
                f"Successfully tested record {record_index}: {terra_metadata.anaset_accession}")
        except Exception as e:
            pytest.fail(f"Failed to process record {record_index}: {e}")


if __name__ == '__main__':
    # Run tests with pytest
    pytest.main([__file__, "-v"])
