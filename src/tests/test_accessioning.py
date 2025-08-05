import pytest
import pandas as pd
import json
import os
import tempfile
from unittest.mock import Mock, patch, MagicMock, mock_open
from pathlib import Path

# Import the modules to test
import sc_pipe_management.accession.igvf_payloads as igvf_payloads
import sc_pipe_management.accession.parse_terra_metadata as parse_terra


class TestUtilityFunctions:
    """Test utility functions in igvf_payloads module."""

    def test_dump_json(self):
        """Test _dump_json function."""
        test_json = {"key1": "value1", "key2": "value2"}
        analysis_set_acc = "IGVFDS123ABC"

        with tempfile.TemporaryDirectory() as temp_dir:
            result_path = igvf_payloads._dump_json(
                input_json=test_json,
                analysis_set_acc=analysis_set_acc,
                output_root_dir=temp_dir
            )

            # Check that file was created
            assert os.path.exists(result_path)

            # Check file content
            with open(result_path, 'r') as f:
                saved_json = json.load(f)
            assert saved_json == test_json

            # Check file path structure
            assert analysis_set_acc in result_path
            assert result_path.endswith('.json')

    def test_get_file_aliases(self):
        """Test _get_file_aliases function."""
        col_header = "rna_kb_h5ad"
        lab = "/labs/test-lab/"
        terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC'
        })
        terra_uuids = parse_terra.TerraJobUUIDs(
            gcloud_bucket='bucket',
            submission_id='sub123',
            workflow_id='wf456',
            subworkflow_id='subwf789'
        )

        result = igvf_payloads._get_file_aliases(
            col_header=col_header,
            lab=lab,
            terra_data_record=terra_data_record,
            terra_uuids=terra_uuids
        )

        expected = [
            'test-lab:IGVFDS123ABC_bucket_sub123_wf456_subwf789_rna_kb_h5ad_uniform-pipeline']
        assert result == expected

    @patch('subprocess.run')
    def test_download_qc_file_from_gcp_success(self, mock_subprocess):
        """Test successful GCP file download."""
        # Mock successful subprocess
        mock_subprocess.return_value.returncode = 0

        gs_file_path = "gs://bucket/path/file.json"

        with tempfile.TemporaryDirectory() as temp_dir:
            result = igvf_payloads._download_qc_file_from_gcp(
                gs_file_path=gs_file_path,
                downloaded_dir=temp_dir
            )

            expected_path = os.path.join(temp_dir, "file.json")
            assert result == expected_path

            # Check subprocess was called correctly
            mock_subprocess.assert_called_once_with(
                ['gcloud', 'storage', 'cp', gs_file_path, expected_path],
                capture_output=True
            )

    @patch('subprocess.run')
    def test_download_qc_file_from_gcp_failure(self, mock_subprocess):
        """Test failed GCP file download."""
        # Mock failed subprocess
        mock_subprocess.return_value.returncode = 1
        mock_subprocess.return_value.stderr.decode.return_value = "ERROR: File not found\n"

        gs_file_path = "gs://bucket/path/nonexistent.json"

        with tempfile.TemporaryDirectory() as temp_dir:
            with pytest.raises(Exception, match="GCP download error"):
                igvf_payloads._download_qc_file_from_gcp(
                    gs_file_path=gs_file_path,
                    downloaded_dir=temp_dir
                )

    def test_download_qc_file_invalid_path(self):
        """Test download with invalid GCP path."""
        invalid_path = "/local/path/file.json"

        with tempfile.TemporaryDirectory() as temp_dir:
            with pytest.raises(Exception, match="not a valid GCP file path"):
                igvf_payloads._download_qc_file_from_gcp(
                    gs_file_path=invalid_path,
                    downloaded_dir=temp_dir
                )


class TestMatrixFilePayload:
    """Test MatrixFilePayload class."""

    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance."""
        mock_metadata = Mock()
        mock_metadata.terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC',
            'rna_kb_h5ad': 'gs://bucket/file.h5ad'
        })
        mock_metadata.taxa = 'Homo sapiens'
        mock_metadata.anaset_accession = 'IGVFDS123ABC'

        # Mock UUIDs
        mock_uuids = Mock()
        mock_uuids.aliases.return_value = 'bucket_sub123_wf456_subwf789'
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        # Mock input file accessions
        mock_input_files = Mock()
        mock_input_files.get_derived_from.return_value = [
            'IGVFFF001AAA', 'IGVFFF002BBB']
        mock_input_files.reference_files = ['IGVFFF003CCC']
        mock_metadata._get_input_file_accs_from_table.return_value = mock_input_files

        return mock_metadata

    @pytest.fixture
    @patch('sc_pipe_management.accession.igvf_payloads.const.OUTPUT_SUBMITTER_INFO',
           {'lab': '/labs/test-lab/', 'award': '/awards/test-award/'})
    @patch('sc_pipe_management.accession.igvf_payloads.const.MATRIX_FILETYPES')
    @patch('sc_pipe_management.accession.igvf_payloads.const.GENOME_ASSEMBLY_INFO')
    def matrix_payload(self, mock_genome_info, mock_matrix_types, mock_terra_metadata):
        """Create MatrixFilePayload instance for testing."""
        # Mock file metadata
        mock_file_metadata = Mock()
        mock_file_metadata.assay_type = 'rna'
        mock_file_metadata.analysis_step_version = '/analysis-step-versions/test/'
        mock_file_metadata.content_type = 'intermediate file'
        mock_file_metadata.file_format = 'h5ad'
        mock_file_metadata.description = 'Test matrix file'
        mock_file_metadata.file_format_specifications = [
            '/file-format-specifications/test/']

        mock_matrix_types.__getitem__.return_value = mock_file_metadata
        mock_genome_info.get.return_value = 'GRCh38'

        return igvf_payloads.MatrixFilePayload(mock_terra_metadata, 'rna_kb_h5ad')

    def test_init(self, matrix_payload, mock_terra_metadata):
        """Test MatrixFilePayload initialization."""
        assert matrix_payload.terra_output_name == 'rna_kb_h5ad'
        assert matrix_payload.lab == '/labs/test-lab/'
        assert matrix_payload.award == '/awards/test-award/'
        assert matrix_payload.terra_metadata == mock_terra_metadata

    def test_aliases_property(self, matrix_payload):
        """Test aliases property."""
        with patch('sc_pipe_management.accession.igvf_payloads._get_file_aliases') as mock_aliases:
            mock_aliases.return_value = [
                'test-lab:IGVFDS123ABC_bucket_sub123_wf456_subwf789_rna_kb_h5ad_uniform-pipeline']

            result = matrix_payload.aliases
            assert len(result) == 1
            assert 'test-lab:' in result[0]
            assert 'IGVFDS123ABC' in result[0]

    def test_submitted_file_name_property(self, matrix_payload):
        """Test submitted_file_name property."""
        result = matrix_payload.submitted_file_name
        assert result == 'gs://bucket/file.h5ad'

    @patch('sc_pipe_management.accession.igvf_payloads.api_tools.calculate_gsfile_hex_hash')
    def test_md5sum_property(self, mock_hash, matrix_payload):
        """Test md5sum property."""
        mock_hash.return_value = 'abc123def456'

        result = matrix_payload.md5sum
        assert result == 'abc123def456'
        mock_hash.assert_called_once_with(file_path='gs://bucket/file.h5ad')

    def test_get_payload(self, matrix_payload):
        """Test get_payload method."""
        with patch('sc_pipe_management.accession.igvf_payloads._get_file_aliases') as mock_aliases, \
                patch('sc_pipe_management.accession.igvf_payloads.api_tools.calculate_gsfile_hex_hash') as mock_hash:

            mock_aliases.return_value = ['test-alias']
            mock_hash.return_value = 'abc123'

            payload = matrix_payload.get_payload()

            # Check required fields
            assert payload['award'] == '/awards/test-award/'
            assert payload['lab'] == '/labs/test-lab/'
            assert payload['aliases'] == ['test-alias']
            assert payload['md5sum'] == 'abc123'
            assert payload['file_set'] == 'IGVFDS123ABC'
            assert payload['submitted_file_name'] == 'gs://bucket/file.h5ad'
            assert payload['_profile'] == 'matrix_file'
            assert payload['principal_dimension'] == 'cell'
            assert payload['secondary_dimensions'] == ['gene']
            assert payload['filtered'] is False


class TestAlignmentFilePayload:
    """Test AlignmentFilePayload class."""

    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance for alignment files."""
        mock_metadata = Mock()
        mock_metadata.terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC',
            'atac_bam': 'gs://bucket/alignment.bam'
        })

        # Mock UUIDs
        mock_uuids = Mock()
        mock_uuids.aliases.return_value = 'bucket_sub123_wf456_subwf789'
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        # Mock input file accessions
        mock_input_files = Mock()
        mock_input_files.get_derived_from.return_value = ['IGVFFF001AAA']
        mock_input_files.reference_files = ['IGVFFF002BBB']
        mock_input_files.sequence_files = ['IGVFFF003CCC']
        mock_metadata._get_input_file_accs_from_table.return_value = mock_input_files

        return mock_metadata

    @pytest.fixture
    def mock_igvf_api(self):
        """Mock IGVF API."""
        mock_api = Mock()
        # Mock sequence file response
        mock_seq_file = Mock()
        mock_seq_file.controlled_access = False
        mock_api.get_by_id.return_value.actual_instance = mock_seq_file
        return mock_api

    @pytest.fixture
    @patch('sc_pipe_management.accession.igvf_payloads.const.OUTPUT_SUBMITTER_INFO',
           {'lab': '/labs/test-lab/', 'award': '/awards/test-award/'})
    @patch('sc_pipe_management.accession.igvf_payloads.const.ALIGNMENT_FILETYPES')
    def alignment_payload(self, mock_alignment_types, mock_terra_metadata, mock_igvf_api):
        """Create AlignmentFilePayload instance for testing."""
        # Mock file metadata
        mock_file_metadata = Mock()
        mock_file_metadata.assay_type = 'atac'
        mock_file_metadata.analysis_step_version = '/analysis-step-versions/test/'
        mock_file_metadata.file_format = 'bam'
        mock_file_metadata.content_type = 'alignments'
        mock_file_metadata.description = 'Test alignment file'

        mock_alignment_types.__getitem__.return_value = mock_file_metadata

        return igvf_payloads.AlignmentFilePayload(mock_terra_metadata, mock_igvf_api)

    def test_init(self, alignment_payload, mock_terra_metadata, mock_igvf_api):
        """Test AlignmentFilePayload initialization."""
        assert alignment_payload.terra_output_name == 'atac_bam'
        assert alignment_payload.lab == '/labs/test-lab/'
        assert alignment_payload.award == '/awards/test-award/'
        assert alignment_payload.igvf_api == mock_igvf_api

    def test_get_access_status_open_access(self, alignment_payload):
        """Test _get_access_status returns False for open access files."""
        result = alignment_payload._get_access_status()
        assert result is False

    def test_get_access_status_controlled_access(self, alignment_payload):
        """Test _get_access_status returns True when any file is controlled access."""
        # Mock one controlled access file
        mock_seq_file = Mock()
        mock_seq_file.controlled_access = True
        alignment_payload.igvf_api.get_by_id.return_value.actual_instance = mock_seq_file

        result = alignment_payload._get_access_status()
        assert result is True

    def test_get_payload(self, alignment_payload):
        """Test get_payload method."""
        with patch('sc_pipe_management.accession.igvf_payloads._get_file_aliases') as mock_aliases, \
                patch('sc_pipe_management.accession.igvf_payloads.api_tools.calculate_gsfile_hex_hash') as mock_hash:

            mock_aliases.return_value = ['test-alias']
            mock_hash.return_value = 'abc123'

            payload = alignment_payload.get_payload()

            # Check required fields
            assert payload['award'] == '/awards/test-award/'
            assert payload['lab'] == '/labs/test-lab/'
            assert payload['aliases'] == ['test-alias']
            assert payload['md5sum'] == 'abc123'
            assert payload['file_set'] == 'IGVFDS123ABC'
            assert payload['_profile'] == 'alignment_file'
            assert payload['filtered'] is False
            assert payload['redacted'] is False
            assert 'controlled_access' in payload


class TestFragmentFilePayload:
    """Test FragmentFilePayload class."""

    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance for fragment files."""
        mock_metadata = Mock()
        mock_metadata.terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC',
            'atac_fragments': 'gs://bucket/fragments.tsv.gz'
        })
        mock_metadata.taxa = 'Homo sapiens'

        # Mock UUIDs
        mock_uuids = Mock()
        mock_uuids.aliases.return_value = 'bucket_sub123_wf456_subwf789'
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        # Mock input file accessions
        mock_input_files = Mock()
        mock_input_files.get_derived_from.return_value = ['IGVFFF001AAA']
        mock_input_files.description = 'Test fragment file'
        mock_metadata._get_input_file_accs_from_table.return_value = mock_input_files

        return mock_metadata

    @pytest.fixture
    def mock_igvf_api(self):
        """Mock IGVF API."""
        return Mock()

    @pytest.fixture
    @patch('sc_pipe_management.accession.igvf_payloads.const.OUTPUT_SUBMITTER_INFO',
           {'lab': '/labs/test-lab/', 'award': '/awards/test-award/'})
    @patch('sc_pipe_management.accession.igvf_payloads.const.TABULAR_FILETYPES')
    @patch('sc_pipe_management.accession.igvf_payloads.const.GENOME_ASSEMBLY_INFO')
    def fragment_payload(self, mock_genome_info, mock_tabular_types, mock_terra_metadata, mock_igvf_api):
        """Create FragmentFilePayload instance for testing."""
        # Mock file metadata
        mock_file_metadata = Mock()
        mock_file_metadata.assay_type = 'atac'
        mock_file_metadata.analysis_step_version = '/analysis-step-versions/test/'
        mock_file_metadata.content_type = 'fragments'
        mock_file_metadata.file_format = 'tsv'
        mock_file_metadata.file_format_specifications = [
            '/file-format-specifications/test/']

        mock_tabular_types.__getitem__.return_value = mock_file_metadata
        mock_genome_info.get.return_value = 'GRCh38'

        return igvf_payloads.FragmentFilePayload(mock_terra_metadata, mock_igvf_api)

    def test_init(self, fragment_payload, mock_igvf_api):
        """Test FragmentFilePayload initialization."""
        assert fragment_payload.terra_output_name == 'atac_fragments'
        assert fragment_payload.lab == '/labs/test-lab/'
        assert fragment_payload.award == '/awards/test-award/'
        assert fragment_payload.igvf_api == mock_igvf_api

    def test_get_payload(self, fragment_payload):
        """Test get_payload method."""
        with patch('sc_pipe_management.accession.igvf_payloads._get_file_aliases') as mock_aliases, \
                patch('sc_pipe_management.accession.igvf_payloads.api_tools.calculate_gsfile_hex_hash') as mock_hash:

            mock_aliases.return_value = ['test-alias']
            mock_hash.return_value = 'abc123'

            payload = fragment_payload.get_payload()

            # Check required fields
            assert payload['award'] == '/awards/test-award/'
            assert payload['lab'] == '/labs/test-lab/'
            assert payload['aliases'] == ['test-alias']
            assert payload['md5sum'] == 'abc123'
            assert payload['file_set'] == 'IGVFDS123ABC'
            assert payload['_profile'] == 'tabular_file'
            assert payload['file_format_type'] == 'bed3+'
            assert payload['controlled_access'] is False
            assert payload['filtered'] is False
            assert payload['assembly'] == 'GRCh38'


class TestIndexFilePayload:
    """Test IndexFilePayload class."""

    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance for index files."""
        mock_metadata = Mock()
        mock_metadata.terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC',
            'atac_bam_index': 'gs://bucket/alignment.bam.bai'
        })

        # Mock UUIDs
        mock_uuids = Mock()
        mock_uuids.aliases.return_value = 'bucket_sub123_wf456_subwf789'
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        return mock_metadata

    @pytest.fixture
    def mock_igvf_api(self):
        """Mock IGVF API."""
        return Mock()

    @pytest.fixture
    @patch('sc_pipe_management.accession.igvf_payloads.const.OUTPUT_SUBMITTER_INFO',
           {'lab': '/labs/test-lab/', 'award': '/awards/test-award/'})
    @patch('sc_pipe_management.accession.igvf_payloads.const.INDEX_FILETYPES')
    def index_payload(self, mock_index_types, mock_terra_metadata, mock_igvf_api):
        """Create IndexFilePayload instance for testing."""
        # Mock file metadata
        mock_file_metadata = Mock()
        mock_file_metadata.analysis_step_version = '/analysis-step-versions/test/'
        mock_file_metadata.controlled_access = False
        mock_file_metadata.file_format = 'bai'
        mock_file_metadata.description = 'Test index file'

        mock_index_types.__getitem__.return_value = mock_file_metadata

        derived_from = ['IGVFFF001AAA']
        return igvf_payloads.IndexFilePayload(
            mock_terra_metadata, derived_from, 'atac_bam_index', mock_igvf_api
        )

    def test_init(self, index_payload):
        """Test IndexFilePayload initialization."""
        assert index_payload.terra_output_name == 'atac_bam_index'
        assert index_payload.lab == '/labs/test-lab/'
        assert index_payload.award == '/awards/test-award/'
        assert index_payload.derived_from == ['IGVFFF001AAA']

    def test_get_payload(self, index_payload):
        """Test get_payload method."""
        with patch('sc_pipe_management.accession.igvf_payloads._get_file_aliases') as mock_aliases, \
                patch('sc_pipe_management.accession.igvf_payloads.api_tools.calculate_gsfile_hex_hash') as mock_hash:

            mock_aliases.return_value = ['test-alias']
            mock_hash.return_value = 'abc123'

            payload = index_payload.get_payload()

            # Check required fields
            assert payload['award'] == '/awards/test-award/'
            assert payload['lab'] == '/labs/test-lab/'
            assert payload['aliases'] == ['test-alias']
            assert payload['md5sum'] == 'abc123'
            assert payload['file_set'] == 'IGVFDS123ABC'
            assert payload['_profile'] == 'index_file'
            assert payload['content_type'] == 'index'
            assert payload['controlled_access'] is False
            assert payload['derived_from'] == ['IGVFFF001AAA']


class TestQCMetricsPayload:
    """Test QCMetricsPayload class."""

    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance for QC metrics."""
        mock_metadata = Mock()
        mock_metadata.terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC',
            'rna_kb_library_qc_metrics_json': 'gs://bucket/qc_metrics.json',
            'rna_kb_parameters_json': 'gs://bucket/parameters.json'
        })

        # Mock UUIDs
        mock_uuids = Mock()
        mock_uuids.aliases.return_value = 'bucket_sub123_wf456_subwf789'
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        return mock_metadata

    @pytest.fixture
    def mock_qc_info_map(self):
        """Mock QC info map."""
        mock_info = Mock()
        mock_info.analysis_step_version = '/analysis-step-versions/test/'
        mock_info.description = 'Test QC metric'
        mock_info.object_type = 'single_cell_rna_seq_quality_metric'
        mock_info.__getitem__ = Mock(side_effect=lambda key: {
            'metadata': ['rna_kb_library_qc_metrics_json'],
            'attachment': {'rnaseq_kb_info': 'rna_kb_parameters_json'},
            'metadata_map': {'numRecords': 'n_records', 'numReads': 'n_reads'}
        }[key])
        return mock_info

    @pytest.fixture
    def mock_igvf_api(self):
        """Mock IGVF API."""
        return Mock()

    @pytest.fixture
    @patch('sc_pipe_management.accession.igvf_payloads.const.OUTPUT_SUBMITTER_INFO',
           {'lab': '/labs/test-lab/', 'award': '/awards/test-award/'})
    def qc_payload(self, mock_terra_metadata, mock_qc_info_map, mock_igvf_api):
        """Create QCMetricsPayload instance for testing."""
        with tempfile.TemporaryDirectory() as temp_dir:
            return igvf_payloads.QCMetricsPayload(
                terra_metadata=mock_terra_metadata,
                qc_info_map=mock_qc_info_map,
                qc_prefix='gene_count',
                qc_of=['IGVFFF001AAA'],
                igvf_api=mock_igvf_api,
                root_output_dir=temp_dir
            )

    def test_init(self, qc_payload, mock_terra_metadata, mock_igvf_api):
        """Test QCMetricsPayload initialization."""
        assert qc_payload.terra_output_name == 'gene_count_metrics'
        assert qc_payload.lab == '/labs/test-lab/'
        assert qc_payload.award == '/awards/test-award/'
        assert qc_payload.qc_of == ['IGVFFF001AAA']
        assert qc_payload.igvf_api == mock_igvf_api

    def test_submitted_file_name_property(self, qc_payload):
        """Test submitted_file_name property returns None."""
        assert qc_payload.submitted_file_name is None

    def test_md5sum_property(self, qc_payload):
        """Test md5sum property returns None."""
        assert qc_payload.md5sum is None

    def test_aliases_property(self, qc_payload):
        """Test aliases property."""
        result = qc_payload.aliases
        assert len(result) == 1
        assert 'test-lab:' in result[0]
        assert 'IGVFDS123ABC' in result[0]
        assert 'gene_count_metrics' in result[0]

    @patch('sc_pipe_management.accession.igvf_payloads._download_qc_file_from_gcp')
    def test_get_qc_files(self, mock_download, qc_payload):
        """Test _get_qc_files method."""
        # Mock successful downloads
        mock_download.side_effect = [
            '/tmp/qc_metrics.json', '/tmp/parameters.json']

        result = qc_payload._get_qc_files()

        assert isinstance(result, igvf_payloads.QCFileDownloadInfo)
        assert len(result.paths_of_metadata_files) == 1
        assert len(result.paths_of_attachment_files) == 1
        assert result.paths_of_metadata_files[0] == '/tmp/qc_metrics.json'

    def test_read_json_file(self):
        """Test _read_json_file static method."""
        test_data = {"numRecords": 1000, "numReads": 50000}

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_data, f)
            temp_file = f.name

        try:
            result = igvf_payloads.QCMetricsPayload._read_json_file([
                                                                    temp_file])
            assert result == test_data
        finally:
            os.unlink(temp_file)


class TestPipelineParamsInfo:
    """Test PipelineParamsInfo class."""

    @pytest.fixture
    def sample_terra_datable(self):
        """Sample Terra data table."""
        return pd.DataFrame({
            'analysis_set_acc': ['IGVFDS123ABC', 'IGVFDS456DEF'],
            'rna_kb_h5ad': [
                'gs://bucket/sub123/wf456/subwf789/file1.h5ad',
                'gs://bucket/sub456/wf789/subwf012/file2.h5ad'
            ]
        })

    @pytest.fixture
    def pipeline_params(self, sample_terra_datable):
        """Create PipelineParamsInfo instance."""
        with tempfile.TemporaryDirectory() as temp_dir:
            return igvf_payloads.PipelineParamsInfo(
                terra_datable=sample_terra_datable,
                output_root_dir=temp_dir
            )

    def test_init(self, pipeline_params, sample_terra_datable):
        """Test PipelineParamsInfo initialization."""
        assert pipeline_params.terra_namespace == 'DACC_ANVIL'
        assert pipeline_params.terra_workspace == 'IGVF Single-Cell Data Processing'
        assert pipeline_params.terra_datable.equals(sample_terra_datable)

    def test_init_with_custom_params(self, sample_terra_datable):
        """Test initialization with custom parameters."""
        with tempfile.TemporaryDirectory() as temp_dir:
            params = igvf_payloads.PipelineParamsInfo(
                terra_datable=sample_terra_datable,
                terra_namespace='CUSTOM_NAMESPACE',
                terra_workspace='Custom Workspace',
                output_root_dir=temp_dir
            )

            assert params.terra_namespace == 'CUSTOM_NAMESPACE'
            assert params.terra_workspace == 'Custom Workspace'
            assert params.output_root_dir == temp_dir

    @patch('sc_pipe_management.accession.igvf_payloads.fapi.get_workflow_metadata')
    @patch('sc_pipe_management.accession.igvf_payloads._dump_json')
    def test_get_single_input_params_success(self, mock_dump_json, mock_get_metadata, pipeline_params):
        """Test successful input params retrieval."""
        # Mock API response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            'inputs': {'param1': 'value1', 'param2': 'value2'}
        }
        mock_get_metadata.return_value = mock_response

        # Mock file dump
        mock_dump_json.return_value = '/path/to/config.json'

        # Mock terra metadata
        mock_terra_metadata = Mock()
        mock_terra_metadata.analysis_set_acc = 'IGVFDS123ABC'
        mock_uuids = Mock()
        mock_uuids.submission_id = 'sub123'
        mock_uuids.workflow_id = 'wf456'
        mock_terra_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        result = pipeline_params._get_single_input_params(mock_terra_metadata)

        assert result == '/path/to/config.json'
        mock_get_metadata.assert_called_once_with(
            namespace='DACC_ANVIL',
            workspace='IGVF Single-Cell Data Processing',
            submission_id='sub123',
            workflow_id='wf456'
        )

    @patch('sc_pipe_management.accession.igvf_payloads.fapi.get_workflow_metadata')
    def test_get_single_input_params_api_error(self, mock_get_metadata, pipeline_params):
        """Test API error handling."""
        # Mock failed API response
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get_metadata.return_value = mock_response

        mock_terra_metadata = Mock()
        mock_terra_metadata.analysis_set_acc = 'IGVFDS123ABC'

        with pytest.raises(Exception):  # FireCloudServerError
            pipeline_params._get_single_input_params(mock_terra_metadata)


class TestDocumentPayload:
    """Test DocumentPayload class."""

    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance."""
        mock_metadata = Mock()
        mock_metadata.anaset_accession = 'IGVFDS123ABC'

        # Mock UUIDs
        mock_uuids = Mock()
        mock_uuids.input_param_aliases.return_value = 'bucket_sub123_wf456_subwf789'
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        return mock_metadata

    @pytest.fixture
    def mock_pipeline_params_info(self):
        """Mock pipeline params info."""
        return {'IGVFDS123ABC': '/path/to/config.json'}

    @pytest.fixture
    def mock_igvf_api(self):
        """Mock IGVF API."""
        return Mock()

    @pytest.fixture
    @patch('sc_pipe_management.accession.igvf_payloads.const.OUTPUT_SUBMITTER_INFO',
           {'lab': '/labs/test-lab/', 'award': '/awards/test-award/'})
    def document_payload(self, mock_terra_metadata, mock_pipeline_params_info, mock_igvf_api):
        """Create DocumentPayload instance."""
        return igvf_payloads.DocumentPayload(
            terra_metadata=mock_terra_metadata,
            pipeline_params_info=mock_pipeline_params_info,
            igvf_api=mock_igvf_api
        )

    def test_init(self, document_payload, mock_terra_metadata, mock_igvf_api):
        """Test DocumentPayload initialization."""
        assert document_payload.terra_output_name == 'pipeline_parameters'
        assert document_payload.lab == '/labs/test-lab/'
        assert document_payload.award == '/awards/test-award/'
        assert document_payload.input_params_file_path == '/path/to/config.json'
        assert document_payload.igvf_api == mock_igvf_api

    def test_mk_doc_aliases(self, document_payload):
        """Test _mk_doc_aliases method."""
        result = document_payload._mk_doc_aliases()
        assert len(result) == 1
        assert 'igvf-dacc-processing-pipeline:' in result[0]
        assert 'bucket_sub123_wf456_subwf789' in result[0]
        assert '_pipeline_config' in result[0]

    def test_get_payload(self, document_payload):
        """Test get_payload method."""
        payload = document_payload.get_payload()

        # Check required fields
        assert payload['award'] == '/awards/test-award/'
        assert payload['lab'] == '/labs/test-lab/'
        assert payload['content_type'] == 'application/json'
        assert payload['document_type'] == 'pipeline parameters'
        assert payload['file_format'] == 'json'
        assert payload['_profile'] == 'document'
        assert payload['attachment']['path'] == '/path/to/config.json'
        assert 'description' in payload
        assert len(payload['aliases']) == 1


class TestAnalysisSetPatchingPayload:
    """Test AnalysisSetPatchingPayload class."""

    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance."""
        mock_metadata = Mock()
        mock_metadata.analysis_set_acc = 'IGVFDS123ABC'
        return mock_metadata

    @pytest.fixture
    def mock_igvf_utils_api(self):
        """Mock IGVF utils API."""
        mock_api = Mock()
        mock_api.IGVFID_KEY = '@id'
        return mock_api

    @pytest.fixture
    def patching_payload(self, mock_terra_metadata, mock_igvf_utils_api):
        """Create AnalysisSetPatchingPayload instance."""
        return igvf_payloads.AnalysisSetPatchingPayload(
            terra_metadata=mock_terra_metadata,
            input_params_doc_uuid='doc-uuid-123',
            igvf_utils_api=mock_igvf_utils_api
        )

    def test_init(self, patching_payload, mock_terra_metadata, mock_igvf_utils_api):
        """Test AnalysisSetPatchingPayload initialization."""
        assert patching_payload.analysis_set_acc == 'IGVFDS123ABC'
        assert patching_payload.input_params_doc_uuid == 'doc-uuid-123'
        assert patching_payload.igvf_utils_api == mock_igvf_utils_api

    def test_get_existing_analysis_set_docs_no_docs(self, patching_payload):
        """Test _get_existing_analysis_set_docs when no documents exist."""
        # Mock analysis set with no documents
        patching_payload.igvf_utils_api.get.return_value = {'documents': None}

        result = patching_payload._get_existing_analysis_set_docs()
        assert result == []

    def test_get_existing_analysis_set_docs_with_docs(self, patching_payload):
        """Test _get_existing_analysis_set_docs with existing documents."""
        # Mock analysis set with documents
        patching_payload.igvf_utils_api.get.side_effect = [
            {'documents': ['/documents/doc1/', '/documents/doc2/']},
            {'uuid': 'uuid1'},
            {'uuid': 'uuid2'}
        ]

        result = patching_payload._get_existing_analysis_set_docs()
        assert len(result) == 2
        assert 'uuid1' in result
        assert 'uuid2' in result

    def test_get_patch_payload_doc_already_exists(self, patching_payload):
        """Test _get_patch_payload when document already exists."""
        with patch.object(patching_payload, '_get_existing_analysis_set_docs',
                          return_value=['doc-uuid-123', 'other-uuid']):
            result = patching_payload._get_patch_payload()
            assert result is None

    def test_get_patch_payload_new_doc(self, patching_payload):
        """Test _get_patch_payload when document is new."""
        with patch.object(patching_payload, '_get_existing_analysis_set_docs',
                          return_value=['other-uuid']):
            result = patching_payload._get_patch_payload()

            assert result is not None
            assert result['@id'] == '/analysis-sets/IGVFDS123ABC/'
            assert result['uniform_pipeline_status'] == 'completed'
            assert result['_profile'] == 'analysis_set'
            assert 'documents' in result


if __name__ == '__main__':
    pytest.main([__file__, "-v"])
