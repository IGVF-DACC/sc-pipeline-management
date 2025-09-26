import pytest
from unittest.mock import Mock, MagicMock, patch
import json
import tempfile
import os
import sys
import hashlib


# Add the absolute path to the 'src' directory to sys.path
src_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.wrangler_utils.check_accession_results as qa_run
import sc_pipe_management.igvf_and_terra_api_tools as api_tools


def compute_md5sum(file_path):
    """Compute MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


@pytest.fixture
def mock_igvf_client():
    return Mock()


@pytest.fixture
def mock_file_obj():
    mock_file = Mock()
    mock_file.type = ['MatrixFile']
    mock_file.accession = 'IGVFDS123456'
    mock_file.analysis_step_version = '/analysis-step-versions/test-version/'
    mock_file.upload_status = 'validated'
    mock_file.controlled_access = True
    mock_file.reference_files = ['/reference-files/test-ref/']
    mock_file.submitted_file_name = 'IGVFDS123456_test_file.h5ad'
    mock_file.input_file_for = []
    mock_file.quality_metrics = ['/quality-metrics/test-qm/']
    return mock_file


@pytest.fixture
def mock_analysis_set_obj():
    mock_analysis_set = Mock()
    mock_analysis_set.files = ['/matrix-files/test1/', '/matrix-files/test2/']
    mock_analysis_set.pipeline_parameters = ['/documents/test-doc/']
    return mock_analysis_set


class TestUtilityFunctions:

    def test_get_raw_data_access_level_mixed_controlled_access(self, mock_igvf_client, mock_analysis_set_obj):
        """Test _check_if_controlled_access returns True when at least one sequence file has controlled access."""
        # Create sequence files with mixed controlled_access status
        mock_seq_file_1 = Mock()
        mock_seq_file_1.controlled_access = False  # Not controlled

        mock_seq_file_2 = Mock()
        mock_seq_file_2.controlled_access = True   # Controlled

        mock_seq_file_3 = Mock()
        mock_seq_file_3.controlled_access = False  # Not controlled

        # Create file set with multiple sequence files
        mock_file_set = Mock()
        mock_file_set.files = [
            '/sequence-files/test-seq-1/',
            '/sequence-files/test-seq-2/',
            '/sequence-files/test-seq-3/'
        ]

        # Mock responses for different file IDs
        mock_responses = {
            '/measurement-sets/test-set/': Mock(actual_instance=mock_file_set),
            '/sequence-files/test-seq-1/': Mock(actual_instance=mock_seq_file_1),
            '/sequence-files/test-seq-2/': Mock(actual_instance=mock_seq_file_2),
            '/sequence-files/test-seq-3/': Mock(actual_instance=mock_seq_file_3)
        }
        mock_igvf_client.get_by_id.side_effect = lambda x: mock_responses[x]

        # Mock analysis set object with input file sets
        mock_analysis_set_obj.input_file_sets = ['/measurement-sets/test-set/']

        result = qa_run._check_if_controlled_access(
            mock_analysis_set_obj, mock_igvf_client)

        # Should return True because at least one sequence file has controlled_access=True
        assert result is True

    def test_get_raw_data_access_level_all_false(self, mock_igvf_client, mock_analysis_set_obj):
        """Test _check_if_controlled_access returns False when no sequence files have controlled access."""
        # Create sequence files with all controlled_access = False
        mock_seq_file_1 = Mock()
        mock_seq_file_1.controlled_access = False

        mock_seq_file_2 = Mock()
        mock_seq_file_2.controlled_access = False

        mock_seq_file_3 = Mock()
        mock_seq_file_3.controlled_access = False

        # Create file set with multiple sequence files
        mock_file_set = Mock()
        mock_file_set.files = [
            '/sequence-files/test-seq-1/',
            '/sequence-files/test-seq-2/',
            '/sequence-files/test-seq-3/'
        ]

        # Mock responses for different file IDs
        mock_responses = {
            '/measurement-sets/test-set/': Mock(actual_instance=mock_file_set),
            '/sequence-files/test-seq-1/': Mock(actual_instance=mock_seq_file_1),
            '/sequence-files/test-seq-2/': Mock(actual_instance=mock_seq_file_2),
            '/sequence-files/test-seq-3/': Mock(actual_instance=mock_seq_file_3)
        }
        mock_igvf_client.get_by_id.side_effect = lambda x: mock_responses[x]

        # Mock analysis set object with input file sets
        mock_analysis_set_obj.input_file_sets = ['/measurement-sets/test-set/']

        result = qa_run._check_if_controlled_access(
            mock_analysis_set_obj, mock_igvf_client)

        # Should return False because no sequence files have controlled_access=True
        assert result is False


class TestBaseFileChecker:

    def test_check_basic_file_properties_success(self, mock_igvf_client, mock_file_obj):
        """Test _check_basic_file_properties with valid file."""
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)
        # Should return empty list (no errors)
        errors = checker._check_basic_file_properties(
            mock_file_obj, '/analysis-step-versions/test-version/')
        assert errors == []

    def test_check_basic_file_properties_no_asv(self, mock_igvf_client, mock_file_obj):
        """Test _check_basic_file_properties with missing ASV."""
        mock_file_obj.analysis_step_version = None
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_basic_file_properties(
            mock_file_obj, '/analysis-step-versions/test-version/')

        # Should return list with error message
        assert len(errors) > 0
        assert any("has no analysis step version" in error for error in errors)

    def test_check_basic_file_properties_wrong_asv(self, mock_igvf_client, mock_file_obj):
        """Test _check_basic_file_properties with wrong ASV."""
        mock_file_obj.analysis_step_version = '/analysis-step-versions/wrong-version/'
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_basic_file_properties(
            mock_file_obj, '/analysis-step-versions/test-version/')

        # Should return list with error message
        assert len(errors) > 0
        assert any(
            "has unexpected analysis step version" in error for error in errors)

    def test_check_basic_file_properties_invalid_upload_status(self, mock_igvf_client, mock_file_obj):
        """Test _check_basic_file_properties with invalid upload status."""
        mock_file_obj.upload_status = 'pending'
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_basic_file_properties(
            mock_file_obj, '/analysis-step-versions/test-version/')

        # Should return list with error message
        assert len(errors) > 0
        assert any("has invalid upload status" in error for error in errors)

    def test_check_basic_file_properties_controlled_access_mismatch(self, mock_igvf_client, mock_file_obj):
        """Test _check_basic_file_properties with controlled access mismatch."""
        mock_file_obj.controlled_access = False
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_basic_file_properties(
            mock_file_obj, '/analysis-step-versions/test-version/', True)

        # Should return list with error message
        assert len(errors) > 0
        assert any(
            "should be controlled access True" in error for error in errors)

    def test_check_basic_file_properties_multiple_errors(self, mock_igvf_client, mock_file_obj):
        """Test _check_basic_file_properties with multiple errors."""
        mock_file_obj.analysis_step_version = None
        mock_file_obj.upload_status = 'pending'
        mock_file_obj.controlled_access = False

        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_basic_file_properties(
            mock_file_obj, '/analysis-step-versions/test-version/', True)

        # Should return list with multiple error messages
        assert len(errors) == 3
        assert any("has no analysis step version" in error for error in errors)
        assert any("has invalid upload status" in error for error in errors)
        assert any(
            "should be controlled access True" in error for error in errors)

    def test_check_reference_files_success(self, mock_igvf_client, mock_file_obj):
        """Test _check_reference_files with valid reference files."""
        reference_list = ['/reference-files/test-ref/',
                          '/reference-files/another-ref/']
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_reference_files(mock_file_obj, reference_list)
        assert errors == []

    def test_check_reference_files_none(self, mock_igvf_client, mock_file_obj):
        """Test _check_reference_files with None reference files."""
        mock_file_obj.reference_files = None
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_reference_files(
            mock_file_obj, ['/reference-files/test-ref/'])

        assert len(errors) > 0
        assert any("has no reference files" in error for error in errors)

    def test_check_file_name_success(self, mock_igvf_client, mock_file_obj):
        """Test _check_file_name with valid file name."""
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_file_name(mock_file_obj, 'IGVFDS123456')
        assert errors == []

    def test_check_file_name_invalid(self, mock_igvf_client, mock_file_obj):
        """Test _check_file_name with invalid file name."""
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        errors = checker._check_file_name(mock_file_obj, 'WRONG_NAME')

        assert len(errors) > 0
        assert any("has unexpected file name" in error for error in errors)


class TestDocumentChecker:

    def test_check_pipeline_parameters_document_success(self, mock_igvf_client):
        """Test check_pipeline_parameters_document with valid document."""
        mock_doc = Mock()
        mock_doc.document_type = 'pipeline parameters'

        # Create a proper Attachment mock object
        mock_attachment = Mock()
        mock_attachment.download = 'IGVFDS1333CJLG_run_config.json'
        mock_attachment.href = '@@download/attachment/IGVFDS1333CJLG_run_config.json'
        mock_attachment.type = 'application/json'
        mock_attachment.md5sum = 'ccf30880b1b7f10953490b6101e02e32'
        mock_attachment.size = None
        mock_attachment.width = None
        mock_attachment.height = None

        mock_doc.attachment = mock_attachment

        mock_response = Mock()
        mock_response.actual_instance = mock_doc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.DocumentChecker('IGVFDS123456', mock_igvf_client)
        errors = checker.check_pipeline_parameters_document(
            ['/documents/test-doc/'])
        assert errors == []

    def test_check_pipeline_parameters_document_wrong_type(self, mock_igvf_client):
        """Test check_pipeline_parameters_document with wrong document type."""
        mock_doc = Mock()
        mock_doc.document_type = 'other document'

        # Create a proper Attachment mock object
        mock_attachment = Mock()
        mock_attachment.download = 'IGVFDS1333CJLG_run_config.json'
        mock_attachment.href = '@@download/attachment/IGVFDS1333CJLG_run_config.json'
        mock_attachment.type = 'application/json'
        mock_attachment.md5sum = 'ccf30880b1b7f10953490b6101e02e32'
        mock_attachment.size = None
        mock_attachment.width = None
        mock_attachment.height = None

        mock_doc.attachment = mock_attachment

        mock_response = Mock()
        mock_response.actual_instance = mock_doc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.DocumentChecker('IGVFDS123456', mock_igvf_client)

        errors = checker.check_pipeline_parameters_document(
            ['/documents/test-doc/'])

        assert len(errors) > 0
        assert any(
            "Unexpected pipeline parameters document type" in error for error in errors)

    def test_check_pipeline_parameters_document_no_docs(self, mock_igvf_client):
        """Test check_pipeline_parameters_document with no documents."""
        checker = qa_run.DocumentChecker('IGVFDS123456', mock_igvf_client)

        errors = checker.check_pipeline_parameters_document([])

        assert len(errors) > 0
        assert any(
            "No pipeline parameters document found" in error for error in errors)


class TestQCMetricsChecker:

    def test_check_rna_qc_metrics_success(self, mock_igvf_client, mock_file_obj):
        """Test check_rna_qc_metrics with valid QC metrics."""
        mock_qc = Mock()
        mock_qc.type = ['SingleCellRnaSeqQualityMetric']
        mock_qc.analysis_step_version = '/analysis-step-versions/test-version/'
        mock_qc.rnaseq_kb_info = 'some-attachment'
        mock_qc.n_reads = 1000000

        mock_response = Mock()
        mock_response.actual_instance = mock_qc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.QCMetricsChecker('IGVFDS123456', mock_igvf_client)
        errors = checker.check_rna_qc_metrics(mock_file_obj)
        assert errors == []

    def test_check_atac_qc_metrics_alignment_success(self, mock_igvf_client, mock_file_obj):
        """Test check_atac_qc_metrics with valid alignment file QC metrics."""
        mock_file_obj.type = ['AlignmentFile']

        mock_qc = Mock()
        mock_qc.type = ['SingleCellAtacSeqQualityMetric']
        mock_qc.analysis_step_version = '/analysis-step-versions/test-version/'
        mock_qc.atac_bam_summary_stats = 'some-attachment'

        mock_response = Mock()
        mock_response.actual_instance = mock_qc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.QCMetricsChecker('IGVFDS123456', mock_igvf_client)
        errors = checker.check_atac_qc_metrics(mock_file_obj)
        assert errors == []


class TestIndexFileChecker:

    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_active_file_objs')
    def test_check_index_file_presence_and_format_alignment(self, mock_get_active, mock_igvf_client):
        """Test check_index_file_presence_and_format for AlignmentFile."""
        mock_file_obj = Mock()
        mock_file_obj.type = ['AlignmentFile']
        mock_file_obj.input_file_for = ['/index-files/test.bai']
        mock_file_obj.controlled_access = True

        mock_index = Mock()
        mock_index.file_format = 'bai'
        mock_index.controlled_access = True
        mock_index.submitted_file_name = 'IGVFDS123456.bam.bai'
        mock_index.accession = 'IGVFFI123456'

        mock_get_active.return_value = [mock_index]

        checker = qa_run.IndexFileChecker('IGVFDS123456', mock_igvf_client)
        errors = checker.check_index_file_presence_and_format(mock_file_obj)
        assert errors == []


class TestRNASeqFileChecker:

    def test_check_rna_file_count_success(self, mock_igvf_client):
        """Test check_rna_file_count with exactly 2 files."""
        checker = qa_run.RNASeqFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock(), Mock()]
        errors = checker.check_rna_file_count(mock_files)
        assert errors == []

    def test_check_rna_file_count_wrong_count(self, mock_igvf_client):
        """Test check_rna_file_count with wrong number of files."""
        checker = qa_run.RNASeqFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock()]

        errors = checker.check_rna_file_count(mock_files)

        assert len(errors) > 0
        assert any(
            "does not have exactly 2 non-deprecated RNA-seq output data" in error for error in errors)

    def test_check_rna_file_format_and_content_h5ad(self, mock_igvf_client):
        """Test check_rna_file_format_and_content with h5ad format."""
        mock_file = Mock()
        mock_file.file_format = 'h5ad'
        mock_file.content_type = 'sparse gene count matrix'
        mock_file.submitted_file_name = 'IGVFDS123456.h5ad'
        mock_file.type = ['MatrixFile']
        mock_file.accession = 'IGVFFI123456'

        checker = qa_run.RNASeqFileChecker('IGVFDS123456', mock_igvf_client)
        errors = checker.check_rna_file_format_and_content(mock_file)
        assert errors == []


class TestATACFileChecker:

    def test_check_atac_file_count_success(self, mock_igvf_client):
        """Test check_atac_file_count with exactly 1 file."""
        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock()]
        errors = checker.check_atac_file_count(mock_files, "alignment")
        assert errors == []

    def test_check_atac_file_count_wrong_count(self, mock_igvf_client):
        """Test check_atac_file_count with wrong number of files."""
        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock(), Mock()]

        errors = checker.check_atac_file_count(mock_files, "alignment")

        assert len(errors) > 0
        assert any("has more than 1 alignment file" in error for error in errors)

    def test_check_alignment_file_specifics(self, mock_igvf_client):
        """Test check_alignment_file_specifics with valid alignment file."""
        mock_file = Mock()
        mock_file.file_format = 'bam'
        mock_file.content_type = 'alignments'
        mock_file.submitted_file_name = 'IGVFDS123456.bam'
        mock_file.type = ['AlignmentFile']
        mock_file.accession = 'IGVFFI123456'

        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)
        errors = checker.check_alignment_file_specifics(mock_file)
        assert errors == []


class TestQualityCheckAnalysisSet:

    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_active_file_objs')
    @patch('sc_pipe_management.wrangler_utils.check_accession_results._check_if_controlled_access')
    def test_init(self, mock_get_access, mock_get_active, mock_igvf_client, mock_analysis_set_obj):
        """Test QualityCheckAnalysisSet initialization."""
        mock_response = Mock()
        mock_response.actual_instance = mock_analysis_set_obj
        mock_igvf_client.get_by_id.return_value = mock_response
        mock_get_access.return_value = False
        mock_get_active.return_value = []

        qc_checker = qa_run.QualityCheckAnalysisSet(
            'IGVFDS123456', mock_igvf_client)

        assert qc_checker.analysis_set_acc == 'IGVFDS123456'
        assert qc_checker.analysis_set_obj == mock_analysis_set_obj

    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_active_file_objs')
    @patch('sc_pipe_management.wrangler_utils.check_accession_results._check_if_controlled_access')
    def test_run_all_checks_success(self, mock_get_access, mock_get_active, mock_igvf_client, mock_analysis_set_obj):
        """Test run_all_checks with all checks passing."""
        mock_response = Mock()
        mock_response.actual_instance = mock_analysis_set_obj
        mock_igvf_client.get_by_id.return_value = mock_response
        mock_get_access.return_value = False
        mock_get_active.return_value = []

        # Mock the analysis set to have files (not empty) to avoid "has no files" error
        mock_analysis_set_obj.files = [
            '/matrix-files/test1/', '/matrix-files/test2/']

        qc_checker = qa_run.QualityCheckAnalysisSet(
            'IGVFDS123456', mock_igvf_client)

        # Mock all checkers to return empty error lists
        qc_checker.doc_checker.check_pipeline_parameters_document = Mock(
            return_value=[])
        qc_checker.rna_checker.check_rna_file_count = Mock(return_value=[])
        qc_checker.rna_checker._check_basic_file_properties = Mock(
            return_value=[])
        qc_checker.rna_checker._check_reference_files = Mock(return_value=[])
        qc_checker.qc_checker.check_rna_qc_metrics = Mock(return_value=[])
        qc_checker.rna_checker.check_rna_file_format_and_content = Mock(
            return_value=[])

        result = qc_checker.run_all_checks()
        assert result == "All checks passed."

    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_active_file_objs')
    @patch('sc_pipe_management.wrangler_utils.check_accession_results._check_if_controlled_access')
    def test_run_all_checks_with_errors(self, mock_get_access, mock_get_active, mock_igvf_client, mock_analysis_set_obj):
        """Test run_all_checks with errors."""
        mock_response = Mock()
        mock_response.actual_instance = mock_analysis_set_obj
        mock_igvf_client.get_by_id.return_value = mock_response
        mock_get_access.return_value = False
        mock_get_active.return_value = []

        qc_checker = qa_run.QualityCheckAnalysisSet(
            'IGVFDS123456', mock_igvf_client)

        # Mock the document checker to return errors
        qc_checker.doc_checker.check_pipeline_parameters_document = Mock(
            return_value=["Test error"])

        result = qc_checker.run_all_checks()

        # Check if result is a dictionary with errors
        assert isinstance(result, dict)
        assert "Error_1" in result
        assert "Test error" in result["Error_1"]


class TestMainFunction:

    @patch('sc_pipe_management.wrangler_utils.check_accession_results.QualityCheckAnalysisSet')
    def test_main_function(self, mock_qc_class, mock_igvf_client):
        """Test main function with basic functionality."""

        # Mock the QualityCheckAnalysisSet to return expected result
        mock_qc_instance = Mock()
        mock_qc_instance.run_all_checks.return_value = "All checks passed."
        mock_qc_class.return_value = mock_qc_instance

        # Create a temporary file for output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            temp_output_file = f.name

        try:
            result = qa_run.main(
                list_of_analysis_set_acc=['IGVFDS123456'],
                igvf_client_api=mock_igvf_client,
                output_file_path=temp_output_file
            )

            expected = {'IGVFDS123456': 'All checks passed.'}
            assert result == expected

            # Verify file was written
            assert os.path.exists(temp_output_file)

            # Verify file contents
            with open(temp_output_file, 'r') as f:
                file_contents = json.load(f)
            assert file_contents == expected

        finally:
            # Clean up temp file
            if os.path.exists(temp_output_file):
                os.unlink(temp_output_file)

    def test_main_function_no_output_file(self):
        """Test main function with real sandbox cases."""
        test_output_file = 'src/tests/test_files/test_qa_result.json'
        ref_output_file = 'src/tests/test_files/reference_qa_result.json'

        igvf_endpoint = 'sandbox'
        igvf_keys = api_tools.set_up_api_keys(igvf_endpoint=igvf_endpoint)
        real_igvf_client = api_tools.get_igvf_client_auth(
            igvf_keys, igvf_endpoint)

        qa_run.main(list_of_analysis_set_acc=['TSTDS33660419'],
                    igvf_client_api=real_igvf_client,
                    output_file_path=test_output_file)

        # Compare JSON content directly
        with open(ref_output_file, 'r') as f:
            expected = json.load(f)

        with open(test_output_file, 'r') as f:
            actual = json.load(f)

        assert actual == expected
