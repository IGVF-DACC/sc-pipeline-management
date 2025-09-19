import pytest
from unittest.mock import Mock, MagicMock, patch
import json
import tempfile
import os
import sys


# Add the absolute path to the 'src' directory to sys.path
src_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.wrangler_utils.check_accession_results as qa_run
import sc_pipe_management.igvf_and_terra_api_tools as api_tools


@pytest.fixture
def mock_igvf_client():
    """Create a mock IGVF client."""
    client = Mock()
    return client


@pytest.fixture
def mock_file_obj():
    """Create a mock file object."""
    file_obj = Mock()
    file_obj.type = ['MatrixFile']
    file_obj.accession = 'IGVFFI123456'
    file_obj.analysis_step_version = '/analysis-step-versions/test-version/'
    file_obj.upload_status = 'validated'
    file_obj.controlled_access = False
    file_obj.reference_files = ['/reference-files/test-ref/']
    file_obj.submitted_file_name = 'test_file.h5ad'
    file_obj.status = 'released'
    file_obj.quality_metrics = ['/quality-metrics/test-qc/']
    file_obj.input_file_for = []
    file_obj.file_format = 'h5ad'
    file_obj.content_type = 'sparse gene count matrix'
    return file_obj


@pytest.fixture
def mock_analysis_set_obj():
    """Create a mock analysis set object."""
    analysis_set = Mock()
    analysis_set.files = ['/matrix-files/test1/', '/alignment-files/test2/']
    analysis_set.pipeline_parameters = ['/documents/test-doc/']
    analysis_set.input_file_sets = ['/measurement-sets/test-set/']
    return analysis_set


class TestUtilityFunctions:
    """Test utility functions."""

    def test_get_active_file_objs_empty_list(self, mock_igvf_client):
        """Test _get_active_file_objs with empty file_ids."""
        result = qa_run._get_active_file_objs([], mock_igvf_client)
        assert result == []

    def test_get_active_file_objs_with_active_files(self, mock_igvf_client, mock_file_obj):
        """Test _get_active_file_objs with active files."""
        mock_response = Mock()
        mock_response.actual_instance = mock_file_obj
        mock_igvf_client.get_by_id.return_value = mock_response

        result = qa_run._get_active_file_objs(['file1'], mock_igvf_client)
        assert len(result) == 1
        assert result[0] == mock_file_obj

    def test_get_active_file_objs_filters_inactive_files(self, mock_igvf_client):
        """Test _get_active_file_objs filters out inactive files."""
        inactive_file = Mock()
        inactive_file.status = 'deleted'

        mock_response = Mock()
        mock_response.actual_instance = inactive_file
        mock_igvf_client.get_by_id.return_value = mock_response

        result = qa_run._get_active_file_objs(['file1'], mock_igvf_client)
        assert result == []

    def test_get_raw_data_access_level_true(self, mock_igvf_client, mock_analysis_set_obj):
        """Test _get_raw_data_access_level returns True when controlled access."""
        mock_seq_file = Mock()
        mock_seq_file.controlled_access = True

        mock_file_set = Mock()
        mock_file_set.files = ['/sequence-files/test-seq/']

        mock_responses = {
            '/measurement-sets/test-set/': Mock(actual_instance=mock_file_set),
            '/sequence-files/test-seq/': Mock(actual_instance=mock_seq_file)
        }
        mock_igvf_client.get_by_id.side_effect = lambda x: mock_responses[x]

        result = qa_run._get_raw_data_access_level(
            mock_analysis_set_obj, mock_igvf_client)
        assert result is True

    def test_get_raw_data_access_level_mixed_controlled_access(self, mock_igvf_client, mock_analysis_set_obj):
        """Test _get_raw_data_access_level returns True when at least one sequence file has controlled access."""
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

        result = qa_run._get_raw_data_access_level(
            mock_analysis_set_obj, mock_igvf_client)

        # Should return True because at least one sequence file has controlled_access=True
        assert result is True

    def test_is_sublist_in_nested_simple_list(self):
        """Test _is_sublist_in_nested with simple list."""
        sublist = ['a', 'b']
        main_list = ['a', 'b', 'c', 'd']
        assert qa_run._is_sublist_in_nested(sublist, main_list) is True

    def test_is_sublist_in_nested_nested_list(self):
        """Test _is_sublist_in_nested with nested list."""
        sublist = ['a', 'b']
        main_list = [['a', 'b', 'c'], ['d', 'e', 'f']]
        assert qa_run._is_sublist_in_nested(sublist, main_list) is True

    def test_is_sublist_in_nested_not_found(self):
        """Test _is_sublist_in_nested when sublist not found."""
        sublist = ['x', 'y']
        main_list = [['a', 'b', 'c'], ['d', 'e', 'f']]
        assert qa_run._is_sublist_in_nested(sublist, main_list) is False


class TestBaseFileChecker:
    """Test BaseFileChecker class."""

    def test_init(self, mock_igvf_client):
        """Test BaseFileChecker initialization."""
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)
        assert checker.analysis_set_acc == 'IGVFDS123456'
        assert checker.igvf_client_api == mock_igvf_client

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
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)
        reference_list = ['/reference-files/test-ref/',
                          '/reference-files/other-ref/']
        # Should not raise any exception
        checker._check_reference_files(mock_file_obj, reference_list)

    def test_check_reference_files_none(self, mock_igvf_client, mock_file_obj):
        """Test _check_reference_files with None reference files."""
        mock_file_obj.reference_files = None
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="has no reference files"):
            checker._check_reference_files(
                mock_file_obj, ['/reference-files/test-ref/'])

    def test_check_file_name_success(self, mock_igvf_client, mock_file_obj):
        """Test _check_file_name with valid file name."""
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker._check_file_name(mock_file_obj, 'test_file')

    def test_check_file_name_invalid(self, mock_igvf_client, mock_file_obj):
        """Test _check_file_name with invalid file name."""
        checker = qa_run.BaseFileChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="has unexpected file name"):
            checker._check_file_name(mock_file_obj, 'wrong_name')


class TestDocumentChecker:
    """Test DocumentChecker class."""

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
        # Should not raise any exception
        checker.check_pipeline_parameters_document(['/documents/test-doc/'])

    def test_check_pipeline_parameters_document_empty(self, mock_igvf_client):
        """Test check_pipeline_parameters_document with empty list."""
        checker = qa_run.DocumentChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="No pipeline parameters document found"):
            checker.check_pipeline_parameters_document([])

    def test_check_pipeline_parameters_document_wrong_type(self, mock_igvf_client):
        """Test check_pipeline_parameters_document with wrong document type."""
        mock_doc = Mock()
        mock_doc.document_type = 'other document'

        # Create a proper Attachment mock object (instead of simple dict)
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

        with pytest.raises(ValueError, match="Unexpected pipeline parameters document type"):
            checker.check_pipeline_parameters_document(
                ['/documents/test-doc/'])

    def test_check_pipeline_parameters_document_no_attachment(self, mock_igvf_client):
        """Test check_pipeline_parameters_document with no attachment."""
        mock_doc = Mock()
        mock_doc.document_type = 'pipeline parameters'
        mock_doc.attachment = None

        mock_response = Mock()
        mock_response.actual_instance = mock_doc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.DocumentChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="has no attachment"):
            checker.check_pipeline_parameters_document(
                ['/documents/test-doc/'])


class TestQCMetricsChecker:
    """Test QCMetricsChecker class."""

    def test_check_qc_metrics_basic_success(self, mock_igvf_client, mock_file_obj):
        """Test check_qc_metrics_basic with valid QC metrics."""
        mock_qc = Mock()
        mock_qc.analysis_step_version = '/analysis-step-versions/test-version/'

        mock_response = Mock()
        mock_response.actual_instance = mock_qc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.QCMetricsChecker('IGVFDS123456', mock_igvf_client)
        result = checker.check_qc_metrics_basic(mock_file_obj)
        assert result == mock_qc

    def test_check_qc_metrics_basic_no_metrics(self, mock_igvf_client, mock_file_obj):
        """Test check_qc_metrics_basic with no quality metrics."""
        mock_file_obj.quality_metrics = None
        checker = qa_run.QCMetricsChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="has no quality metrics"):
            checker.check_qc_metrics_basic(mock_file_obj)

    def test_check_qc_metrics_basic_multiple_metrics(self, mock_igvf_client, mock_file_obj):
        """Test check_qc_metrics_basic with multiple quality metrics."""
        mock_file_obj.quality_metrics = ['/qc1/', '/qc2/']
        checker = qa_run.QCMetricsChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="has more than 1 quality metrics"):
            checker.check_qc_metrics_basic(mock_file_obj)

    def test_check_rna_qc_metrics_success(self, mock_igvf_client, mock_file_obj):
        """Test check_rna_qc_metrics with valid RNA QC metrics."""
        mock_qc = Mock()
        mock_qc.analysis_step_version = '/analysis-step-versions/test-version/'
        mock_qc.type = ['SingleCellRnaSeqQualityMetric',
                        'QualityMetrics', 'Item']
        mock_qc.rnaseq_kb_info = {'href': '/test-attachment/'}
        mock_qc.n_reads = 1000000

        mock_response = Mock()
        mock_response.actual_instance = mock_qc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.QCMetricsChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker.check_rna_qc_metrics(mock_file_obj)

    def test_check_atac_qc_metrics_alignment_file(self, mock_igvf_client, mock_file_obj):
        """Test check_atac_qc_metrics with AlignmentFile."""
        mock_file_obj.type = ['AlignmentFile', 'File', 'Item']
        mock_qc = Mock()
        mock_qc.analysis_step_version = '/analysis-step-versions/test-version/'
        mock_qc.type = ['SingleCellAtacSeqQualityMetric',
                        'QualityMetrics', 'Item']
        mock_qc.atac_bam_summary_stats = Mock()
        mock_qc.atac_bam_summary_stats.href = '@@download/atac_bam_summary_stats/qc.barcode.summary.csv'

        mock_response = Mock()
        mock_response.actual_instance = mock_qc
        mock_igvf_client.get_by_id.return_value = mock_response

        checker = qa_run.QCMetricsChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker.check_atac_qc_metrics(mock_file_obj)


class TestIndexFileChecker:
    """Test IndexFileChecker class."""

    def test_get_index_file_ids(self, mock_igvf_client, mock_file_obj):
        """Test _get_index_file_ids."""
        mock_file_obj.input_file_for = [
            '/index-files/test1/', '/other-files/test2/']
        checker = qa_run.IndexFileChecker('IGVFDS123456', mock_igvf_client)

        result = checker._get_index_file_ids(mock_file_obj)
        assert result == ['/index-files/test1/']

    def test_validate_index_file_count_success(self, mock_igvf_client, mock_file_obj):
        """Test _validate_index_file_count with exactly one index file."""
        checker = qa_run.IndexFileChecker('IGVFDS123456', mock_igvf_client)
        mock_index = Mock()
        # Should not raise any exception
        checker._validate_index_file_count(mock_file_obj, [mock_index])

    def test_validate_index_file_count_no_files(self, mock_igvf_client, mock_file_obj):
        """Test _validate_index_file_count with no index files."""
        checker = qa_run.IndexFileChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="is missing an index file"):
            checker._validate_index_file_count(mock_file_obj, [])

    def test_validate_index_file_count_multiple_files(self, mock_igvf_client, mock_file_obj):
        """Test _validate_index_file_count with multiple index files."""
        checker = qa_run.IndexFileChecker('IGVFDS123456', mock_igvf_client)
        mock_index1, mock_index2 = Mock(), Mock()

        with pytest.raises(ValueError, match="has more than 1 active index file"):
            checker._validate_index_file_count(
                mock_file_obj, [mock_index1, mock_index2])

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

        mock_get_active.return_value = [mock_index]

        checker = qa_run.IndexFileChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker.check_index_file_presence_and_format(mock_file_obj)


class TestRNASeqFileChecker:
    """Test RNASeqFileChecker class."""

    def test_check_rna_file_count_success(self, mock_igvf_client):
        """Test check_rna_file_count with exactly 2 files."""
        checker = qa_run.RNASeqFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock(), Mock()]
        # Should not raise any exception
        checker.check_rna_file_count(mock_files)

    def test_check_rna_file_count_wrong_count(self, mock_igvf_client):
        """Test check_rna_file_count with wrong number of files."""
        checker = qa_run.RNASeqFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock()]

        with pytest.raises(ValueError, match="does not have exactly 2 non-deprecated RNA-seq"):
            checker.check_rna_file_count(mock_files)

    def test_check_rna_file_format_and_content_h5ad(self, mock_igvf_client):
        """Test check_rna_file_format_and_content with h5ad format."""
        mock_file = Mock()
        mock_file.file_format = 'h5ad'
        mock_file.content_type = 'sparse gene count matrix'
        mock_file.submitted_file_name = 'IGVFDS123456.h5ad'

        checker = qa_run.RNASeqFileChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker.check_rna_file_format_and_content(mock_file)

    def test_check_rna_file_format_and_content_tar(self, mock_igvf_client):
        """Test check_rna_file_format_and_content with tar format."""
        mock_file = Mock()
        mock_file.file_format = 'tar'
        mock_file.content_type = 'kallisto single cell RNAseq output'
        mock_file.submitted_file_name = 'IGVFDS123456.tar.gz'

        checker = qa_run.RNASeqFileChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker.check_rna_file_format_and_content(mock_file)


class TestATACFileChecker:
    """Test ATACFileChecker class."""

    def test_check_atac_file_count_success(self, mock_igvf_client):
        """Test check_atac_file_count with exactly 1 file."""
        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock()]
        # Should not raise any exception
        checker.check_atac_file_count(mock_files, "AlignmentFile")

    def test_check_atac_file_count_wrong_count(self, mock_igvf_client):
        """Test check_atac_file_count with wrong number of files."""
        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)
        mock_files = [Mock(), Mock()]

        with pytest.raises(ValueError, match="has more than 1 AlignmentFile"):
            checker.check_atac_file_count(mock_files, "AlignmentFile")

    def test_check_alignment_file_specifics(self, mock_igvf_client):
        """Test check_alignment_file_specifics with valid alignment file."""
        mock_file = Mock()
        mock_file.file_format = 'bam'
        mock_file.content_type = 'alignments'
        mock_file.submitted_file_name = 'IGVFDS123456.bam'

        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker.check_alignment_file_specifics(mock_file)

    def test_check_fragment_file_specifics(self, mock_igvf_client):
        """Test check_fragment_file_specifics with valid fragment file."""
        mock_file = Mock()
        mock_file.file_format = 'bed'
        mock_file.file_format_type = 'bed3+'
        mock_file.content_type = 'fragments'
        mock_file.submitted_file_name = 'IGVFDS123456.fragments.tsv.gz'

        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)
        # Should not raise any exception
        checker.check_fragment_file_specifics(mock_file)

    def test_check_fragment_controlled_access(self, mock_igvf_client):
        """Test check_fragment_controlled_access with invalid controlled access."""
        mock_file = Mock()
        mock_file.controlled_access = True
        mock_file.accession = 'IGVFFI123456'

        checker = qa_run.ATACFileChecker('IGVFDS123456', mock_igvf_client)

        with pytest.raises(ValueError, match="should not be controlled_access files"):
            checker.check_fragment_controlled_access(mock_file)


class TestQualityCheckAnalysisSet:
    """Test QualityCheckAnalysisSet class."""

    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_active_file_objs')
    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_raw_data_access_level')
    def test_init(self, mock_get_access, mock_get_active, mock_igvf_client, mock_analysis_set_obj):
        """Test QualityCheckAnalysisSet initialization."""
        mock_response = Mock()
        mock_response.actual_instance = mock_analysis_set_obj
        mock_igvf_client.get_by_id.return_value = mock_response
        mock_get_access.return_value = False

        qc_checker = qa_run.QualityCheckAnalysisSet(
            'IGVFDS123456', mock_igvf_client)

        assert qc_checker.analysis_set_acc == 'IGVFDS123456'
        assert qc_checker.igvf_client_api == mock_igvf_client
        assert qc_checker.analysis_set_obj == mock_analysis_set_obj
        assert qc_checker.controlled_access is False

    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_active_file_objs')
    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_raw_data_access_level')
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

        # Mock all checkers to not raise any errors
        qc_checker.doc_checker.check_pipeline_parameters_document = Mock()
        qc_checker.rna_checker.check_rna_file_count = Mock()
        qc_checker.rna_checker._check_basic_file_properties = Mock()
        qc_checker.rna_checker.check_rna_controlled_access = Mock()
        qc_checker.rna_checker._check_reference_files = Mock()
        qc_checker.qc_checker.check_rna_qc_metrics = Mock()
        qc_checker.rna_checker.check_rna_file_format_and_content = Mock()

        result = qc_checker.run_all_checks()
        assert result == "All checks passed."

    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_active_file_objs')
    @patch('sc_pipe_management.wrangler_utils.check_accession_results._get_raw_data_access_level')
    def test_run_all_checks_with_errors(self, mock_get_access, mock_get_active, mock_igvf_client, mock_analysis_set_obj):
        """Test run_all_checks with errors."""
        mock_response = Mock()
        mock_response.actual_instance = mock_analysis_set_obj
        mock_igvf_client.get_by_id.return_value = mock_response
        mock_get_access.return_value = False
        mock_get_active.return_value = []

        qc_checker = qa_run.QualityCheckAnalysisSet(
            'IGVFDS123456', mock_igvf_client)

        # Mock the document checker to raise an error
        qc_checker.doc_checker.check_pipeline_parameters_document = Mock(
            side_effect=ValueError("Test error")
        )

        result = qc_checker.run_all_checks()

        # Check if result is a dictionary with errors
        if isinstance(result, dict):
            # Check if "Test error" is in any of the error values
            assert any("Test error" in str(error) for error in result.values())
        else:
            # If it's still a string format
            assert "Test error" in result


class TestMainFunction:
    """Test main function."""

    @patch('sc_pipe_management.wrangler_utils.check_accession_results.QualityCheckAnalysisSet')
    def test_main_function(self, mock_qc_class, mock_igvf_client):
        """Test main function."""
        # Mock QualityCheckAnalysisSet instance
        mock_qc_instance = Mock()
        mock_qc_instance.run_all_checks.return_value = "All checks passed."
        mock_qc_class.return_value = mock_qc_instance

        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            temp_file = f.name

        try:
            result = qa_run.main(['IGVFDS123456'], mock_igvf_client, temp_file)

            # Check results
            assert result == {'IGVFDS123456': 'All checks passed.'}

            # Check that file was written
            assert os.path.exists(temp_file)
            with open(temp_file, 'r') as f:
                data = json.load(f)
            assert data == {'IGVFDS123456': 'All checks passed.'}

        finally:
            # Clean up
            if os.path.exists(temp_file):
                os.unlink(temp_file)

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


if __name__ == "__main__":
    pytest.main([__file__])
