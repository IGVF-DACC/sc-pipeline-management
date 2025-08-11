import pytest
from unittest.mock import MagicMock, patch
import os
import sys

# Add the absolute path to the 'src' directory to sys.path
src_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing
import sc_pipe_management.input_params_generation.constant as const


@pytest.fixture
def mock_igvf_api():
    """Fixture to mock the IGVF API client."""
    mock_api = MagicMock()
    return mock_api


def test_construct_full_href_url():
    """Test the construct_full_href_url function."""
    with patch.object(const, 'BASE_IGVF_PORTAL_URL', "https://example.com/"):
        igvf_href = "/files/IGVFFI123456/@@download/IGVFFI123456.fastq.gz"
        expected_url = "https://example.com/files/IGVFFI123456/@@download/IGVFFI123456.fastq.gz"
        assert portal_parsing.construct_full_href_url(
            igvf_href) == expected_url


class TestGetSeqFileMetadata:
    """Test class for GetSeqFileMetadata methods."""

    @pytest.fixture
    def mock_seqfile_obj(self):
        """Mock sequence file object."""
        mock_obj = MagicMock()
        mock_obj.accession = "IGVFFI123456"
        mock_obj.file_set = "FileSet1"
        mock_obj.illumina_read_type = "Read1"
        mock_obj.sequencing_run = 1
        mock_obj.lane = 2
        mock_obj.flowcell_id = "Flowcell123"
        mock_obj.href = "/sequence-files/IGVFFI123456/@@download/IGVFFI123456.fastq.gz"
        mock_obj.read_names = ["Read1", "Read2"]
        mock_obj.seqspecs = ["/configuration-files/SEQSPEC123/"]
        return mock_obj

    @pytest.fixture
    def mock_seqspec_obj(self):
        """Mock seqspec object."""
        mock_obj = MagicMock()
        mock_obj.href = "/configuration-files/SEQSPEC123/@@download/SEQSPEC123.yaml.gz"
        mock_obj.status = "released"
        return mock_obj

    def test_init(self, mock_igvf_api):
        """Test GetSeqFileMetadata initialization."""
        seqfile_id = "/sequence-files/IGVFFI123456/"
        instance = portal_parsing.GetSeqFileMetadata(
            seqfile_id=seqfile_id, igvf_api=mock_igvf_api)
        assert instance.seqfile_id == seqfile_id
        assert instance.igvf_api == mock_igvf_api

    def test_get_seqfile_metadata(self, mock_igvf_api, mock_seqfile_obj, mock_seqspec_obj):
        """Test get_seqfile_metadata method."""
        # Mock API behavior
        mock_igvf_api.get_by_id.side_effect = lambda x: MagicMock(
            actual_instance=mock_seqspec_obj if x == "/configuration-files/SEQSPEC123/" else mock_seqfile_obj
        )

        with patch.object(const, 'BASE_IGVF_PORTAL_URL', "https://example.com/"):
            instance = portal_parsing.GetSeqFileMetadata(
                seqfile_id="/sequence-files/IGVFFI123456/", igvf_api=mock_igvf_api)
            seqfile_metadata = instance.get_seqfile_metadata()

            assert seqfile_metadata.file_accession == "IGVFFI123456"
            assert seqfile_metadata.file_set == "FileSet1"
            assert seqfile_metadata.illumina_read_type == "Read1"
            assert seqfile_metadata.sequencing_run == 1
            assert seqfile_metadata.lane == 2
            assert seqfile_metadata.flowcell_id == "Flowcell123"
            assert seqfile_metadata.file_url == "https://example.com/sequence-files/IGVFFI123456/@@download/IGVFFI123456.fastq.gz"
            assert seqfile_metadata.seqspec_urls == [
                "https://example.com/configuration-files/SEQSPEC123/@@download/SEQSPEC123.yaml.gz"]
            assert seqfile_metadata.read_names == ["Read1", "Read2"]

    def test_get_seqfile_metadata_no_seqspecs(self, mock_igvf_api, mock_seqfile_obj):
        """Test get_seqfile_metadata method when no seqspecs are present."""
        mock_seqfile_obj.seqspecs = []
        mock_igvf_api.get_by_id.return_value.actual_instance = mock_seqfile_obj

        with patch.object(const, 'BASE_IGVF_PORTAL_URL', "https://example.com/"):
            instance = portal_parsing.GetSeqFileMetadata(
                seqfile_id="/sequence-files/IGVFFI123456/", igvf_api=mock_igvf_api)
            seqfile_metadata = instance.get_seqfile_metadata()

            assert seqfile_metadata.seqspec_urls == []

    def test_get_seqfile_metadata_deprecated_seqspecs(self, mock_igvf_api, mock_seqfile_obj, mock_seqspec_obj):
        """Test get_seqfile_metadata method when seqspecs are deprecated."""
        mock_seqspec_obj.status = "revoked"
        mock_igvf_api.get_by_id.side_effect = lambda x: MagicMock(
            actual_instance=mock_seqspec_obj if x == "/configuration-files/SEQSPEC123/" else mock_seqfile_obj
        )

        with patch.object(const, 'BASE_IGVF_PORTAL_URL', "https://example.com/"):
            instance = portal_parsing.GetSeqFileMetadata(
                seqfile_id="/sequence-files/IGVFFI123456/", igvf_api=mock_igvf_api)
            seqfile_metadata = instance.get_seqfile_metadata()

            assert seqfile_metadata.seqspec_urls == []


class TestGetMeasurementSetMetadata:
    """Test class for GetMeasurementSetMetadata methods."""

    @pytest.fixture
    def mock_measet_obj(self):
        """Mock measurement set object."""
        mock_obj = MagicMock()
        mock_obj.accession = "MEASET123"
        mock_obj.assay_term = "/assay-terms/OBI_0003109/"
        mock_obj.preferred_assay_titles = ["SHARE-seq"]
        mock_obj.files = ["/sequence-files/IGVFFI123456"]
        mock_obj.onlist_method = "no combination"
        mock_obj.onlist_files = [
            "/tabular-files/file1/"]
        mock_obj.barcode_replacement_file = "/tabular-files/barcodefile1/"
        return mock_obj

    @pytest.fixture
    def mock_seqfile_metadata(self):
        """Mock sequence file metadata."""
        return portal_parsing.SeqFileMetadata(
            file_accession="IGVFFI123456",
            file_set="FileSet1",
            illumina_read_type="Read1",
            sequencing_run=1,
            lane=2,
            flowcell_id="Flowcell123",
            file_url="https://example.com/sequence-files/IGVFFI123456/@@download/IGVFFI123456.fastq.gz",
            seqspec_urls=[
                "https://example.com/configuration-files/SEQSPEC123/@@download/SEQSPEC123.yaml.gz"],
            read_names=["Read1", "Read2"]
        )

    def test_init(self, mock_igvf_api):
        """Test GetMeasurementSetMetadata initialization."""
        measet_id = "MEASET123"
        instance = portal_parsing.GetMeasurementSetMetadata(
            measet_id=measet_id, igvf_api=mock_igvf_api)
        assert instance.measet_id == measet_id
        assert instance.igvf_api == mock_igvf_api

    def test_get_measurement_set_metadata_onlist_mapping_false(self, mock_igvf_api, mock_measet_obj, mock_seqfile_metadata):
        """Test get_measurement_set_metadata method."""
        mock_igvf_api.get_by_id.return_value.actual_instance = mock_measet_obj

        # Patch barcode replacement file's href to be a string
        mock_igvf_api.get_by_id.return_value.actual_instance.barcode_replacement_file = "/tabular-files/barcodefile1/"
        mock_igvf_api.get_by_id.return_value.actual_instance.href = "/tabular-files/barcodefile1/@@download/barcodefile1"

        with patch.object(const, 'ASSAY_NAMES_CONVERSION_REF', {"/assay-terms/OBI_0003109/": "rna"}):
            with patch('sc_pipe_management.input_params_generation.portal_metadata_parsing.GetSeqFileMetadata') as mock_get_seqfile:
                mock_get_seqfile.return_value.get_seqfile_metadata.return_value = mock_seqfile_metadata

                instance = portal_parsing.GetMeasurementSetMetadata(
                    measet_id="MEASET123", igvf_api=mock_igvf_api)
                measet_metadata = instance.get_measurement_set_metadata()

                assert measet_metadata.measet_acc == "MEASET123"
                assert measet_metadata.assay_type == "rna"
                assert measet_metadata.onlist_mapping is False
                assert len(measet_metadata.seqfiles) == 1
                assert measet_metadata.seqfiles[0].file_accession == "IGVFFI123456"
                assert measet_metadata.onlist_method == "no combination"
                assert measet_metadata.onlist_files == [
                    "/tabular-files/file1/"]
                assert measet_metadata.barcode_replacement_file == "https://api.data.igvf.org/tabular-files/barcodefile1/@@download/barcodefile1"

    def test_get_measurement_set_metadata_onlist_mapping_true(self, mock_igvf_api, mock_measet_obj, mock_seqfile_metadata):
        """Test get_measurement_set_metadata method when onlist_mapping is True."""
        mock_measet_obj.preferred_assay_titles = ['10x multiome']
        mock_igvf_api.get_by_id.return_value.actual_instance = mock_measet_obj

        # Patch barcode replacement file's href to be a string
        mock_igvf_api.get_by_id.return_value.actual_instance.barcode_replacement_file = "/tabular-files/barcodefile1/"
        mock_igvf_api.get_by_id.return_value.actual_instance.href = "/tabular-files/barcodefile1/"

        with patch.object(const, 'ASSAY_NAMES_CONVERSION_REF', {"/assay-terms/OBI_0003109/": "rna"}):
            with patch('sc_pipe_management.input_params_generation.portal_metadata_parsing.GetSeqFileMetadata') as mock_get_seqfile:
                mock_get_seqfile.return_value.get_seqfile_metadata.return_value = mock_seqfile_metadata

                instance = portal_parsing.GetMeasurementSetMetadata(
                    measet_id="MEASET123", igvf_api=mock_igvf_api)
                measet_metadata = instance.get_measurement_set_metadata()

                assert measet_metadata.onlist_mapping is True

    def test_get_measurement_set_metadata_no_barcode_replacement(self, mock_igvf_api, mock_measet_obj, mock_seqfile_metadata):
        """Test get_measurement_set_metadata method when barcode_replacement_file is not present."""
        mock_measet_obj.barcode_replacement_file = None
        mock_igvf_api.get_by_id.return_value.actual_instance = mock_measet_obj

        with patch.object(const, 'ASSAY_NAMES_CONVERSION_REF', {"/assay-terms/OBI_0003109/": "rna"}):
            with patch('sc_pipe_management.input_params_generation.portal_metadata_parsing.GetSeqFileMetadata') as mock_get_seqfile:
                mock_get_seqfile.return_value.get_seqfile_metadata.return_value = mock_seqfile_metadata

                instance = portal_parsing.GetMeasurementSetMetadata(
                    measet_id="MEASET123", igvf_api=mock_igvf_api)
                measet_metadata = instance.get_measurement_set_metadata()

                assert measet_metadata.barcode_replacement_file is None


class TestGetAnalysisSetMetadata:
    """Test class for GetAnalysisSetMetadata methods."""

    @pytest.fixture
    def mock_sample_obj(self):
        """Mock sample object."""
        mock_obj = MagicMock()
        mock_obj.taxa = 'Homo sapiens'
        mock_obj.accession = 'SAMPLE123'
        return mock_obj

    @pytest.fixture
    def mock_analysis_set_obj(self):
        """Mock analysis set object."""
        mock_obj = MagicMock()
        mock_obj.samples = ['/samples/SAMPLE123/']
        mock_obj.input_file_sets = [
            "/measurement-sets/MEASET123/", "/measurement-sets/MEASET456/"]
        return mock_obj

    @pytest.fixture
    def mock_rna_measet_metadata(self):
        """Mock RNA measurement set metadata."""
        return portal_parsing.MeasurementSetMetadata(
            measet_acc="MEASET123",
            assay_type="rna",
            onlist_mapping=False,
            seqfiles=[],
            onlist_method="no combination",
            onlist_files=["/tabular-files/file1/"],
            barcode_replacement_file=None
        )

    @pytest.fixture
    def mock_atac_measet_metadata(self):
        """Mock ATAC measurement set metadata."""
        return portal_parsing.MeasurementSetMetadata(
            measet_acc="MEASET456",
            assay_type="atac",
            onlist_mapping=False,
            seqfiles=[],
            onlist_method="product",
            onlist_files=["/tabular-files/file3/", "/tabular-files/file4/"],
            barcode_replacement_file=None
        )

    def test_init(self, mock_igvf_api):
        """Test GetAnalysisSetMetadata initialization."""
        analysis_set_accession = "ANASET123"
        instance = portal_parsing.GetAnalysisSetMetadata(
            analysis_set_accession=analysis_set_accession, igvf_api=mock_igvf_api)
        assert instance.analysis_set_accession == analysis_set_accession
        assert instance.igvf_api == mock_igvf_api

    def test_get_input_analysis_set_metadata(self, mock_igvf_api, mock_analysis_set_obj, mock_sample_obj,
                                             mock_rna_measet_metadata, mock_atac_measet_metadata):
        """Test get_input_analysis_set_metadata method."""
        def mock_get_by_id_side_effect(obj_id):
            if obj_id == "/samples/SAMPLE123/":
                return MagicMock(actual_instance=mock_sample_obj)
            else:
                return MagicMock(actual_instance=mock_analysis_set_obj)

        mock_igvf_api.get_by_id.side_effect = mock_get_by_id_side_effect

        def mock_get_measurement_set_metadata_side_effect(measet_id, igvf_api):
            if measet_id == "/measurement-sets/MEASET123/":
                mock_instance = MagicMock()
                mock_instance.get_measurement_set_metadata.return_value = mock_rna_measet_metadata
                return mock_instance
            elif measet_id == "/measurement-sets/MEASET456/":
                mock_instance = MagicMock()
                mock_instance.get_measurement_set_metadata.return_value = mock_atac_measet_metadata
                return mock_instance

        with patch('sc_pipe_management.input_params_generation.portal_metadata_parsing.GetMeasurementSetMetadata') as mock_get_measet:
            mock_get_measet.side_effect = mock_get_measurement_set_metadata_side_effect

            instance = portal_parsing.GetAnalysisSetMetadata(
                analysis_set_accession="ANASET123", igvf_api=mock_igvf_api)
            analysis_set_metadata = instance.get_input_analysis_set_metadata()

            assert analysis_set_metadata.analysis_set_acc == "ANASET123"

            assert len(analysis_set_metadata.sample_info) == 1
            assert analysis_set_metadata.sample_info[0].taxa == "Homo sapiens"
            assert analysis_set_metadata.sample_info[0].subpool_id == "SAMPLE123"

            assert len(analysis_set_metadata.rna_input_info) == 1
            assert analysis_set_metadata.rna_input_info[0].measet_acc == "MEASET123"

            assert len(analysis_set_metadata.atac_input_info) == 1
            assert analysis_set_metadata.atac_input_info[0].measet_acc == "MEASET456"

    def test_get_input_analysis_set_metadata_only_one_assay_type(self, mock_igvf_api, mock_analysis_set_obj,
                                                                 mock_sample_obj, mock_rna_measet_metadata, mock_atac_measet_metadata):
        """Test get_input_analysis_set_metadata method with only RNA measurement sets."""
        mock_analysis_set_obj.input_file_sets = ["/measurement-sets/MEASET123"]

        def mock_get_by_id_side_effect(obj_id):
            if obj_id == "/samples/SAMPLE123/":
                return MagicMock(actual_instance=mock_sample_obj)
            else:
                return MagicMock(actual_instance=mock_analysis_set_obj)

        mock_igvf_api.get_by_id.side_effect = mock_get_by_id_side_effect

        # Test if only RNA measurement set metadata is returned
        with patch('sc_pipe_management.input_params_generation.portal_metadata_parsing.GetMeasurementSetMetadata') as mock_get_measet:
            mock_instance = MagicMock()
            mock_instance.get_measurement_set_metadata.return_value = mock_rna_measet_metadata
            mock_get_measet.return_value = mock_instance

            instance = portal_parsing.GetAnalysisSetMetadata(
                analysis_set_accession="ANASET123", igvf_api=mock_igvf_api)
            analysis_set_metadata = instance.get_input_analysis_set_metadata()

            assert len(analysis_set_metadata.rna_input_info) == 1
            assert len(analysis_set_metadata.atac_input_info) == 0

        # Test if only ATAC measurement set metadata is returned
        with patch('sc_pipe_management.input_params_generation.portal_metadata_parsing.GetMeasurementSetMetadata') as mock_get_measet:
            mock_instance = MagicMock()
            mock_instance.get_measurement_set_metadata.return_value = mock_atac_measet_metadata
            mock_get_measet.return_value = mock_instance

            instance = portal_parsing.GetAnalysisSetMetadata(
                analysis_set_accession="ANASET123", igvf_api=mock_igvf_api)
            analysis_set_metadata = instance.get_input_analysis_set_metadata()

            assert len(analysis_set_metadata.rna_input_info) == 0
            assert len(analysis_set_metadata.atac_input_info) == 1
