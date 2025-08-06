import pytest
import pandas as pd
from unittest.mock import Mock
import sys
import os

# Add the absolute path to the 'src' directory to sys.path
src_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)


import sc_pipe_management.accession.parse_terra_metadata as parse_terra


class TestUtilityFunctions:
    def test_parse_terra_str_list(self):
        # Typical Terra string list
        s = ["['IGVFFI4773YQEF', 'IGVFFI7241VYRQ']"]
        result = parse_terra._parse_terra_str_list(s)
        assert result == ["IGVFFI4773YQEF", "IGVFFI7241VYRQ"]


class TestTerraMetadataParse:
    @pytest.fixture
    def mock_terra_metadata(self):
        """Mock TerraOutputMetadata instance."""
        mock_metadata = Mock()
        mock_metadata.terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC',
            'rna_kb_h5ad': 'gs://bucket/sub123/wf456/subwf789/file1.h5ad',
            'atac_MeaSetIDs': "['atac1', 'atac2']",
            'rna_MeaSetIDs': "['rna1', 'rna2']"
        })
        mock_metadata.taxa = 'Homo sapiens'
        mock_metadata.anaset_accession = 'IGVFDS123ABC'
        mock_metadata.multiome_types = ['ATACseq', 'RNAseq']

        # Mock UUIDs
        mock_uuids = Mock()
        mock_uuids.gcloud_bucket = 'bucket'
        mock_uuids.submission_id = 'sub123'
        mock_uuids.workflow_id = 'wf456'
        mock_uuids.subworkflow_id = 'subwf789'
        mock_uuids.input_param_aliases.return_value = ['sub123_wf456']
        mock_uuids.aliases.return_value = ['bucket_sub123_wf456_subwf789']
        mock_metadata._parse_workflow_uuids_from_gs_path.return_value = mock_uuids

        # Mock input file accessions
        mock_input_files = Mock()
        mock_input_files.sequence_files = ['IGVFFF001AAA', 'IGVFFF002BBB']
        mock_input_files.seqspec_files = ['IGVFFF003CCC', 'IGVFFF004DDD']
        mock_input_files.rna_barcode_replacement = 'IGVFFF005EEE'
        mock_input_files.reference_files = ['IGVFFF006FFF']

        mock_input_files.get_derived_from.return_value = [
            'IGVFFF001AAA', 'IGVFFF002BBB', 'IGVFFF003CCC', 'IGVFFF005EEE']

        mock_metadata._get_input_file_accs_from_table.return_value = mock_input_files

        return mock_metadata

    @pytest.fixture
    def mock_igvf_api(self):
        return Mock()

    def test_terra_job_uuids_init(self, mock_terra_metadata):
        uuids = mock_terra_metadata._parse_workflow_uuids_from_gs_path()
        assert uuids.gcloud_bucket == 'bucket'
        assert uuids.submission_id == 'sub123'
        assert uuids.workflow_id == 'wf456'
        assert uuids.subworkflow_id == 'subwf789'

        aliases = uuids.aliases()
        assert aliases[0] == 'bucket_sub123_wf456_subwf789'

        config_aliases = uuids.input_param_aliases()
        assert config_aliases[0] == 'sub123_wf456'

    def test_parse_workflow_uuids_from_gs_path(self, gs_path):
        uuids = parse_terra.parse_workflow_uuids_from_gs_path(gs_path)
        assert uuids.gcloud_bucket == 'bucket'
        assert uuids.submission_id == 'sub123'
        assert uuids.workflow_id == 'wf456'
        assert uuids.subworkflow_id == 'subwf789'

    def test_terra_output_metadata_init(self, mock_terra_metadata):
        assert mock_terra_metadata.terra_data_record['analysis_set_acc'] == 'IGVFDS123ABC'
        assert mock_terra_metadata.terra_data_record[
            'rna_kb_h5ad'] == 'gs://bucket/sub123/wf456/subwf789/file1.h5ad'

    def test_get_input_file_accs_from_table(self, mock_terra_metadata):
        input_files = mock_terra_metadata._get_input_file_accs_from_table()
        assert input_files.sequence_files == ['IGVFFF001AAA', 'IGVFFF002BBB']
        assert input_files.seqspec_files == ['IGVFFF003CCC', 'IGVFFF004DDD']
        assert input_files.rna_barcode_replacement == 'IGVFFF005EEE'
        assert input_files.reference_files == ['IGVFFF006FFF']
        assert sorted(input_files.get_derived_from()) == sorted([
            'IGVFFF001AAA', 'IGVFFF002BBB', 'IGVFFF003CCC', 'IGVFFF005EEE'
        ])

    def test_get_multiome_types(self, mock_terra_metadata):
        multiome_types = parse_terra._get_multiome_types(
            mock_terra_metadata.terra_data_record)
        assert sorted(multiome_types) == sorted(['ATACseq', 'RNAseq'])


if __name__ == '__main__':
    pytest.main([__file__, "-v"])
