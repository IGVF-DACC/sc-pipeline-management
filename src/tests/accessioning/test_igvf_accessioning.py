import pytest
import pandas as pd
import os
import sys
from unittest.mock import Mock, patch

# Add the absolute path to the 'src' directory to sys.path
src_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

# Import the modules to test
import sc_pipe_management.accession.terra_to_portal_posting as t2p_posting


class TestTerraToPortalPosting:
    @patch('sc_pipe_management.accession.terra_to_portal_posting.terra_parse.TerraOutputMetadata')
    @patch('sc_pipe_management.accession.terra_to_portal_posting.igvf_payloads.PipelineParamsInfo')
    @patch('sc_pipe_management.accession.terra_to_portal_posting.IGVFAccessioning')
    def test_post_single_pipeline_run_integration(
        self,
        mock_IGVFAccessioning,
        mock_PipelineParamsInfo,
        mock_TerraOutputMetadata
    ):
        """Test the post_single_pipeline_run function with integration test."""

        # Prepare a mock terra_data_record (single row)
        terra_data_record = pd.Series({
            'analysis_set_acc': 'IGVFDS123ABC',
            'atac_MeaSetIDs': ['IGVFDS5678DF'],
            'rna_MeaSetIDs': ['IGVFDS9012GH'],
            'taxa': 'Homo sapiens',
            'some_other_field': 'value'
        })

        # Prepare mock pipeline_params_info
        pipeline_params_info = mock_PipelineParamsInfo.return_value

        # Prepare mock IGVF client and utils APIs
        mock_igvf_client_api = Mock()
        mock_igvf_utils_api = Mock()

        # Prepare mock TerraOutputMetadata instance
        mock_terra_metadata = mock_TerraOutputMetadata.return_value
        mock_terra_metadata.terra_data_record = terra_data_record
        mock_terra_metadata.multiome_types = ['RNAseq', 'ATACseq']
        mock_terra_metadata.anaset_accession = 'IGVFDS123ABC'

        # Prepare mock IGVFAccessioning instance and its methods
        mock_accessioning = mock_IGVFAccessioning.return_value
        mock_accessioning.post_all_rnaseq_output.return_value = [
            Mock(col_header='RNAseq', uuid='uuid1', error=None)]
        mock_accessioning.post_all_atac_alignment_output.return_value = [
            Mock(col_header='ATACseq_alignment', uuid='uuid2', error=None)]
        mock_accessioning.post_all_atac_fragment_output.return_value = [
            Mock(col_header='ATACseq_fragment', uuid='uuid3', error=None)]
        mock_accessioning.post_document.return_value = Mock(
            col_header='Document', uuid='doc-uuid', error=None)
        mock_accessioning.patch_analysis_set.return_value = Mock(
            col_header='AnalysisSetPatch', uuid='patch-uuid', error=None)

        # Call the function
        results = t2p_posting.post_single_pipeline_run(
            terra_data_record=terra_data_record,
            pipeline_params_info=pipeline_params_info,
            igvf_client_api=mock_igvf_client_api,
            igvf_utils_api=mock_igvf_utils_api,
            output_root_dir='/tmp/test_output',
            upload_file=False,
            resumed_posting=False
        )

        # Assertions
        # Should include all expected post results
        assert isinstance(results, list)
        assert any(r.col_header == 'RNAseq' for r in results)
        assert any(r.col_header == 'ATACseq_alignment' for r in results)
        assert any(r.col_header == 'ATACseq_fragment' for r in results)
        assert any(r.col_header == 'Document' for r in results)
        assert any(r.col_header == 'AnalysisSetPatch' for r in results)

        # Check that the accessioning methods were called
        mock_accessioning.post_all_rnaseq_output.assert_called_once()
        mock_accessioning.post_all_atac_alignment_output.assert_called_once()
        mock_accessioning.post_all_atac_fragment_output.assert_called_once()
        mock_accessioning.post_document.assert_called_once()
        mock_accessioning.patch_analysis_set.assert_called_once_with(
            document_uuid='doc-uuid')


if __name__ == '__main__':
    pytest.main([__file__, "-v"])
