import os
import sys
from unittest.mock import Mock, patch
import unittest
import pytest
import dataclasses

# Add the absolute path to the 'src' directory to sys.path
src_path = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.wrangler_utils.analysis_set_setup_utils as utils


class TestAnalysisSetUtils(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.mock_igvf_client = Mock()
        self.mock_igvf_utils = Mock()

        # Sample test data
        self.test_input_file_sets = [
            '/measurement-sets/IGVFDS123456/', '/measurement-sets/IGVFDS789012/']
        self.test_lab = '/labs/j-michael-cherry/'
        self.test_award = '/awards/HG012012/'


class TestCheckIsScpipeline(TestAnalysisSetUtils):
    def test_check_is_scpipeline_not_analysis_set(self):
        """Test that non-analysis-set ID returns None."""
        # Assertion check: not single cell
        result = utils.check_is_scpipeline(
            '/measurement-sets/IGVFDS123456/', self.mock_igvf_client)
        self.assertFalse(result)

    def test_check_is_scpipeline_with_scpipe_alias(self):
        """Test that file set with scpipe alias returns True."""
        mock_file_set_obj = Mock()
        mock_file_set_obj.aliases = ['igvf:test-scpipe-analysis']
        mock_file_set_obj.description = None
        mock_file_set_obj.workflows = []

        # Mock the HTTP response object
        mock_response = Mock()
        mock_response.actual_instance = mock_file_set_obj

        # Set what get_by_id() returns
        self.mock_igvf_client.get_by_id.return_value = mock_response

        # Assertion check: single cell
        result = utils.check_is_scpipeline(
            '/analysis-sets/IGVFAS123456/', self.mock_igvf_client)
        self.assertTrue(result)

    def test_check_is_scpipeline_with_single_cell_description(self):
        """Test that file set with single cell description returns True."""
        mock_file_set_obj = Mock()
        mock_file_set_obj.aliases = []
        mock_file_set_obj.description = 'single cell uniform pipeline analysis'
        mock_file_set_obj.workflows = []

        # Mock the HTTP response object
        mock_response = Mock()
        mock_response.actual_instance = mock_file_set_obj

        # Set what get_by_id() returns
        self.mock_igvf_client.get_by_id.return_value = mock_response

        # Assertion check: single cell
        result = utils.check_is_scpipeline(
            '/analysis-sets/IGVFAS123456/', self.mock_igvf_client)
        self.assertTrue(result)

    def test_check_is_scpipeline_with_sc_pipeline_alias(self):
        """Test that file set with sc-pipeline alias returns True."""
        mock_file_set_obj = Mock()
        mock_file_set_obj.aliases = ['igvf:my-sc-pipeline-test']
        mock_file_set_obj.description = None
        mock_file_set_obj.workflows = []

        # Mock the HTTP response object
        mock_response = Mock()
        mock_response.actual_instance = mock_file_set_obj

        # Set what get_by_id() returns
        self.mock_igvf_client.get_by_id.return_value = mock_response

        # Assertion check: single cell
        result = utils.check_is_scpipeline(
            '/analysis-sets/IGVFAS123456/', self.mock_igvf_client)
        self.assertTrue(result)

    def test_check_is_scpipeline_false_case(self):
        """Test that file set without pipeline keywords returns False."""
        mock_file_set_obj = Mock()
        mock_file_set_obj.aliases = [
            'igvf:other-analysis', 'igvf:non-pipeline-test']
        mock_file_set_obj.description = 'some other analysis description'
        mock_file_set_obj.workflows = []
        mock_file_set_obj.uniform_pipeline_status = None

        # Mock the HTTP response object
        mock_response = Mock()
        mock_response.actual_instance = mock_file_set_obj

        # Set what get_by_id() returns
        self.mock_igvf_client.get_by_id.return_value = mock_response

        # Assertion check: no single cell
        result = utils.check_is_scpipeline(
            '/analysis-sets/IGVFAS123456/', self.mock_igvf_client)
        self.assertFalse(result)

    def test_check_is_scpipeline_no_aliases_no_description(self):
        """Test that file set with no aliases or description returns False."""
        mock_file_set_obj = Mock()
        mock_file_set_obj.aliases = []
        mock_file_set_obj.description = None
        mock_file_set_obj.workflows = []
        mock_file_set_obj.uniform_pipeline_status = None

        # Mock the HTTP response object
        mock_response = Mock()
        mock_response.actual_instance = mock_file_set_obj

        # Set what get_by_id() returns
        self.mock_igvf_client.get_by_id.return_value = mock_response

        # Assertion check: no single cell
        result = utils.check_is_scpipeline(
            '/analysis-sets/IGVFAS123456/', self.mock_igvf_client)
        self.assertFalse(result)

    def test_check_is_scpipeline_with_uniform_workflows(self):
        """Test that file set with workflows returns True (i.e., has files)."""
        # Create a mock workflow object with uniform_pipeline property
        mock_workflow = Mock()
        mock_workflow.uniform_pipeline = True
        mock_workflow.id = 'IGVFWF7365VWQV'

        # Set workflows to contain the mock workflow object
        mock_file_set_obj = Mock()
        mock_file_set_obj.workflows = [mock_workflow.id]
        mock_file_set_obj.uniform_pipeline_status = None

        # Mock the HTTP response object
        mock_response = Mock()
        mock_response.actual_instance = mock_file_set_obj

        # Set what get_by_id() returns
        self.mock_igvf_client.get_by_id.return_value = mock_response

        # Assertion check: single cell
        result = utils.check_is_scpipeline(
            '/analysis-sets/IGVFAS123456/', self.mock_igvf_client)
        self.assertTrue(result)

    def test_check_is_scpipeline_with_pipeline_status(self):
        """Test that file set with pipeline status returns True."""
        mock_file_set_obj = Mock()
        mock_file_set_obj.aliases = []
        mock_file_set_obj.description = None
        mock_file_set_obj.workflows = []
        mock_file_set_obj.uniform_pipeline_status = 'processing'

        # Mock the HTTP response object
        mock_response = Mock()
        mock_response.actual_instance = mock_file_set_obj

        # Set what get_by_id() returns
        self.mock_igvf_client.get_by_id.return_value = mock_response

        # Assertion check: single cell
        result = utils.check_is_scpipeline(
            '/analysis-sets/IGVFAS123456/', self.mock_igvf_client)
        self.assertTrue(result)


class TestCalcInputFileSets(TestAnalysisSetUtils):

    def test_calc_input_file_sets_single_basic(self):
        """Test basic input file set calculation."""
        # Mock query result
        mock_query_res = Mock()
        mock_query_res.id = '/measurement-sets/IGVFDS123456/'
        mock_query_res.related_multiome_datasets = [
            '/measurement-sets/IGVFDS789012/']
        mock_query_res.barcode_replacement_file = None

        # Mock check_is_duped_for_all to return False
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        expected = ['/measurement-sets/IGVFDS123456/',
                    '/measurement-sets/IGVFDS789012/']
        self.assertEqual(sorted(result), sorted(expected))

    def test_calc_input_file_sets_single_parse(self):
        """Test basic input file set calculation of Parse SPLiT-seq."""
        # Mock query result
        mock_query_res = Mock()
        mock_query_res.id = '/measurement-sets/IGVFDS123456/'
        mock_query_res.related_multiome_datasets = [
            '/measurement-sets/IGVFDS789012/']
        mock_query_res.barcode_replacement_file = '/tabular-files/IGVFFI111111/'

        # Mock barcode replacement file response
        mock_brf_response = Mock()
        mock_brf_obj = Mock()
        mock_brf_obj.file_set = '/curated-sets/IGVFDS333333/'
        mock_brf_response.actual_instance = mock_brf_obj
        self.mock_igvf_client.get_by_id.return_value = mock_brf_response

        # Mock check_is_duped_for_all to return False
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        expected = ['/measurement-sets/IGVFDS123456/',
                    '/curated-sets/IGVFDS333333/', '/measurement-sets/IGVFDS789012/']
        self.assertEqual(sorted(result), sorted(expected))

    def test_calc_input_file_sets_single_no_related_datasets(self):
        """Test input file set calculation without related datasets."""
        mock_query_res = Mock()
        mock_query_res.id = '/measurement-sets/IGVFDS123456/'
        mock_query_res.related_multiome_datasets = None
        mock_query_res.barcode_replacement_file = None

        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        expected = ['/measurement-sets/IGVFDS123456/']
        self.assertEqual(result, expected)

    def test_calc_input_file_sets_single_is_duplicate(self):
        """Test that duplicate analysis set returns None."""
        mock_query_res = Mock()
        mock_query_res.id = '/measurement-sets/IGVFDS123456/'

        with patch.object(utils, 'check_is_duped_for_all', return_value=True):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        self.assertIsNone(result)


# REMOVE the entire TestGetSampleAccessions class since get_sample_accessions function no longer exists


class TestGetSampleMetadata:
    """Test the GetSampleMetadata class."""

    @pytest.fixture
    def mock_measurement_set(self):
        """Mock measurement set object."""
        mock_measet = Mock()
        mock_measet.samples = ['/samples/IGVFSM001/', '/samples/IGVFSM002/']
        return mock_measet

    @pytest.fixture
    def mock_sample_with_subpool(self):
        """Mock sample object with cellular sub pool."""
        mock_sample = Mock()
        mock_sample.cellular_sub_pool = 'pool-1'
        return mock_sample

    @pytest.fixture
    def mock_sample_without_subpool(self):
        """Mock sample object without cellular sub pool."""
        mock_sample = Mock()
        mock_sample.cellular_sub_pool = None
        return mock_sample

    @pytest.fixture
    def mock_igvf_client(self, mock_measurement_set, mock_sample_with_subpool, mock_sample_without_subpool):
        """Mock IGVF client API."""
        mock_client = Mock()

        # Mock responses for different object types
        mock_responses = {
            '/measurement-sets/IGVFDS001/': Mock(actual_instance=mock_measurement_set),
            '/samples/IGVFSM001/': Mock(actual_instance=mock_sample_with_subpool),
            '/samples/IGVFSM002/': Mock(actual_instance=mock_sample_without_subpool)
        }

        mock_client.get_by_id.side_effect = lambda obj_id: mock_responses.get(
            obj_id)
        return mock_client

    def test_init(self, mock_igvf_client):
        """Test GetSampleMetadata initialization."""
        input_file_sets = ['/measurement-sets/IGVFDS001/',
                           '/other-file-sets/IGVFFS001/']

        getter = utils.GetSampleMetadata(
            input_file_sets=input_file_sets,
            igvf_client_api=mock_igvf_client
        )

        assert getter.input_file_sets == input_file_sets
        assert getter.igvf_client_api == mock_igvf_client

    def test_get_sample_ids(self, mock_igvf_client):
        """Test _get_sample_ids method."""
        input_file_sets = ['/measurement-sets/IGVFDS001/',
                           '/other-file-sets/IGVFFS001/']

        getter = utils.GetSampleMetadata(
            input_file_sets=input_file_sets,
            igvf_client_api=mock_igvf_client
        )

        sample_ids = getter._get_sample_ids()

        # Should only process measurement sets, not other file sets
        expected_sample_ids = ['/samples/IGVFSM001/', '/samples/IGVFSM002/']
        assert sample_ids == expected_sample_ids

    def test_get_sample_subpool_ids_with_subpools(self, mock_igvf_client):
        """Test _get_sample_subpool_ids method when subpools exist."""
        input_file_sets = ['/measurement-sets/IGVFDS001/']

        getter = utils.GetSampleMetadata(
            input_file_sets=input_file_sets,
            igvf_client_api=mock_igvf_client
        )

        # Mock the _get_sample_ids method call
        with patch.object(getter, '_get_sample_ids', return_value=['/samples/IGVFSM001/', '/samples/IGVFSM002/']):
            subpool_ids = getter._get_sample_subpool_ids()

        # Should return sorted list of unique subpool IDs
        expected_subpool_ids = ['pool-1']  # Only one sample has a subpool
        assert subpool_ids == expected_subpool_ids

    def test_get_sample_subpool_ids_no_subpools(self, mock_igvf_client):
        """Test _get_sample_subpool_ids method when no subpools exist."""
        input_file_sets = ['/measurement-sets/IGVFDS001/']

        # Mock samples without subpools
        mock_sample_no_subpool = Mock()
        mock_sample_no_subpool.cellular_sub_pool = None

        mock_client = Mock()
        mock_client.get_by_id.return_value = Mock(
            actual_instance=mock_sample_no_subpool)

        getter = utils.GetSampleMetadata(
            input_file_sets=input_file_sets,
            igvf_client_api=mock_client
        )

        # Mock the _get_sample_ids method call
        with patch.object(getter, '_get_sample_ids', return_value=['/samples/IGVFSM001/']):
            subpool_ids = getter._get_sample_subpool_ids()

        # Should return None when no subpools exist
        assert subpool_ids is None

    def test_get_samples_metadata(self, mock_igvf_client):
        """Test get_samples_metadata method returns SampleMetadata dataclass."""
        input_file_sets = ['/measurement-sets/IGVFDS001/']

        getter = utils.GetSampleMetadata(
            input_file_sets=input_file_sets,
            igvf_client_api=mock_igvf_client
        )

        # Don't mock the internal methods, let them run with proper fixtures
        metadata = getter.get_samples_metadata()

        # Should return SampleMetadata instance with accession strings (extracted from paths)
        assert isinstance(metadata, utils.SampleMetadata)
        assert metadata.sample_accessions == ['IGVFSM001', 'IGVFSM002']
        assert metadata.subpool_ids == ['pool-1']


class TestSampleMetadata:
    """Test the SampleMetadata dataclass."""

    def test_generate_alias_body_with_both(self):
        """Test generate_alias_body with both sample accessions and subpool IDs."""
        metadata = utils.SampleMetadata(
            sample_accessions=['IGVFSM001', 'IGVFSM002'],
            subpool_ids=['pool-1', 'pool-2']
        )

        # Updated to match actual implementation with AnalysisSet_ prefix and single-cell-uniform-pipeline suffix
        expected = 'AnalysisSet_IGVFSM001-IGVFSM002_pool-1-pool-2_single-cell-uniform-pipeline'
        assert metadata.generate_alias_body() == expected

    def test_generate_alias_body_no_subpools(self):
        """Test generate_alias_body with no subpool IDs."""
        metadata = utils.SampleMetadata(
            sample_accessions=['IGVFSM001', 'IGVFSM002'],
            subpool_ids=None
        )

        expected = 'AnalysisSet_IGVFSM001-IGVFSM002_single-cell-uniform-pipeline'
        assert metadata.generate_alias_body() == expected

    def test_generate_alias_body_no_samples(self):
        """Test generate_alias_body with no sample accessions."""
        metadata = utils.SampleMetadata(
            sample_accessions=None,
            subpool_ids=['pool-1', 'pool-2']
        )

        expected = 'AnalysisSet_pool-1-pool-2_single-cell-uniform-pipeline'
        assert metadata.generate_alias_body() == expected

    def test_generate_alias_body_empty_lists(self):
        """Test generate_alias_body with empty lists."""
        metadata = utils.SampleMetadata(
            sample_accessions=[],
            subpool_ids=[]
        )

        expected = 'AnalysisSet_single-cell-uniform-pipeline'
        assert metadata.generate_alias_body() == expected

    def test_frozen_dataclass(self):
        """Test that SampleMetadata is frozen (immutable)."""
        metadata = utils.SampleMetadata(
            sample_accessions=['IGVFSM001'],
            subpool_ids=['pool-1']
        )

        # Should raise an error when trying to modify
        with pytest.raises(dataclasses.FrozenInstanceError):
            metadata.sample_accessions = ['IGVFSM002']


class TestCreateAnalysisSetPayload(TestAnalysisSetUtils):

    def test_create_analysis_set_payload(self):
        """Test creating analysis set payload."""
        input_file_sets = ['/measurement-sets/IGVFDS123456/',
                           '/measurement-sets/IGVFDS789012/']

        # Mock GetSampleMetadata class and its methods
        mock_sample_metadata = Mock()
        # Fixed: Use IGVFSM (sample) accessions instead of IGVFDS (measurement set) accessions
        mock_sample_metadata.generate_alias_body.return_value = 'AnalysisSet_IGVFSM111111-IGVFSM222222_single-cell-uniform-pipeline'

        with patch.object(utils, 'GetSampleMetadata') as mock_get_metadata_class:
            mock_getter = Mock()
            mock_getter.get_samples_metadata.return_value = mock_sample_metadata
            mock_get_metadata_class.return_value = mock_getter

            result = utils.create_analysis_set_payload(
                input_file_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        expected = {
            '_profile': 'analysis_set',
            'input_file_sets': input_file_sets,
            'lab': self.test_lab,
            'award': self.test_award,
            'aliases': ['j-michael-cherry:AnalysisSet_IGVFSM111111-IGVFSM222222_single-cell-uniform-pipeline'],
            'file_set_type': 'intermediate analysis',
            'uniform_pipeline_status': 'preprocessing'
        }

        self.assertEqual(result, expected)


class TestIntegrationScenarios(TestAnalysisSetUtils):
    """Integration tests combining multiple functions."""

    def test_full_workflow_scenario(self):
        """Test a complete workflow from query to payload creation."""
        # Mock measurement set query result
        mock_measet = Mock()
        mock_measet.id = '/measurement-sets/IGVFDS123456/'
        mock_measet.related_multiome_datasets = [
            '/measurement-sets/IGVFDS789012/']
        mock_measet.barcode_replacement_file = None

        # Mock that it's not duplicated
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            # Test calc_input_file_sets_single
            input_sets = utils.calc_input_file_sets_single(
                mock_measet, self.mock_igvf_client)

        expected_input_sets = [
            '/measurement-sets/IGVFDS123456/', '/measurement-sets/IGVFDS789012/']
        self.assertEqual(sorted(input_sets), sorted(expected_input_sets))

        # Mock sample metadata
        mock_sample_metadata = Mock()
        mock_sample_metadata.generate_alias_body.return_value = 'AnalysisSet_SAMPLE1-SAMPLE2_single-cell-uniform-pipeline'

        with patch.object(utils, 'GetSampleMetadata') as mock_get_metadata_class:
            mock_getter = Mock()
            mock_getter.get_samples_metadata.return_value = mock_sample_metadata
            mock_get_metadata_class.return_value = mock_getter

            payload = utils.create_analysis_set_payload(
                input_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        # Verify payload structure
        self.assertEqual(payload['input_file_sets'], input_sets)
        self.assertEqual(payload['lab'], self.test_lab)
        self.assertEqual(payload['award'], self.test_award)
        self.assertIn(
            'j-michael-cherry:AnalysisSet_SAMPLE1-SAMPLE2_single-cell-uniform-pipeline', payload['aliases'])

    def test_full_workflow_scenario_with_barcode_replacement(self):
        """Test a complete workflow with barcode replacement file from query to payload creation."""
        # Mock measurement set query result with barcode replacement file
        mock_measet = Mock()
        mock_measet.id = '/measurement-sets/IGVFDS123456/'
        mock_measet.related_multiome_datasets = [
            '/measurement-sets/IGVFDS789012/']
        mock_measet.barcode_replacement_file = '/tabular-files/IGVFFI111111/'

        # Mock barcode replacement file response
        mock_brf_response = Mock()
        mock_brf_obj = Mock()
        mock_brf_obj.file_set = '/curated-sets/IGVFDS333333/'
        mock_brf_response.actual_instance = mock_brf_obj
        self.mock_igvf_client.get_by_id.return_value = mock_brf_response

        # Mock that it's not duplicated
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            # Test calc_input_file_sets_single
            input_sets = utils.calc_input_file_sets_single(
                mock_measet, self.mock_igvf_client)

        expected_input_sets = [
            '/measurement-sets/IGVFDS123456/',
            '/curated-sets/IGVFDS333333/',
            '/measurement-sets/IGVFDS789012/'
        ]
        self.assertEqual(sorted(input_sets), sorted(expected_input_sets))

        # Verify the API was called to get barcode replacement file info
        self.mock_igvf_client.get_by_id.assert_called_once_with(
            '/tabular-files/IGVFFI111111/')

        # Mock sample metadata
        mock_sample_metadata = Mock()
        mock_sample_metadata.generate_alias_body.return_value = 'AnalysisSet_SAMPLE1-SAMPLE2-SAMPLE3_single-cell-uniform-pipeline'

        with patch.object(utils, 'GetSampleMetadata') as mock_get_metadata_class:
            mock_getter = Mock()
            mock_getter.get_samples_metadata.return_value = mock_sample_metadata
            mock_get_metadata_class.return_value = mock_getter

            payload = utils.create_analysis_set_payload(
                input_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        # Verify payload structure includes all three file sets
        self.assertEqual(payload['input_file_sets'], input_sets)
        self.assertEqual(payload['lab'], self.test_lab)
        self.assertEqual(payload['award'], self.test_award)
        self.assertIn(
            'j-michael-cherry:AnalysisSet_SAMPLE1-SAMPLE2-SAMPLE3_single-cell-uniform-pipeline', payload['aliases'])

        # Verify it includes the curated set from barcode replacement
        self.assertIn('/curated-sets/IGVFDS333333/',
                      payload['input_file_sets'])
        self.assertEqual(len(payload['input_file_sets']), 3)

    def test_full_workflow_scenario_with_subpools(self):
        """Test a complete workflow with subpools."""
        # Mock measurement set query result
        mock_measet = Mock()
        mock_measet.id = '/measurement-sets/IGVFDS123456/'
        mock_measet.related_multiome_datasets = [
            '/measurement-sets/IGVFDS789012/']
        mock_measet.barcode_replacement_file = None

        # Mock that it's not duplicated
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            # Test calc_input_file_sets_single
            input_sets = utils.calc_input_file_sets_single(
                mock_measet, self.mock_igvf_client)

        expected_input_sets = [
            '/measurement-sets/IGVFDS123456/', '/measurement-sets/IGVFDS789012/']
        self.assertEqual(sorted(input_sets), sorted(expected_input_sets))

        # Mock sample metadata with subpools
        mock_sample_metadata = Mock()
        mock_sample_metadata.generate_alias_body.return_value = 'AnalysisSet_IGVFSM001-IGVFSM002_pool-1-pool-2_single-cell-uniform-pipeline'

        with patch.object(utils, 'GetSampleMetadata') as mock_get_metadata_class:
            mock_getter = Mock()
            mock_getter.get_samples_metadata.return_value = mock_sample_metadata
            mock_get_metadata_class.return_value = mock_getter

            payload = utils.create_analysis_set_payload(
                input_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        # Verify payload structure
        self.assertEqual(payload['input_file_sets'], input_sets)
        self.assertEqual(payload['lab'], self.test_lab)
        self.assertEqual(payload['award'], self.test_award)
        self.assertIn(
            'j-michael-cherry:AnalysisSet_IGVFSM001-IGVFSM002_pool-1-pool-2_single-cell-uniform-pipeline', payload['aliases'])

    def test_full_workflow_scenario_without_subpools(self):
        """Test a complete workflow without subpools."""
        # Mock measurement set query result
        mock_measet = Mock()
        mock_measet.id = '/measurement-sets/IGVFDS123456/'
        mock_measet.related_multiome_datasets = [
            '/measurement-sets/IGVFDS789012/']
        mock_measet.barcode_replacement_file = None

        # Mock that it's not duplicated
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            # Test calc_input_file_sets_single
            input_sets = utils.calc_input_file_sets_single(
                mock_measet, self.mock_igvf_client)

        expected_input_sets = [
            '/measurement-sets/IGVFDS123456/', '/measurement-sets/IGVFDS789012/']
        self.assertEqual(sorted(input_sets), sorted(expected_input_sets))

        # Mock sample metadata without subpools
        mock_sample_metadata = Mock()
        mock_sample_metadata.generate_alias_body.return_value = 'AnalysisSet_IGVFSM001-IGVFSM002_single-cell-uniform-pipeline'

        with patch.object(utils, 'GetSampleMetadata') as mock_get_metadata_class:
            mock_getter = Mock()
            mock_getter.get_samples_metadata.return_value = mock_sample_metadata
            mock_get_metadata_class.return_value = mock_getter

            payload = utils.create_analysis_set_payload(
                input_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        # Verify payload structure
        self.assertEqual(payload['input_file_sets'], input_sets)
        self.assertEqual(payload['lab'], self.test_lab)
        self.assertEqual(payload['award'], self.test_award)
        self.assertIn(
            'j-michael-cherry:AnalysisSet_IGVFSM001-IGVFSM002_single-cell-uniform-pipeline', payload['aliases'])


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)
