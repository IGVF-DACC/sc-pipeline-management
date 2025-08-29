import sys
import os

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import unittest
from unittest.mock import Mock, patch, MagicMock
import sc_pipe_management.wranger_utils.analysis_set_setup_utils as utils

# TODO: this will be replaced by AnalysisSet pipeline status once the ticket is in


class TestAnalysisSetUtils(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.mock_igvf_client = Mock()
        self.mock_igvf_utils = Mock()

        # Sample test data
        self.test_input_file_sets = [
            '/measurement-sets/IGVFMS123456/', '/measurement-sets/IGVFMS789012/']
        self.test_lab = '/labs/j-michael-cherry/'
        self.test_award = '/awards/HG012012/'


class TestCheckIsScpipeline(TestAnalysisSetUtils):
    def test_check_is_scpipeline_not_analysis_set(self):
        """Test that non-analysis-set ID returns None."""
        # Assertion check: not single cell
        result = utils.check_is_scpipeline(
            '/measurement-sets/IGVFMS123456/', self.mock_igvf_client)
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
        mock_query_res.id = '/measurement-sets/IGVFMS123456/'
        mock_query_res.related_multiome_datasets = [
            '/measurement-sets/IGVFMS789012/']
        mock_query_res.barcode_replacement_file = None

        # Mock check_is_duped_for_all to return False
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        expected = ['/measurement-sets/IGVFMS123456/',
                    '/measurement-sets/IGVFMS789012/']
        self.assertEqual(sorted(result), sorted(expected))

    def test_calc_input_file_sets_single_parse(self):
        """Test basic input file set calculation of Parse SPLiT-seq."""
        # Mock query result
        mock_query_res = Mock()
        mock_query_res.id = '/measurement-sets/IGVFMS123456/'
        mock_query_res.related_multiome_datasets = [
            '/measurement-sets/IGVFMS789012/']
        mock_query_res.barcode_replacement_file = '/tabular-files/IGVFFI111111/'

        # Mock barcode replacement file response
        mock_brf_response = Mock()
        mock_brf_obj = Mock()
        mock_brf_obj.file_set = '/curated-sets/IGVFMS333333/'
        mock_brf_response.actual_instance = mock_brf_obj
        self.mock_igvf_client.get_by_id.return_value = mock_brf_response

        # Mock check_is_duped_for_all to return False
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        expected = ['/measurement-sets/IGVFMS123456/',
                    '/curated-sets/IGVFMS333333/', '/measurement-sets/IGVFMS789012/']
        self.assertEqual(sorted(result), sorted(expected))

    def test_calc_input_file_sets_single_no_related_datasets(self):
        """Test input file set calculation without related datasets."""
        mock_query_res = Mock()
        mock_query_res.id = '/measurement-sets/IGVFMS123456/'
        mock_query_res.related_multiome_datasets = None
        mock_query_res.barcode_replacement_file = None

        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        expected = ['/measurement-sets/IGVFMS123456/']
        self.assertEqual(result, expected)

    def test_calc_input_file_sets_single_is_duplicate(self):
        """Test that duplicate analysis set returns None."""
        mock_query_res = Mock()
        mock_query_res.id = '/measurement-sets/IGVFMS123456/'

        with patch.object(utils, 'check_is_duped_for_all', return_value=True):
            result = utils.calc_input_file_sets_single(
                mock_query_res, self.mock_igvf_client)

        self.assertIsNone(result)


class TestGetSampleAccessions(TestAnalysisSetUtils):

    def test_get_sample_accessions(self):
        """Test getting sample accessions from input file sets."""
        input_file_sets = ['/measurement-sets/IGVFMS123456/',
                           '/measurement-sets/IGVFMS789012/']

        # Mock the first measurement set
        mock_file_set_1 = Mock()
        mock_file_set_1.samples = ['/samples/IGVFSM111111/']

        # Mock the second measurement set
        mock_file_set_2 = Mock()
        mock_file_set_2.samples = ['/samples/IGVFSM333333/']  # One duplicate

        self.mock_igvf_client.get_by_id.side_effect = [
            mock_file_set_1, mock_file_set_2]

        result = utils.get_sample_accessions(
            input_file_sets, self.mock_igvf_client)

        # Should return sorted unique sample accessions joined by underscore
        expected_samples = ['IGVFSM111111', 'IGVFSM333333']
        expected_result = '_'.join(expected_samples)
        self.assertEqual(result, expected_result)

    def test_get_sample_accessions_non_measurement_sets(self):
        """Test that non-measurement-set IDs are skipped."""
        input_file_sets = ['/analysis-sets/IGVFAS123456/',
                           '/measurement-sets/IGVFMS789012/']

        # Mock only the measurement set
        mock_file_set = Mock()
        mock_file_set.samples = ['/samples/IGVFSM111111/']
        self.mock_igvf_client.get_by_id.return_value = mock_file_set

        result = utils.get_sample_accessions(
            input_file_sets, self.mock_igvf_client)

        # Should only process the measurement set
        self.assertEqual(result, 'IGVFSM111111')
        self.mock_igvf_client.get_by_id.assert_called_once_with(
            '/measurement-sets/IGVFMS789012/')


class TestCreateAnalysisSetPayload(TestAnalysisSetUtils):

    def test_create_analysis_set_payload(self):
        """Test creating analysis set payload."""
        input_file_sets = ['/measurement-sets/IGVFMS123456/',
                           '/measurement-sets/IGVFMS789012/']

        with patch.object(utils, 'get_sample_accessions', return_value='IGVFSM111111_IGVFSM222222'):
            result = utils.create_analysis_set_payload(
                input_file_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        expected = {
            '_profile': 'analysis_set',
            'input_file_sets': input_file_sets,
            'lab': self.test_lab,
            'award': self.test_award,
            'aliases': ['j-michael-cherry:Sample-IGVFSM111111_IGVFSM222222_single-cell-uniform-pipeline'],
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
        mock_measet.id = '/measurement-sets/IGVFMS123456/'
        mock_measet.related_multiome_datasets = [
            '/measurement-sets/IGVFMS789012/']
        mock_measet.barcode_replacement_file = None

        # Mock that it's not duplicated
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            # Test calc_input_file_sets_single
            input_sets = utils.calc_input_file_sets_single(
                mock_measet, self.mock_igvf_client)

        expected_input_sets = [
            '/measurement-sets/IGVFMS123456/', '/measurement-sets/IGVFMS789012/']
        self.assertEqual(sorted(input_sets), sorted(expected_input_sets))

        # Mock sample accessions
        with patch.object(utils, 'get_sample_accessions', return_value='SAMPLE1_SAMPLE2'):
            payload = utils.create_analysis_set_payload(
                input_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        # Verify payload structure
        self.assertEqual(payload['input_file_sets'], input_sets)
        self.assertEqual(payload['lab'], self.test_lab)
        self.assertEqual(payload['award'], self.test_award)
        self.assertIn(
            'j-michael-cherry:Sample-SAMPLE1_SAMPLE2_single-cell-uniform-pipeline', payload['aliases'])

    def test_full_workflow_scenario_with_barcode_replacement(self):
        """Test a complete workflow with barcode replacement file from query to payload creation."""
        # Mock measurement set query result with barcode replacement file
        mock_measet = Mock()
        mock_measet.id = '/measurement-sets/IGVFMS123456/'
        mock_measet.related_multiome_datasets = [
            '/measurement-sets/IGVFMS789012/']
        mock_measet.barcode_replacement_file = '/tabular-files/IGVFFI111111/'

        # Mock barcode replacement file response
        mock_brf_response = Mock()
        mock_brf_obj = Mock()
        mock_brf_obj.file_set = '/curated-sets/IGVFMS333333/'
        mock_brf_response.actual_instance = mock_brf_obj
        self.mock_igvf_client.get_by_id.return_value = mock_brf_response

        # Mock that it's not duplicated
        with patch.object(utils, 'check_is_duped_for_all', return_value=False):
            # Test calc_input_file_sets_single
            input_sets = utils.calc_input_file_sets_single(
                mock_measet, self.mock_igvf_client)

        expected_input_sets = [
            '/measurement-sets/IGVFMS123456/',
            '/curated-sets/IGVFMS333333/',
            '/measurement-sets/IGVFMS789012/'
        ]
        self.assertEqual(sorted(input_sets), sorted(expected_input_sets))

        # Verify the API was called to get barcode replacement file info
        self.mock_igvf_client.get_by_id.assert_called_once_with(
            '/tabular-files/IGVFFI111111/')

        # Mock sample accessions
        with patch.object(utils, 'get_sample_accessions', return_value='SAMPLE1_SAMPLE2_SAMPLE3'):
            payload = utils.create_analysis_set_payload(
                input_sets, self.test_lab, self.test_award, self.mock_igvf_client
            )

        # Verify payload structure includes all three file sets
        self.assertEqual(payload['input_file_sets'], input_sets)
        self.assertEqual(payload['lab'], self.test_lab)
        self.assertEqual(payload['award'], self.test_award)
        self.assertIn(
            'j-michael-cherry:Sample-SAMPLE1_SAMPLE2_SAMPLE3_single-cell-uniform-pipeline',
            payload['aliases'])

        # Verify it includes the curated set from barcode replacement
        self.assertIn('/curated-sets/IGVFMS333333/',
                      payload['input_file_sets'])
        self.assertEqual(len(payload['input_file_sets']), 3)


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)
