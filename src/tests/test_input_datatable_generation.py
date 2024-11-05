import sys
import os

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from sc_pipe_management.portal_to_terra_input_from_anaset import (
    generate_pipeline_input_table,
    generate_ordered_read_ids,
    seqspec_index_get,
    generate_finalinclusion_list,
    get_today_mmddyyyy
)
from sc_pipe_management.igvf_and_terra_api_tools import set_up_api_keys, get_igvf_client_auth
import unittest
from pandas.testing import assert_series_equal, assert_frame_equal
import pandas as pd
import os
import hashlib


# Seqspec files for different assays (Splitseq is not released, needs to download when running the test)
SEQSPEC_FILES_BY_ASSAY_TITLES = {
    '10x multiome': {'atac': 'src/tests/test_files/IGVFFI1365BOUI.yaml.gz',
                     'rna': 'src/tests/test_files/IGVFFI2433YUEJ.yaml.gz'},  # 10x multiome RNA seqspec
    'parse splitseq': {'atac': None,
                       'rna': 'src/tests/test_files/IGVFFI2264BQQD.yaml.gz'},  # SPLiT-seq RNA seqspec
    'shareseq': {'atac': 'src/tests/test_files/IGVFFI8012OCZQ.yaml.gz',  # ShareSeq does not have atac seqspec
                 'rna': 'src/tests/test_files/IGVFFI5825ATCM.yaml.gz'}  # ShareSeq RNA seqspec
}

# Expected outputs from seqspec
SEQSPEC_OUTPUT_BY_ASSAY_TITLES = {
    '10x multiome': {
        'atac': {'read_index': 'bc:8:23:-,r1:0:49,r2:0:49',
                 'index_input': 'IGVFFI4665EVGC,IGVFFI4986VMDU,IGVFFI0391NHGA',
                 'onlist_input': 'IGVFFI0391NHGA',
                 'onlist_final_list': '59ed4e3b7c64af921bde8b8497924c02',
                 'onlist_method': 'no combination'
                 },
        'rna': {'read_index': '0,0,16:0,16,28:1,0,90',
                'index_input': 'IGVFFI6948UZJO,IGVFFI7951DAQB',
                'onlist_input': 'IGVFFI7951DAQB',
                'onlist_final_list': '5b223bbf2d45ba5e8bb55787872321c8',
                'onlist_method': 'no combination'
                }
    },
    'parse splitseq': {
        'atac': None,
        'rna': {'read_index': '1,10,18,1,48,56,1,78,86:1,0,10:0,0,140',
                'index_input': 'IGVFFI4633WXZI.fastq.gz,IGVFFI9320MADV.fastq.gz',
                'onlist_input': 'IGVFFI9320MADV.fastq.gz',
                'onlist_final_list': '7f1e51743561883977e5a905bb345ac8',
                'onlist_method': 'product'
                }
    },
    'shareseq': {
        'atac': {'read_index': 'bc:115:122,bc:153:160,bc:191:198,r1:0:99,r2:0:99',
                 'index_input': 'IGVFFI1918YZDJ.fastq.gz,IGVFFI4773YQEF.fastq.gz',
                 'onlist_input': 'IGVFFI4773YQEF.fastq.gz',
                 'onlist_final_list': 'f7eabd710222e594459698a6c624de7b',
                 'onlist_method': 'product'
                 },
        'rna': {'read_index': '1,115,123,1,153,161,1,191,199:1,0,10:0,0,100',
                'index_input': 'IGVFFI0777CMJH.fastq.gz,IGVFFI7330UNCB.fastq.gz',
                'onlist_input': 'IGVFFI7330UNCB.fastq.gz',
                'onlist_final_list': 'f7eabd710222e594459698a6c624de7b',
                'onlist_method': 'product'
                }
    }
}


def compute_md5sum(file_path: str) -> str:
    """Compute the MD5 checksum of a file.

    Args:
        file_path (str): File of interest to compute the md5sum for.

    Returns:
        str: MD5 checksum of the file in hexadecimal format.
    """
    md5_hash = hashlib.md5()
    with open(file_path, 'rb') as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b''):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()


def import_input_datatable_tsv(tsv_file_path: str) -> pd.DataFrame:
    """Import TSV tables

    Args:
        tsv_file_path (str): fixed paths to ref or generated TSV file.

    Returns:
        pd.DataFrame: final input table
    """
    return pd.read_csv(tsv_file_path, sep='\t', index_col=0)


def test_generate_pipeline_input_table() -> pd.DataFrame:
    """Generate the test pipeline input table for unittesting.

    Returns:
        pd.DataFrame: Test table to compare to the reference table.
    """
    # Generate the final input datatable
    query_analysis_set_accs = ['IGVFDS9664YLAD',    # 10X multiome
                               'IGVFDS0223KLTB',    # Parse SPLiT-seq
                               'IGVFDS0657NHTA'     # ShareSeq
                               ]
    pipeline_input_table = generate_pipeline_input_table(query_analysis_set_accs=query_analysis_set_accs,
                                                         igvf_api=IGVF_PROD_CLIENT_API,
                                                         terra_etype='unittest_pipeline_tester',
                                                         local_barcode_file_dir=os.path.join(
                                                             os.getcwd(), 'final_barcode_list/', get_today_mmddyyyy()),
                                                         gs_barcode_list_bucket='gs://unittest_mock_bucket/submissions/final_barcode_onlist/'
                                                         )
    pipeline_input_table.to_csv(TEST_TABLES_PATHS['generated'], sep='\t')
    return import_input_datatable_tsv(tsv_file_path=TEST_TABLES_PATHS['generated'])


# API tool
IGVF_ENDPOINT = 'production'  # or 'sandbox'
IGVF_PROD_API_KEYS = set_up_api_keys(igvf_endpoint=IGVF_ENDPOINT)
IGVF_PROD_CLIENT_API = get_igvf_client_auth(igvf_api_keys=IGVF_PROD_API_KEYS,
                                            igvf_site=IGVF_ENDPOINT)


# Relative paths to the test input datatable generation files
TEST_TABLES_PATHS = {'ref': 'src/tests/test_files/test_input_datatable_generation_reference.tsv',
                     'generated': 'src/tests/test_files/test_input_datatable_generation_generated.tsv'
                     }

# Reference tables
IMPORTED_INPUT_TABLES = {
    'ref': pd.read_csv(TEST_TABLES_PATHS['ref'], sep='\t', index_col=0),
    'generated': test_generate_pipeline_input_table()
}


class TestInputDatatableGeneration(unittest.TestCase):
    def setUp(self):
        """Get the values into the test class for comparison.
        """
        self.ref_table = IMPORTED_INPUT_TABLES['ref']
        self.generated_table = IMPORTED_INPUT_TABLES['generated']

    def test_generate_read_id_input(self):
        """Test if read id input -i arg is generated correctly.
        """
        for assay_title, seqspec_paths in SEQSPEC_FILES_BY_ASSAY_TITLES.items():
            for assay_type, seqspec_path in seqspec_paths.items():
                if seqspec_path is None:
                    # Skip if the seqspec path is None (e.g., SPLiT-seq for atac)
                    print(
                        f'Skipping {assay_title} for {assay_type} as it has no seqspec.')
                    continue
                for read_id_methd in ['onlist', 'index']:
                    print(
                        f'Testing {assay_title} {assay_type} with read_id_methd: {read_id_methd}.')
                    print(seqspec_path)
                    curr_read_input = generate_ordered_read_ids(seqspec_file_path=seqspec_path,
                                                                assay_type=assay_type,
                                                                usage_purpose=read_id_methd,
                                                                igvf_api=IGVF_PROD_CLIENT_API
                                                                )
                    curr_error_detail = (f'Generated index read ids do not match expected for {assay_title} {assay_type} with method {read_id_methd}. '
                                         f'Expected: {SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type][f"{read_id_methd}_input"]}, '
                                         f' Got: {curr_read_input}.'
                                         )
                    self.assertEqual(
                        first=curr_read_input, second=SEQSPEC_OUTPUT_BY_ASSAY_TITLES[
                            assay_title][assay_type][f'{read_id_methd}_input'],
                        msg=curr_error_detail)
        print('All read id input tests done.\n')

    def test_read_index_generation(self):
        """Check if rna and atac read index are correct.
        """
        for assay_title, seqspec_paths in SEQSPEC_FILES_BY_ASSAY_TITLES.items():
            for assay_type, seqspec_path in seqspec_paths.items():
                if seqspec_path is None:
                    # Skip if the seqspec path is None (e.g., SPLiT-seq for atac)
                    print(
                        f'Skipping {assay_title} for {assay_type} as it has no seqspec.')
                    continue
                print(
                    f'Testing {assay_title} {assay_type} for read index generation.')
                # Get the current read index from seqspec\
                print(seqspec_path)
                curr_index_input = seqspec_index_get(
                    seqspec_file_path=seqspec_path, assay_type=assay_type, igvf_api=IGVF_PROD_CLIENT_API)
                curr_error_detail = (f'Generated {assay_type} index read index does not match expected for {assay_title} {assay_type}. '
                                     f'Expected: {SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type]["read_index"]}, '
                                     f' Got: {curr_index_input}.'
                                     )
                self.assertEqual(
                    first=curr_index_input,
                    second=SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type]['read_index'],
                    msg=curr_error_detail
                )
        print('All read index tests done.\n')

    def test_finalinclusion_list_generation(self):
        """Check if the final inclusion list file has correct md5sum for each assay type.
        """
        for assay_title, seqspec_paths in SEQSPEC_FILES_BY_ASSAY_TITLES.items():
            for assay_type, seqspec_path in seqspec_paths.items():
                if seqspec_path is None:
                    # Skip if the seqspec path is None (e.g., SPLiT-seq for atac)
                    print(
                        f'Skipping {assay_title} for {assay_type} as it has no seqspec.')
                    continue
                print(
                    f'Testing final inclusion list generation for {assay_title} {assay_type}.')
                print(seqspec_path)
                output_file_path = generate_finalinclusion_list(
                    seqspec_file_path=seqspec_path,
                    assay_type=assay_type,
                    onlist_method=SEQSPEC_OUTPUT_BY_ASSAY_TITLES[assay_title][assay_type]['onlist_method'],
                    final_inclusion_list_path=os.path.join(
                        'src/tests/test_files/', f'{assay_title}_{assay_type}_final_barcode_list.txt'),
                    igvf_api=IGVF_PROD_CLIENT_API
                )
                curr_md5sum = compute_md5sum(
                    file_path=output_file_path
                )
                # Compare the md5sum of the generated file with the expected md5sum
                expected_md5sum = SEQSPEC_OUTPUT_BY_ASSAY_TITLES[
                    assay_title][assay_type]['onlist_final_list']
                curr_error_detail = (f'Generated final inclusion list md5sum does not match expected for {assay_title} {assay_type}. '
                                     f'Expected: {expected_md5sum}, '
                                     f' Got: {curr_md5sum}.'
                                     )
                self.assertEqual(
                    first=curr_md5sum,
                    second=expected_md5sum,
                    msg=curr_error_detail
                )
        print('All final inclusion list tests done.\n')

    def test_file_md5sums_comparison(self):
        """Compare the generated input file md5sums to ensure they match the reference file.
        """
        print('Checking table file md5sums...')
        self.assertEqual(first=compute_md5sum(TEST_TABLES_PATHS['ref']),
                         second=compute_md5sum(TEST_TABLES_PATHS['ref']),
                         msg='MD5 checksum of the generated file does not match the reference file.'
                         )
        print('MD5sums check done.\n')

    def test_df_shapes_comparison(self):
        """Check if the shapes of the generated input datatable match the reference input datatable
        """
        print('Checking datatable shapes...')
        try:
            self.assertEqual(first=self.ref_table.shape,
                             second=self.generated_table.shape, msg='Shapes do not match.')
        except AssertionError as e:
            print(e)
        print('Shapes comparison done.\n')

    def test_column_comparison(self):
        """Check if the columns of the generated input datatable match the reference input datatable
        """
        print('Checking datatable columns...')
        try:
            self.assertEqual(first=sorted(self.ref_table.columns), second=sorted(
                self.generated_table.columns), msg='Columns do not match.')
        except AssertionError as e:
            print(e)
        print('Columns comparison done.\n')

    def test_row_comparison(self):
        """Check if the rows of the generated input datatable match the reference input datatable
        """
        print('Checking datatable rows...')
        try:
            self.assertEqual(first=sorted(self.ref_table.index), second=sorted(
                self.generated_table.index), msg='Rows do not match.')
        except AssertionError as e:
            print(e)
        print('Rows comparison done.\n')

    def test_column_value_comparison(self):
        # Check if the values of the generated input datatable match the reference input datatable
        print('Checking datatable values...')
        try:
            for col in self.generated_table.columns:
                assert_series_equal(
                    left=self.ref_table[col], right=self.generated_table[col])
        except AssertionError as e:
            print(e)
        print('Column-wise value comparison done.\n')

    def test_full_table_comparison(self):
        """Check if the generated input datatable is equal to the reference input datatable
        """
        try:
            print('Checking full datatables...')
            assert_frame_equal(left=self.ref_table,
                               right=self.generated_table,
                               check_dtype=True,
                               check_index_type='equiv',
                               check_column_type=True
                               )
        except AssertionError as e:
            print(e)
        print('Total full datatable comparison done.\n')


def suite():
    """Set up a test suite for all the tests in this class in a specific order.
    """
    test_suite = unittest.TestSuite()
    test_suite.addTest(TestInputDatatableGeneration(
        'test_generate_read_id_input'))
    test_suite.addTest(TestInputDatatableGeneration(
        'test_read_index_generation'))
    test_suite.addTest(TestInputDatatableGeneration(
        'test_finalinclusion_list_generation'))
    test_suite.addTest(TestInputDatatableGeneration(
        'test_file_md5sums_comparison'))
    test_suite.addTest(TestInputDatatableGeneration(
        'test_df_shapes_comparison'))
    test_suite.addTest(TestInputDatatableGeneration('test_column_comparison'))
    test_suite.addTest(TestInputDatatableGeneration('test_row_comparison'))
    test_suite.addTest(TestInputDatatableGeneration(
        'test_column_value_comparison'))
    test_suite.addTest(TestInputDatatableGeneration(
        'test_full_table_comparison'))
    return test_suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
