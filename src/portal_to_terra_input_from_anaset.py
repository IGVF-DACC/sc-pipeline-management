import urllib.parse
import os
import pandas as pd
import itertools
import urllib
from datetime import datetime
import subprocess
import requests


# TODO:
# 1) A SeqFile may have multiple seqspec files but only one is released. Once that is sorted out, need to
#   update the seqspec URLs to be the released one.
# 2) Make a func for generating full hrefs

# Read_names to read types
READ_NAME_TO_READ_TYPE_MAP = {'Read 1': 'read1',
                              'Read 2': 'read2',
                              'Barcode index': 'barcode'
                              }

# Matching assay term name to assay type
# NOTE: This list may need to be periodically updated.
ASSAY_NAMES_CONVERSION_REF = {
    '/assay-terms/OBI_0002764/': 'atac',    # single-cell ATAC-seq
    '/assay-terms/OBI_0002762/': 'atac',     # single-nucleus ATAC-seq
    '/assay-terms/OBI_0003109/': 'rna',     # single-nucleus RNA sequencing assay
    '/assay-terms/OBI_0002631/': 'rna',  # single-cell RNA sequencing assay
    # in vitro CRISPR screen using single-cell RNA-seq
    '/assay-terms/OBI_0003660/': 'rna'
}

# Preferred assay titles that will require ATAC-RNA onlist mapping
ASSAYS_NEED_ONLIST_MAPPING = [
    '10x multiome', '10x multiome with MULTI-seq']

# Taxa to genome ref files
TAXA_TO_GENOME_REF_FILES = {'Homo sapiens': ['GRCh38', '<human genome tsv>'],
                            'Mus musculus': ['GRCm39', '<mouse genome tsv>']
                            }

# For URL request
BASE_IGVF_PORTAL_URL = 'https://api.data.igvf.org'


def construct_full_href_url(igvf_href: str, base_url: str = BASE_IGVF_PORTAL_URL) -> str:
    """_summary_

    Args:
        igvf_href (str): _description_
        base_url (str, optional): _description_. Defaults to BASE_IGVF_PORTAL_URL.

    Returns:
        str: _description_
    """
    return urllib.parse.urljoin(base=base_url, url=igvf_href)


# Get sample level info (taxa, accession for subpool ID)
def get_sample_obj(sample_id: str, igvf_api):
    """Get the analysis set object from the IGVF API.

    Args:
        sample_id (str): e.g., /in-vitro-system/IGVFxxxxx
        igvf_api (_type_): _description_

    Returns:
        _type_: _description_
    """
    return igvf_api.get_by_id(sample_id).actual_instance


def get_sample_taxa(sample_object) -> str:
    """Get sample taxa from the sample object.

    Args:
        sample_object (_type_): Sample object

    Returns:
        str: Homo sapiens, Mus musculus, etc.
    """
    return sample_object.taxa


def get_sample_accession(sample_object) -> str:
    """Get sample accession.

    Args:
        sample_object (_type_): Sample object

    Returns:
        str: IGVFxxxx
    """
    return sample_object.accession


# Get analysis set level info
def get_analysis_set_obj(analysis_set_accessions: list, igvf_api):
    """Get the analysis set object from the IGVF API."""
    return igvf_api.analysis_sets(accession=analysis_set_accessions).graph[0]


def get_measurement_ids(analysis_set_obj) -> list:
    """Extract measurement set IDs from the analysis set object."""
    return [measet_id for measet_id in analysis_set_obj.input_file_sets if 'measurement-sets' in measet_id]


def get_and_sort_measurement_sets_by_assays(measurement_set_ids: list, igvf_api) -> list:
    """Sort a list of measurement set items based on assay terms.

    Args:
        measurement_set_ids (list): A list of measurement set items generated by the APi.
        igvf_api (_type_): _description_

    Returns:
        list: Same content as the input list but sorted by the assay terms.
    """
    return sorted([igvf_api.get_by_id(measet_id).actual_instance for measet_id in measurement_set_ids], key=lambda x: x.assay_term)


def get_seqfile_items(measet_items: list, igvf_api) -> list:
    """Get sequence file items and the preferred assay title.

    Args:
        measet_items (list): A list of measurement set items sorted by assay terms.
        igvf_api (_type_): _description_

    Returns:
        list: [sequence file items]
    """
    file_items = []
    # For each measurement set of the same assay type, look for seqfiles
    for measet_item in measet_items:
        for file in measet_item.files:
            if 'sequence-files' in file:
                file_items.append(igvf_api.get_by_id(file).actual_instance)
    return file_items


def seqfile_sort_func(x):
    """Set up sorting sorter for sequence file item lists.
    """
    return (x.file_set, x.illumina_read_type, x.sequencing_run, x.lane)


def get_seqfile_readnames(seqfile_item, assay_type: str) -> str:
    """HACK: Use sequence file metadata to infer FASTQ file read types.

    Args:
        seqfile_item (_type_): Return item from the API
        preferred_assay_title (str): _description_
        assay_type (str): translated from assay term into atac or rna

    Returns:
        str: e.g., atac_read1, rna_read2
    """
    read_file_titles = []
    curr_read_names = seqfile_item.read_names
    if not curr_read_names:
        raise Exception(
            'Error: No read names found in the sequence file item.')
    for curr_read_name in curr_read_names:
        read_file_titles.append(
            '_'.join([assay_type, READ_NAME_TO_READ_TYPE_MAP[curr_read_name]]))
    return read_file_titles


def get_seqspec_hrefs(seqspec_ids: list[str], igvf_api) -> set:
    # NOTE: See TO-DO #1
    """Get the seqspec URLs from the seqspec IDs.

    Args:
        seqspec_ids (list[str]): A list of seqspec IDs
        igvf_api (_type_): _description_

    Returns:
        set: _description_
    """
    seqspec_urls = set()
    for seqspec_id in seqspec_ids:
        curr_seqspec_item = igvf_api.get_by_id(seqspec_id).actual_instance
        seqspec_urls.add(construct_full_href_url(
            igvf_href=curr_seqspec_item.href))
    return seqspec_urls


def get_onlist_method(measurement_set_ids: list, igvf_api) -> str:
    """Get the onlist method of the measurement set items.

    Args:
        measurement_set_ids (list): _description_

    Returns:
        str: onlist method
    """
    curr_onlist_methods = set()
    for measet_id in measurement_set_ids:
        curr_measet_item = igvf_api.get_by_id(measet_id).actual_instance
        curr_onlist_methods.add(curr_measet_item.onlist_method)
    if len(curr_onlist_methods) != 1:
        raise Exception(
            'Error: Measurement sets of the same assay term have different onlist methods.')
    return list(curr_onlist_methods)[0]


def download_file_via_https(igvf_portal_href_url: str) -> str:
    """Download file via href url on the portal

    Args:
        igvf_portal_href_url (str): _description_
        save_file_path (str): _description_

    Returns:
        str: _description_
    """
    username = os.getenv('IGVF_API_KEY')
    password = os.getenv('IGVF_SECRET_KEY')
    session = requests.Session()
    session.auth = (username, password)
    response = session.get(igvf_portal_href_url)
    if response.status_code == 200:
        curr_output_file = os.path.join(
            './temp', igvf_portal_href_url.split('/')[-1])
        with open(curr_output_file, 'wb') as file:
            file.write(response.content)
        return curr_output_file
    else:
        raise Exception(
            f'Error: Download failed with status code: {response.status_code}.')


def check_emptyinclusion_list(final_inclusion_list_path: str) -> bool:
    """Check if the final inclusion list created is empty.

    Args:
        final_inclusion_list_path (str): _description_

    Returns:
        bool: _description_
    """
    with open(final_inclusion_list_path, 'r') as file_obj:
        first_char = file_obj.read(1)
        if not first_char.isalpha():
            return False


# TODO: This always fails now because the final_inclusion_list_path always returns ./
def generate_finalinclusion_list(seqspec_file_path: str, assay_type: str, onlist_method: str, final_inclusion_list_path: str) -> str:
    """Generate the final barcode inclusion list txt file.

    Args:
        seqspec_file_path (str): _description_
        assay_type (str): _description_
        onlist_method (str): no combination, product, multi
        final_inclusion_list_path (str): _description_

    Returns:
        str: _description_
    """
    # If no combination, seqspec onlist just repeats the name/url of the only onlist file listed. It needs to be downloaded and renamed.
    if onlist_method == 'no combination':
        curr_run_log = subprocess.run(['seqspec', 'onlist', '-m', assay_type, '-s',
                                      'region-type', '-i', 'barcode', '-o', final_inclusion_list_path, seqspec_file_path], capture_output=True)
    # If combinatorial, seqspec onlist will generated a new file
    else:
        curr_run_log = subprocess.run(['seqspec', 'onlist', '-m', assay_type, '-s', 'region-type', '-i',
                                      'barcode', '-f', onlist_method, '-o', final_inclusion_list_path, seqspec_file_path], capture_output=True)
    if curr_run_log.returncode != 0:
        raise Exception('Error: Seqspec onlist generation command error.')
    if check_emptyinclusion_list(final_inclusion_list_path=final_inclusion_list_path) is False:
        raise Exception('Error: Empty final inclusion list file generated.')
    return final_inclusion_list_path


def seqspec_onlist_safetychk_and_get(seqspec_file_paths: list, assay_type: str, onlist_method: str, final_inclusion_list_path: str) -> str:
    """Do a pre-flight safety check for seqspec onlist to ensure all seqspec files for the same pipeline have the same onlist info.

    Args:
        seqspec_file_paths (list): a list of locally saved seqspec file paths
        assay_type (str): rna or atac

    Returns:
        str: Either error messages or the path to the first seqspec file in the list.
    """
    onlist_chk_lst = []
    for file in seqspec_file_paths:
        curr_run_log = subprocess.run(['seqspec', 'file', '-m', assay_type, '-s', 'region-type',
                                      '-i', 'barcode', '-f', 'index', '-k', 'url', file], capture_output=True)
        if curr_run_log.returncode == 0:
            curr_onlist_files = curr_run_log.stdout.decode(
                'utf-8').strip().split(',')
            onlist_chk_lst.append(curr_onlist_files)
        else:
            raise Exception('Error: seqspec file command failed.')
    # If the onlist files are different, return an error message
    if not all(set(sublist) == set(onlist_chk_lst[0]) for sublist in onlist_chk_lst):
        raise Exception(
            'Error: Sequence files have different onlist files in associated seqspecs.')
    else:
        # If there are multiple onlist files, return error if the onlist method is no combination.
        if (len(onlist_chk_lst[0]) > 1) and (onlist_method == 'no combination'):
            raise Exception(
                'Error: Multiple onlist files found but onlist method is no combination.')
        # If there is only one onlist file, return error if the onlist method isn't no combination.
        elif (len(onlist_chk_lst[0]) == 1) and (onlist_method != 'no combination'):
            raise Exception(
                'Error: Only one onlist file found but onlist method is not no combination.')
        # Otherwise, make the final barcode list
        else:
            return generate_finalinclusion_list(seqspec_file_path=str(
                seqspec_file_paths[0]), assay_type=assay_type, onlist_method=onlist_method, final_inclusion_list_path=final_inclusion_list_path)


def seqspec_index_safetychk_and_get(seqspec_file_paths: list, assay_type: str) -> str:
    """Pre-flight safety check for read_format index and generate it.

    Args:
        seqspec_file_paths (list): _description_
        assay_type (str): _description_

    Returns:
        str: _description_
    """
    assay_type_to_tool_format = {'rna': 'kb', 'atac': 'chromap'}
    readformat_chk_lst = []
    for file in seqspec_file_paths:
        curr_run_log = subprocess.run(['seqspec', 'index', '-m', assay_type, '-t',
                                      assay_type_to_tool_format[assay_type], '-s', 'file', file], capture_output=True)
        if curr_run_log.returncode != 0:
            raise Exception('Error: seqspec index command failed.')
        if assay_type == 'rna':
            readformat_chk_lst.append(
                curr_run_log.stdout.decode('utf-8').strip())
        elif assay_type == 'atac':
            readformat_chk_lst.append(
                curr_run_log.stdout.decode('utf-8').strip().split(' ')[-1])
    if len(set(readformat_chk_lst)) != 1:
        raise Exception(
            'Error: Sequence files have different read_format in associated seqspecs.')
    return readformat_chk_lst[0]


def get_onlist_mapping_status(measurement_set_ids: list, igvf_api) -> bool:
    """Get the onlist mapping status of the measurement set items.

    Args:
        measurement_set_ids (list): _description_

    Returns:
        bool: True or False
    """
    curr_onlist_mapping = []
    for measet_id in measurement_set_ids:
        curr_measet_item = igvf_api.get_by_id(measet_id).actual_instance
        if curr_measet_item.preferred_assay_title in ASSAYS_NEED_ONLIST_MAPPING:
            curr_onlist_mapping.append(True)
        else:
            curr_onlist_mapping.append(False)
    if len(set(curr_onlist_mapping)) > 1:
        return 'Error: Measurement sets have different onlist mapping status.'
    else:
        return curr_onlist_mapping[0]


class SingleCellInputBuilder:
    def __init__(self, analysis_set_acc: str, igvf_api):
        self.data = {'analysis_set_acc': '',
                     'atac_MeaSetIDs': [],
                     'rna_MeaSetIDs': [],
                     'subpool_id': '',
                     'taxa': '',
                     'genome_assembly': '',
                     'genome_ref': '',
                     'atac_read1_accessions': [],
                     'atac_read2_accessions': [],
                     'atac_barcode_accessions': [],
                     'rna_read1_accessions': [],
                     'rna_read2_accessions': [],
                     'atac_seqspec_urls': set(),  # NOTE: The seqspec here is temporary for safety check
                     'rna_seqspec_urls': set(),
                     'atac_read1': [],
                     'atac_read2': [],
                     'atac_barcode': [],
                     'rna_read1': [],
                     'rna_read2': [],
                     'atac_barcode_inclusion_list': '',  # TODO: Add codes for running seqspec
                     'atac_read_format': '',
                     'rna_barcode_inclusion_list': '',
                     'rna_read_format': '',
                     'onlist_mapping': None,
                     'possible_errors': []
                     }
        self.analysis_set_acc = analysis_set_acc
        self.igvf_api = igvf_api
        self.analysis_set_obj = get_analysis_set_obj(
            analysis_set_accessions=[analysis_set_acc], igvf_api=igvf_api)
        self.measet_ids = get_measurement_ids(
            analysis_set_obj=self.analysis_set_obj)
        self.onlist_mapping = get_onlist_mapping_status(
            measurement_set_ids=self.measet_ids, igvf_api=self.igvf_api)

    def get_measet_ids_by_assay_type(self):
        """Get measurement set IDs by assay type.
        """
        measets_sorted_by_assays = get_and_sort_measurement_sets_by_assays(
            measurement_set_ids=self.measet_ids, igvf_api=self.igvf_api)
        # Sort the list of measurement set items by the assay term
        for assay_grp, measet_items in itertools.groupby(measets_sorted_by_assays, key=lambda x: x.assay_term):
            # Get a list of seqfile items, preferred assay title, and linked measurement set accessions.
            assay_measet_accs = [
                measet_item.accession for measet_item in measet_items]
            # Infer assay type based on assay terms
            assay_type = ASSAY_NAMES_CONVERSION_REF[assay_grp]
            # Add the linked measurement set accessions to the report dict
            self.data[f'{assay_type}_MeaSetIDs'].extend(assay_measet_accs)

    def seqfile_urls_for_one_assay(self, measet_items: list, assay_type: str):
        """Get sequence files' urls and seqspec urls from measurement sets of the same assay type.

        Args:
            measet_items (list): list of MeaSet objects
            assay_type (str): rna or atac
        """
        seqfile_items = get_seqfile_items(
            measet_items=measet_items, igvf_api=self.igvf_api)
        # Sort the list of seqfile items by file set acc, read type, run, and lane if present, then
        # get S3 URIs and file accessions.
        for seqfile_item in sorted(seqfile_items, key=seqfile_sort_func):
            if not seqfile_item.read_names:
                continue
            read_file_titles = get_seqfile_readnames(
                seqfile_item=seqfile_item, assay_type=assay_type)
            for read_file_title in read_file_titles:
                if (assay_type == 'rna') and ('barcode' in read_file_title):
                    continue
                # Get FASTQ files' accessions and URLs
                self.data[f'{read_file_title}_accessions'].append(
                    seqfile_item.accession)
                self.data[read_file_title].append(
                    construct_full_href_url(seqfile_item.href))
                # Get Seqspec YAML files' accessions and URLs
                self.data[f'{assay_type}_seqspec_urls'].update(get_seqspec_hrefs(
                    seqspec_ids=seqfile_item.seqspecs, igvf_api=self.igvf_api))

    def get_seqfile_uris_and_seqspec_by_assays_and_runs(self):
        """Get s3 uri of sequence files and add to the report dict based on their assay type and read types.
        """
        # Get a list of measurement set items sorted by the assay terms.
        for assay_type, measet_id in [('rna', self.data['rna_MeaSetIDs']), ('atac', self.data['atac_MeaSetIDs'])]:
            # Get a list of measurement set items sorted by the assay terms.
            measet_items = get_and_sort_measurement_sets_by_assays(
                measurement_set_ids=measet_id, igvf_api=self.igvf_api)
            # Get seqfile URLs and seqspec URLs for each assay type
            self.seqfile_urls_for_one_assay(
                measet_items=measet_items, assay_type=assay_type)

    def single_seqspec_preflight_check(self, seqspec_column: str, onlist_method: str):
        """Check on set of seqspec URLs for the same assay type. If they are all the same, download them and generate the seqspec index and final barcode inclusion list.

        Args:
            seqspec_column (str): A column name in the input table (e.g., atac_seqspec_urls)
            onlist_method (str): Onlist method
        """
        curr_assay_type = seqspec_column.split('_')[0]
        curr_seqspec_urls = list(self.data[seqspec_column])
        if not curr_seqspec_urls:
            self.data['possible_errors'].append(
                f'Error: No seqspec URLs found for {curr_assay_type} seqspecs.')
            return
        try:
            curr_seqspec_file_paths = [
                download_file_via_https(igvf_portal_href_url=seqspec_url)
                for seqspec_url in curr_seqspec_urls
            ]
        except Exception as e:
            self.data['possible_errors'].append(
                f'Error: {curr_assay_type} seqspec download error: {str(e)}')
            return
        # Preflight index check and generation
        try:
            curr_seqspec_index_output = seqspec_index_safetychk_and_get(
                seqspec_file_paths=curr_seqspec_file_paths, assay_type=curr_assay_type)
            self.data[f'{curr_assay_type}_read_format'] = curr_seqspec_index_output
        except Exception as e:
            self.data['possible_errors'].append(
                f'Error: {curr_assay_type} seqspec index generation error: {str(e)}')
        # Preflight onlist check and generation
        try:
            curr_inclusion_list_path = f'./final_barcode_list/{self.analysis_set_acc}_{curr_assay_type}_final_barcode_inclusion_list.txt'
            curr_seqspec_onlist_output = seqspec_onlist_safetychk_and_get(
                seqspec_file_paths=curr_seqspec_file_paths, assay_type=curr_assay_type,
                onlist_method=onlist_method, final_inclusion_list_path=curr_inclusion_list_path)
            self.data[f'{curr_assay_type}_barcode_inclusion_list'] = curr_seqspec_onlist_output
        except Exception as e:
            self.data['possible_errors'].append(str(e))

    def seqspec_preflight_check(self):
        """Check if onlist info and read format info are consistent all across input. If so, generate the read index and final barcode inclusion list.
        """
        atac_onlist_method = get_onlist_method(
            measurement_set_ids=self.data['atac_MeaSetIDs'], igvf_api=self.igvf_api)
        rna_onlist_method = get_onlist_method(
            measurement_set_ids=self.data['rna_MeaSetIDs'], igvf_api=self.igvf_api)
        if atac_onlist_method != rna_onlist_method:
            self.data['possible_errors'].append(
                'Error: RNA and ATAC Measurement sets have different onlist methods.')
        for column, curr_onlist_method in [('rna_seqspec_urls', rna_onlist_method), ('atac_seqspec_urls', atac_onlist_method)]:
            self.single_seqspec_preflight_check(
                seqspec_column=column, onlist_method=curr_onlist_method)

    def get_taxa_and_subpool_info(self):
        """Add sample taxa and sample accession as the subpool ID. Assumption is that one pipeline run, one sample.
        """
        # Assumption is that one pipeline run, one sample
        curr_sample_ids = get_analysis_set_obj(
            analysis_set_accessions=[self.analysis_set_acc], igvf_api=self.igvf_api).samples
        if len(curr_sample_ids) > 1:
            return 'Error: More than one sample in the analysis set.'
        curr_sample_obj = get_sample_obj(
            sample_id=curr_sample_ids[0], igvf_api=self.igvf_api)
        self.data['taxa'] = get_sample_taxa(sample_object=curr_sample_obj)
        self.data['subpool_id'] = curr_sample_obj.accession

    def get_genome_assembly_and_ref_info(self):
        """Add genome reference tsv and assembly info to the input table based on taxa.
        """
        self.data['genome_assembly'] = TAXA_TO_GENOME_REF_FILES[self.data['taxa']][0]
        self.data['genome_ref'] = TAXA_TO_GENOME_REF_FILES[self.data['taxa']][1]

    def build_input_dict(self) -> dict:
        """Build the final output report dict
        """
        self.data['analysis_set_acc'] = self.analysis_set_acc
        self.get_measet_ids_by_assay_type()
        self.get_seqfile_uris_and_seqspec_by_assays_and_runs()
        self.seqspec_preflight_check()
        self.data['atac_seqspec_urls'] = list(self.data['atac_seqspec_urls'])
        self.data['rna_seqspec_urls'] = list(self.data['rna_seqspec_urls'])
        self.get_taxa_and_subpool_info()
        self.get_genome_assembly_and_ref_info()
        self.data['onlist_mapping'] = self.onlist_mapping
        return self.data


def generate_pipeline_input_table(query_analysis_set_accs: list, igvf_api) -> pd.DataFrame:
    """Full script to go from a query measurement set accession to a uniform pipeline input TSV.

    Args:
        query_measet_accs (list): A list of MeaSet accessions for the uniform pipeline
        igvf_api (_type_): _description_

    Returns:
        pd.DataFrame: Final uniform pipeline input tsv.
    """
    pipeline_input_list = []
    for anaset_acc in query_analysis_set_accs:
        # curr_samp_id = get_samp_ids_from_measets(measurement_set_accessions=[measet_acc], igvf_api=igvf_api)
        curr_input_builder = SingleCellInputBuilder(
            analysis_set_acc=anaset_acc, igvf_api=igvf_api)
        pipeline_input_list.append(curr_input_builder.build_input_dict())
    pipeline_input_table = pd.DataFrame(pipeline_input_list)
    return pipeline_input_table


def save_pipeline_input_table(pipeline_input_table: pd.DataFrame, output_dir: str):
    """Save the pipeline input table to a TSV file.

    Args:
        pipeline_input_table (pd.DataFrame): The pipeline input table
        output_dir (str): The output directory
    """
    curr_datetime = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    pipeline_input_table.to_csv(os.path.join(
        output_dir, f'single-cell_uniform_pipeline_input_table_{curr_datetime}.tsv'), sep='\t')
