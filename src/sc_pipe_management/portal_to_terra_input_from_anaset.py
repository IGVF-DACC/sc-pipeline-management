import urllib.parse
import os
import pandas as pd
import itertools
import urllib
from datetime import datetime
import subprocess
import requests
from collections import OrderedDict
from itertools import chain
from seqspec.utils import load_spec
import re


class BadDataException(Exception):
    """Custom exception for bad data."""
    pass

# TODO:
# 1) A SeqFile may have multiple seqspec files but only one is released. Once that is sorted out, need to
#   update the seqspec URLs to be the released one.


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
TAXA_TO_GENOME_REF_FILES = {'Homo sapiens': {'chromap_index': 'https://api.data.igvf.org/reference-files/IGVFFI7969JLFC/@@download/IGVFFI7969JLFC.tar.gz',
                                             'genome_fasta': 'https://api.data.igvf.org/reference-files/IGVFFI0653VCGH/@@download/IGVFFI0653VCGH.fasta.gz',
                                             'kb_index': 'https://api.data.igvf.org/reference-files/IGVFFI9561BASO/@@download/IGVFFI9561BASO.tar.gz',
                                             'genome_ref': 'gs://broad-buenrostro-pipeline-genome-annotations/IGVF_human_v43/IGVF_human_v43_Homo_sapiens_genome_files_hg38_v43.tsv'
                                             },
                            'Mus musculus': {'chromap_index': 'https://api.data.igvf.org/reference-files/IGVFFI5593VLWB/@@download/IGVFFI5593VLWB.tar.gz',
                                             'genome_fasta': 'https://api.data.igvf.org/reference-files/IGVFFI9282QLXO/@@download/IGVFFI9282QLXO.fasta.gz',
                                             'kb_index': 'https://api.data.igvf.org/reference-files/IGVFFI5078MNED/@@download/IGVFFI5078MNED.tar.gz',
                                             'genome_ref': 'gs://broad-buenrostro-pipeline-genome-annotations/IGVF_mouse_v32/IGVF_mouse_v32_Mus_musculus_genome_files_mm39_v32.tsv'
                                             }
                            }

# For URL request
BASE_IGVF_PORTAL_URL = 'https://api.data.igvf.org'

# Read ID Regex
READ_ID_REGEX = re.compile(r'(IGVF[A-Z0-9]*)')

# seqspec index tool (kb-single option will output only one cDNA feature read)
ASSAY_TYPE_TO_TOOL_FORMAT = {'rna': 'kb-single', 'atac': 'chromap'}


def terra_str_formatter(input_strs: list) -> str:
    """Convert Python default array of string output to Terra format (in [] with double quotes)

    Args:
        input_strs (list): a list of URLs, accessions, etc

    Returns:
        str: A string version of the input string list with double quotes
    """
    return '[' + ','.join([f'"{v}"' for v in input_strs]) + ']'


def get_today_mmddyyyy() -> str:
    """Get today's date in MM-DD-YYYY format.

    Returns:
        str: Today's date in MM-DD-YYYY format.
    """
    return datetime.now().strftime('%m%d%Y')


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


# Get analysis set level info
def get_analysis_set_obj(analysis_set_accessions: list, igvf_api):
    """Get the analysis set object from the IGVF API."""
    return igvf_api.analysis_sets(accession=analysis_set_accessions).graph[0]


def get_measurement_set_ids(analysis_set_obj) -> list:
    """Extract measurement set IDs from the analysis set object."""
    return [measet_id for measet_id in analysis_set_obj.input_file_sets if 'measurement-sets' in measet_id]


def get_and_sort_measurement_sets_by_assays(measurement_set_ids: list, igvf_api) -> list:
    """Sort a list of measurement set items based on assay terms.

    Args:
        measurement_set_ids (list): A list of measurement set IDs ('/measurement-sets/IGVFxxxx').
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


def get_atac_seqfile_read_titles(read_names: list) -> str:
    """Get ATACseq read titles. Will always have read1, read2, and barcode index even if read2 and barcode may be concatenated.

    Args:
        read_names (list): seqfile.read_names values

    Returns:
        str: assay-type_read-names, e.g., rna_read1
    """
    assay_type = 'atac'
    read_file_titles = []
    for read_name in read_names:
        read_file_titles.append(
            '_'.join([assay_type, READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    return read_file_titles


def get_rna_seqfile_read_titles(read_names: list) -> str:
    """Get RNA seq file read titles. If reads are not concatenated, will always have read1, read2, and barcode index.
    If reads are concatenated, will only have read1 and read2. The barcode fastaq will be ignored.

    Args:
        read_names (list): seqfile.read_names values

    Returns:
        str: assay-type_read-names, e.g., rna_read1
    """
    assay_type = 'rna'
    read_file_titles = []
    if len(read_names) == 1:
        for read_name in read_names:
            read_file_titles.append(
                '_'.join([assay_type, READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    else:
        # NOTE: Hard coded for cases if read2 and barcode are concatenated. Ignore barcode index.
        if sorted(read_names) == sorted(['Read 2', 'Barcode index']):
            read_file_titles.append(
                '_'.join([assay_type, READ_NAME_TO_READ_TYPE_MAP['Read 2']]))
        else:
            read_file_titles.append(
                '_'.join([assay_type, READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    return read_file_titles


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
        raise BadDataException(
            'Error: No read names found in the sequence file item.')
    if assay_type == 'atac':
        return get_atac_seqfile_read_titles(read_names=curr_read_names)
    elif assay_type == 'rna':
        return get_rna_seqfile_read_titles(read_names=curr_read_names)


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


def get_onlist_method(measurement_set_accs: list, igvf_api) -> str:
    """Get the onlist method of the measurement set items.

    Args:
        measurement_set_accs (list): _description_

    Returns:
        str: onlist method
    """
    curr_onlist_methods = set()
    for measet_acc in measurement_set_accs:
        curr_measet_item = igvf_api.get_by_id(
            f'/measurement-sets/{measet_acc}').actual_instance
        curr_onlist_methods.add(curr_measet_item.onlist_method)
    return list(curr_onlist_methods)[0]


def get_onlist_files_from_measet(measurement_set_accs: list, igvf_api) -> list:
    """Get the onlist files from the measurement sets.

    Args:
        measurement_set_accs (list): _description_

    Returns:
        list: [onlist file URLs]
    """
    onlist_file_urls = []
    for measet_acc in measurement_set_accs:
        curr_measet_item = igvf_api.get_by_id(
            f'/measurement-sets/{measet_acc}').actual_instance
        for onlist_file_id in curr_measet_item.onlist_files:
            onlist_file_item = igvf_api.get_by_id(
                onlist_file_id).actual_instance
            onlist_file_urls.append(
                construct_full_href_url(onlist_file_item.href))
    return onlist_file_urls


def download_file_via_https(igvf_portal_href_url: str, output_dir: str = './temp') -> str:
    """Download file via href url on the portal

    Args:
        igvf_portal_href_url (str): _description_
        save_file_path (str): _description_

    Returns:
        str: _description_
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    username = os.getenv('IGVF_API_KEY')
    password = os.getenv('IGVF_SECRET_KEY')
    session = requests.Session()
    session.auth = (username, password)
    response = session.get(igvf_portal_href_url)
    if response.status_code == 200:
        curr_output_file = os.path.join(
            output_dir, igvf_portal_href_url.split('/')[-1])
        with open(curr_output_file, 'wb') as file:
            file.write(response.content)
        return curr_output_file
    else:
        raise BadDataException(
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


def sort_dict_by_order(input_dict, key_order):
    """Sorts a dictionary based on a provided list specifying the order of keys.

    Args:
      input_dict: The dictionary to sort.
      key_order: A list specifying the desired order of keys.

    Returns:
      A new OrderedDict with keys sorted according to key_order.
      Keys not present in key_order will not be included.
    """
    sorted_dict = OrderedDict()
    for key in key_order:
        if key in input_dict:
            sorted_dict[key] = input_dict[key]
    return sorted_dict


def parse_read_id_from_seqspec(seqspec_file_path: str, assay_type: str) -> str:
    """Retrieve read_id values of a specific modality from a seqspec

    Args:
        seqspec_file_path (str): Local seqspec file path
        assay_type (str): rna or atac

    Raises:
        BadDataException: If the read_id does not contain IGVF accession

    Returns:
        str: read_id value
    """
    sequence_specs = load_spec(seqspec_file_path).sequence_spec
    read_ids = []
    for spec in sequence_specs:
        if spec.modality == assay_type:
            read_id = spec.read_id
            if not read_id.startswith('IGVF'):
                raise BadDataException(
                    f'Error: Read ID {read_id} is not an IGVF accession.')
            read_ids.append(read_id)
    return read_ids


def generate_ordered_read_ids(seqspec_file_path: str, assay_type: str, usage_purpose: str, igvf_api) -> str:
    """Get the read1,read2,barcode files order from the seqspec read_id and portal metadata.

    Args:
        seqspec_file_path (str): seqspec file path
        assay_type (str): rna or atac
        usage_purpose (str): Usage purpose of the read files, 'onlist' or 'index'. 'onlist' will take a joined list, and 'index' will return either barcode index accession if a separate file or read2 if barcode index is combined with read2 or not indicated.
        igvf_api (_type_): igvf client api

    Raises:
        BadDataException: If a seqspec file has multiple sequence files for the same read name.

    Returns:
        str: read1_fastaq_accession,read2_fastaq_accession,barcode_fastaq_accession
    """
    # Get read IDs from the seqspec file
    read_ids = parse_read_id_from_seqspec(
        seqspec_file_path=seqspec_file_path, assay_type=assay_type)
    ordered_seqfiles = {}
    for read_id in read_ids:
        # Parse out IGVF accession from read_id if there is some other formatting
        parsed_read_id = READ_ID_REGEX.search(read_id).group(1)
        # Get seq file object
        seqfile_obj = igvf_api.get_by_id(
            f'/sequence-files/{parsed_read_id}/').actual_instance
        if not seqfile_obj.read_names:
            continue
        # Get the read names from the seqfile object
        # If no fastq concatenation, take all 3 fastqs
        if len(seqfile_obj.read_names) == 1:
            for read_name in seqfile_obj.read_names:
                # Generate read1,read2,barcode index
                ordered_seqfiles.setdefault(read_name, []).append(read_id)
        else:
            # If Read2 and Barcode are concatenated, use only Read 2
            if sorted(seqfile_obj.read_names) == sorted(['Read 2', 'Barcode index']):
                ordered_seqfiles.setdefault(
                    'Read 2', []).append(read_id)
    # If a read name has multiple sequence files, raise an error. This should not happen in a well-formed portal metadata.
    for read_name in ordered_seqfiles:
        if len(ordered_seqfiles[read_name]) > 1:
            raise BadDataException(
                f'Error: Multiple sequence files found for the same {read_name}.')
    # Sort the dictionary by the order of keys
    ordered_seqfiles = sort_dict_by_order(input_dict=ordered_seqfiles, key_order=[
                                          'Read 1', 'Read 2', 'Barcode index'])
    # for use with seqspec index, it takes all read ids joined
    if usage_purpose == 'index':
        return ','.join(list(chain(*ordered_seqfiles.values())))
    # for use with seqspec onlist, return read2_read-id if barcode is either combined or not indicated; return barcode_read-id if it is a separate file.
    elif usage_purpose == 'onlist':
        if not ordered_seqfiles.get('Barcode index'):
            return ordered_seqfiles.get('Read 2')[0]
        else:
            return ordered_seqfiles.get('Barcode index')[0]


def generate_finalinclusion_list(seqspec_file_path: str, assay_type: str, onlist_method: str, final_inclusion_list_path: str, igvf_api) -> str:
    """Generate the final barcode inclusion list txt file.

    Args:
        seqspec_file_path (str): local path to the seqspec file
        assay_type (str): rna or atac
        onlist_method (str): no combination, product, multi
        final_inclusion_list_path (str): where the final onlist file will be saved
        igvf_api (_type_): igvf client api

    Returns:
        str: where the final onlist file will be saved
    """
    # Generate the read id used as -i arg
    seqspec_onlist_readid = generate_ordered_read_ids(
        seqspec_file_path=seqspec_file_path, assay_type=assay_type, igvf_api=igvf_api, usage_purpose='onlist')
    # If no combination, seqspec onlist just repeats the name/url of the only onlist file listed. It needs to be downloaded and renamed.
    if onlist_method == 'no combination':
        curr_run_log = subprocess.run(['seqspec', 'onlist', '-m', assay_type, '-s',
                                      'region-type', '-i', 'barcode', '-o', final_inclusion_list_path, seqspec_file_path], capture_output=True)
    # If combinatorial, seqspec onlist will generated a new file
    else:
        # NOTE: Per sc FG mtg, on Mar 10, 2025, all combinatorial onlist files will be generated as product.
        curr_run_log = subprocess.run(['seqspec', 'onlist', '-m', assay_type, '-s', 'read', '-i',
                                      seqspec_onlist_readid, '-f', 'product', '-o', final_inclusion_list_path, seqspec_file_path], capture_output=True)
    if curr_run_log.returncode != 0:
        raise BadDataException(
            'Error: Seqspec onlist generation command error.')
    if check_emptyinclusion_list(final_inclusion_list_path=final_inclusion_list_path) is False:
        raise BadDataException(
            'Error: Empty final inclusion list file generated.')
    return final_inclusion_list_path


def seqspec_onlist_safetychk_and_get(seqspec_file_paths: list, assay_type: str, measet_onlist_files: list, onlist_method: str, final_inclusion_list_path: str, igvf_api) -> str:
    """Do a pre-flight safety check for seqspec onlist to ensure all seqspec files for the same pipeline have the same onlist info.

    Args:
        seqspec_file_paths (list): a list of locally saved seqspec file paths
        assay_type (str): rna or atac
        measet_onlist_files (list): a list of onlist files from the measurement set items
        onlist_method (str): onlist method (no combination, product, multi)
        final_inclusion_list_path (str): path to save the final inclusion list
        igvf_api (_type_): igvf client api

    Returns:
        str: Either error messages or the path to the first seqspec file in the list.
    """
    curr_onlist_files = []
    for file in seqspec_file_paths:
        curr_run_log = subprocess.run(['seqspec', 'file', '-m', assay_type, '-s', 'region-type',
                                      '-i', 'barcode', '-f', 'index', '-k', 'url', file], capture_output=True)
        if curr_run_log.returncode != 0:
            raise BadDataException('Error: seqspec file command failed.')
        curr_onlist_files = curr_run_log.stdout.decode(
            'utf-8').strip().split(',')
        # If the onlist files are different from the portal metadata, return an error message
        if set(curr_onlist_files) != set(measet_onlist_files):
            raise BadDataException(
                f'Error: Seqspec onlist files are different from the {assay_type.upper()} measurement set onlist files.')

    # If there are multiple onlist files, return error if the onlist method is no combination.
    if (len(curr_onlist_files) > 1) and (onlist_method == 'no combination'):
        raise BadDataException(
            'Error: Multiple onlist files found but onlist method is no combination.')
    # If there is only one onlist file, return error if the onlist method isn't no combination.
    elif (len(curr_onlist_files) == 1) and (onlist_method != 'no combination'):
        raise BadDataException(
            'Error: Only one onlist file found but onlist method is not no combination.')
    # Otherwise, make the final barcode list
    else:
        return generate_finalinclusion_list(seqspec_file_path=str(
            seqspec_file_paths[0]), assay_type=assay_type, onlist_method=onlist_method, final_inclusion_list_path=final_inclusion_list_path, igvf_api=igvf_api)


def seqspec_index_get(seqspec_file_path: str, assay_type: str, igvf_api) -> str:
    """Get RNAseq or ATACseq read index

    Args:
        seqspec_file_path (str): local seqspec file path
        assay_type (str): rna or atac
        igvf_api (_type_): igvf client api

    Raises:
        BadDataException: If seqspec index command fails.

    Returns:
        str: RNAseq or ATACseq read index
    """
    seqspec_index_input_files = generate_ordered_read_ids(
        seqspec_file_path=seqspec_file_path, assay_type=assay_type, igvf_api=igvf_api, usage_purpose='index')
    # Run index with specific files
    curr_run_log = subprocess.run(['seqspec', 'index', '-m', assay_type, '-t',
                                   ASSAY_TYPE_TO_TOOL_FORMAT[assay_type], '-s', 'read', '-i', seqspec_index_input_files, seqspec_file_path], capture_output=True)
    if curr_run_log.returncode != 0:
        raise BadDataException('Error: seqspec index command failed.')
    if assay_type == 'rna':
        return curr_run_log.stdout.decode('utf-8').strip()
    elif assay_type == 'atac':
        return curr_run_log.stdout.decode('utf-8').strip().split(' ')[-1]


def seqspec_index_safetychk_and_get(seqspec_file_paths: list, assay_type: str, igvf_api) -> str:
    """Pre-flight safety check for read_format index and generate it.

    Args:
        seqspec_file_paths (list): A path to the seqspec files to check for read_format consistency.
        assay_type (str): rna or atac

    Returns:
        str: A read format string
    """
    readformat_chk_lst = []
    for file in seqspec_file_paths:
        readformat_chk_lst.append(seqspec_index_get(
            seqspec_file_path=file, assay_type=assay_type, igvf_api=igvf_api))
    if len(set(readformat_chk_lst)) != 1:
        raise BadDataException(
            'Error: Sequence files have different read_format in associated seqspecs.')
    return readformat_chk_lst[0]


def get_onlist_mapping_status(measurement_set_ids: list, igvf_api) -> bool:
    """Get the onlist mapping status of the measurement set items.

    Args:
        measurement_set_ids (list): A list of measurement set IDs (e.g., ['/measurement-sets/IGVFxxxx', ...])

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
        return BadDataException('Error: Measurement sets have different onlist mapping status.')
    else:
        return curr_onlist_mapping[0]


class SingleCellInputBuilder:
    def __init__(self, analysis_set_acc: str, igvf_api):
        self.data = {'analysis_set_acc': '',
                     'atac_MeaSetIDs': [],
                     'rna_MeaSetIDs': [],
                     'subpool_id': '',
                     'taxa': '',
                     'atac_read1_accessions': [],
                     'atac_read2_accessions': [],
                     'atac_barcode_accessions': [],
                     'rna_read1_accessions': [],
                     'rna_read2_accessions': [],
                     'rna_barcode_accessions': [],
                     'atac_seqspec_urls': set(),  # NOTE: The seqspec here is temporary for safety check
                     'rna_seqspec_urls': set(),
                     'atac_read1': [],
                     'atac_read2': [],
                     'atac_barcode': [],
                     'rna_read1': [],
                     'rna_read2': [],
                     'rna_barcode': [],
                     'atac_barcode_inclusion_list': '',
                     'atac_read_format': '',
                     'rna_barcode_inclusion_list': '',
                     'rna_read_format': '',
                     'onlist_mapping': False,
                     'barcode_replacement_file': '',
                     'possible_errors': ''
                     }
        self.analysis_set_acc = analysis_set_acc
        self.igvf_api = igvf_api
        self.analysis_set_obj = get_analysis_set_obj(
            analysis_set_accessions=[analysis_set_acc], igvf_api=igvf_api)
        self.measet_ids = get_measurement_set_ids(
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
            if len(read_file_titles) == 0:
                continue
            for read_file_title in read_file_titles:
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
        for assay_type, measet_accs in [('rna', self.data['rna_MeaSetIDs']), ('atac', self.data['atac_MeaSetIDs'])]:
            # Convert measet accessions to measurement set IDs
            measet_ids = [
                f'/measurement-sets/{measet_acc}' for measet_acc in measet_accs]
            # Get a list of measurement set items sorted by the assay terms.
            measet_items = get_and_sort_measurement_sets_by_assays(
                measurement_set_ids=measet_ids, igvf_api=self.igvf_api)
            # Get seqfile URLs and seqspec URLs for each assay type
            self.seqfile_urls_for_one_assay(
                measet_items=measet_items, assay_type=assay_type)

    def single_seqspec_preflight_check(self, seqspec_column: str, onlist_method: str, measurement_set_accs: list, local_barcode_file_dir: str):
        """Check on set of seqspec URLs for the same assay type. If they are all the same, download them and generate the seqspec index and final barcode inclusion list.

        Args:
            seqspec_column (str): A column name in the input table (e.g., atac_seqspec_urls)
            onlist_method (str): Onlist method
            measurement_set_accs (list): A list of measurement set accessions
            local_barcode_file_dir (str, optional): Directory to save the final inclusion list.
        """
        curr_assay_type = seqspec_column.split('_')[0]
        curr_seqspec_urls = list(self.data[seqspec_column])
        if not curr_seqspec_urls:
            self.data['possible_errors'] += f'Error: No seqspec URLs found for {curr_assay_type} seqspecs.'
            return
        try:
            curr_seqspec_file_paths = [
                download_file_via_https(igvf_portal_href_url=seqspec_url)
                for seqspec_url in curr_seqspec_urls
            ]
        except BadDataException as e:
            self.data[
                'possible_errors'] += f'Error: {curr_assay_type} seqspec download error: {str(e)}.'
            return
        # Preflight index check and generation
        try:
            curr_seqspec_index_output = seqspec_index_safetychk_and_get(
                seqspec_file_paths=curr_seqspec_file_paths, assay_type=curr_assay_type, igvf_api=self.igvf_api)
            self.data[f'{curr_assay_type}_read_format'] = curr_seqspec_index_output
        except BadDataException as e:
            self.data[
                'possible_errors'] += f'Error: {curr_assay_type} seqspec index generation error: {str(e)}.'
        # Preflight onlist check and generation
        try:
            curr_measet_onlist_files = get_onlist_files_from_measet(
                measurement_set_accs=measurement_set_accs, igvf_api=self.igvf_api)
        except BadDataException as e:
            self.data['possible_errors'] += str(e)
            return
        try:
            # Dir will be a fixed place with a date folder under the current working directory
            # final_inclusion_list_dir = os.path.join(
            #     os.getcwd(), 'final_barcode_list', get_today_mmddyyyy())
            if not os.path.exists(local_barcode_file_dir):
                os.makedirs(local_barcode_file_dir)
            curr_inclusion_list_path = os.path.join(
                local_barcode_file_dir, f'{self.analysis_set_acc}_{curr_assay_type}_final_barcode_inclusion_list.txt')
            curr_seqspec_onlist_output = seqspec_onlist_safetychk_and_get(
                seqspec_file_paths=curr_seqspec_file_paths, assay_type=curr_assay_type, measet_onlist_files=curr_measet_onlist_files,
                onlist_method=onlist_method, final_inclusion_list_path=curr_inclusion_list_path, igvf_api=self.igvf_api)
            self.data[f'{curr_assay_type}_barcode_inclusion_list'] = curr_seqspec_onlist_output
        except BadDataException as e:
            self.data['possible_errors'] += str(e)

    def seqspec_preflight_check(self, local_barcode_file_dir: str):
        """Check if onlist info and read format info are consistent all across input. If so, generate the read index and final barcode inclusion list.

        Args:
            local_barcode_file_dir (str, optional): Directory to save the final inclusion list.
        """
        for column, curr_measet_accs in [('rna_seqspec_urls', self.data['rna_MeaSetIDs']), ('atac_seqspec_urls', self.data['atac_MeaSetIDs'])]:
            if not curr_measet_accs:
                continue
            onlist_method = get_onlist_method(
                measurement_set_accs=curr_measet_accs, igvf_api=self.igvf_api)
            self.single_seqspec_preflight_check(
                seqspec_column=column, onlist_method=onlist_method, measurement_set_accs=curr_measet_accs, local_barcode_file_dir=local_barcode_file_dir)

    def get_taxa_and_subpool_info(self):
        """Add sample taxa and sample accession as the subpool ID. Assumption is that one pipeline run, one sample.
        """
        # Assumption is that one pipeline run, one sample
        curr_sample_ids = get_analysis_set_obj(
            analysis_set_accessions=[self.analysis_set_acc], igvf_api=self.igvf_api).samples
        # Anticipating multiple samples in an analysis set
        sample_taxa = set()
        sample_subpools = set()
        for curr_sample_id in curr_sample_ids:
            curr_sample_obj = get_sample_obj(
                sample_id=curr_sample_id, igvf_api=self.igvf_api)
            sample_taxa.add(curr_sample_obj.taxa)
            sample_subpools.add(curr_sample_obj.accession)
        # Assumption is that one pipeline run, one taxa
        if len(sample_taxa) > 1:
            self.data['possible_errors'] += 'Error: Multiple taxa found in the same analysis set.'
        self.data['taxa'] = list(sample_taxa)[0]
        # Per Eugenio, it is OK to use hyphen
        self.data['subpool_id'] = '-'.join(sorted(sample_subpools))

    def get_ref_file_info(self):
        """Add genome reference tsv and assembly info to the input table based on taxa.
        """
        for ref_src, igvf_href in TAXA_TO_GENOME_REF_FILES[self.data['taxa']].items():
            self.data[ref_src] = igvf_href

    def get_barcode_replacement_file(self):
        """Get the barcode replacement file if it exists in the measurement set.
        """
        # NOTE: Assumption that all measurement sets in the analysis set have the same barcode replacement file because an audit will show up if not the case.
        # Take one of the RNA measurement set accessions
        # Temporary solution until the pipeline is updated to use the tab file object
        replacement_file_gs_urls = {
            'IGVFFI4834TCML': 'gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/barcode_replacement_files/IGVFFI4834TCML_r1_RT_replace.tsv',
            'IGVFFI4529PWCF': 'gs://fc-secure-de19fd29-2253-41cd-9751-1788cf7ad1a5/barcode_replacement_files/IGVFFI4529PWCF_r1_RT_replace_wt.tsv'
        }
        if self.data['rna_MeaSetIDs']:
            measet_obj = self.igvf_api.get_by_id(
                f'/measurement-sets/{self.data["rna_MeaSetIDs"][0]}').actual_instance
            # If the barcode replacement file exists, get the URL from the tab file object
            if measet_obj.barcode_replacement_file:
                tab_file_obj = self.igvf_api.get_by_id(
                    measet_obj.barcode_replacement_file).actual_instance
                # Once pipeline is updated
                self.data['barcode_replacement_file'] = construct_full_href_url(
                    igvf_href=tab_file_obj.href)
                # self.data['barcode_replacement_file'] = [
                #     replacement_file_gs_urls[tab_file_obj.accession]]

    def reformat_arrays_to_terra_format(self):
        """Reformat all read file URLs arrays to Terra format.
        """
        for col in ['atac_read1', 'atac_read2', 'atac_barcode', 'rna_read1', 'rna_read2', 'rna_barcode']:
            self.data[col] = terra_str_formatter(self.data[col])
        for col in self.data.keys():
            # Per Sid, use no value of any sort for barcode_replacement_file if not applicable
            if col != 'barcode_replacement_file':
                # Empty list: force empty list as '[]', double quotes to avoid JSON parsing issues
                if self.data[col] == []:
                    self.data[col] = "[]"
                # Empty string: force empty string as [''], double quotes to avoid JSON parsing issues
                elif self.data[col] == '':
                    self.data[col] = "None"

    def build_input_dict(self, local_barcode_file_dir: str) -> dict:
        """Build the final output report dict

        Args:
            local_barcode_file_dir (str): Directory to save the final inclusion list.
        """
        # Add analysis set accession
        self.data['analysis_set_acc'] = self.analysis_set_acc
        # Add measurement set accessions by assay type
        self.get_measet_ids_by_assay_type()
        # Add sequence file URLs and seqspec URLs by assay type
        self.get_seqfile_uris_and_seqspec_by_assays_and_runs()
        # Safety check seqspec
        self.seqspec_preflight_check(
            local_barcode_file_dir=local_barcode_file_dir)
        # Clean up seqspec urls
        self.data['atac_seqspec_urls'] = sorted(self.data['atac_seqspec_urls'])
        self.data['rna_seqspec_urls'] = sorted(self.data['rna_seqspec_urls'])
        # Get subpool and taxa info
        self.get_taxa_and_subpool_info()
        # Get ref file info
        self.get_ref_file_info()
        # Get onlist mapping bool value
        self.data['onlist_mapping'] = self.onlist_mapping
        # Get barcode replacement file if it exists
        self.get_barcode_replacement_file()
        # Reformat the table for Terra
        self.reformat_arrays_to_terra_format()
        self.data['kb_strand'] = 'forward'  # currently hard coded
        return self.data


def mod_input_table_for_terra(pipeline_input_table: pd.DataFrame, terra_etype: str, local_barcode_file_dir: str, gs_barcode_list_bucket: str) -> pd.DataFrame:
    """Modify the pipeline input table to have Terra table ID in column name and replace local file path with gs:// bucket path.

    Args:
        pipeline_input_table (pd.DataFrame): The input table generated from IGVF portal without any modification
        terra_etype (str): Terra entity type name
        local_barcode_file_path (str): Final barcode file directory
        gs_barcode_list_bucket (str): Where the barcode list files will be on gcloud

    Returns:
        pd.DataFrame: The reformatted Terra pipeline input table
    """
    pipeline_input_table.index = pipeline_input_table['analysis_set_acc']
    pipeline_input_table.index.name = f'entity:{terra_etype}_id'
    pipeline_input_table_withgs = pipeline_input_table.replace(
        local_barcode_file_dir.rstrip('/'), gs_barcode_list_bucket.rstrip('/'), regex=True)
    return pipeline_input_table_withgs


def generate_pipeline_input_table(query_analysis_set_accs: list, igvf_api, terra_etype: str, local_barcode_file_dir: str, gs_barcode_list_bucket: str) -> pd.DataFrame:
    """Full script to go from a query measurement set accession to a uniform pipeline input TSV.

    Args:
        query_measet_accs (list): A list of MeaSet accessions for the uniform pipeline
        igvf_api (_type_): _description_
        terra_etype (str): Terra entity type name
        local_barcode_file_path (str): Final barcode file directory
        gs_barcode_list_bucket (str): Where the barcode list files will be on gcloud

    Returns:
        pd.DataFrame: Final uniform pipeline input tsv compatible with Terra syntax.
    """
    pipeline_input_list = []
    for anaset_acc in query_analysis_set_accs:
        print('Processing:', anaset_acc)
        curr_input_builder = SingleCellInputBuilder(
            analysis_set_acc=anaset_acc, igvf_api=igvf_api)
        pipeline_input_list.append(curr_input_builder.build_input_dict(
            local_barcode_file_dir=local_barcode_file_dir))
        print('Done:', anaset_acc)
    pipeline_input_table = pd.DataFrame(pipeline_input_list)
    print('Reformatting input table for Terra format...')
    print('Done.')
    return mod_input_table_for_terra(pipeline_input_table=pipeline_input_table,
                                     terra_etype=terra_etype,
                                     local_barcode_file_dir=local_barcode_file_dir,
                                     gs_barcode_list_bucket=gs_barcode_list_bucket)


def save_pipeline_input_table(pipeline_input_table: pd.DataFrame, output_dir: str) -> pd.DataFrame:
    """Save the pipeline input table to a TSV file.

    Args:
        pipeline_input_table (pd.DataFrame): The pipeline input table
        output_dir (str): The output directory
    """
    curr_datetime = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir)
    # Save the pipeline input table to a TSV file
    pipeline_input_table.to_csv(os.path.join(
        output_dir, f'single-cell_uniform_pipeline_input_table_{curr_datetime}_withGSpath.tsv'), sep='\t')
