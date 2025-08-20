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
from seqspec.seqspec_info import seqspec_info_modalities
import seqspec
import re
import logging
import dataclasses
import igvf_client

import sc_pipe_management.input_params_generation.constant as const
import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing
import sc_pipe_management.igvf_and_terra_api_tools as igvf_api_tools


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
            '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP[read_name]]))
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
                '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    else:
        # NOTE: Hard coded for cases if read2 and barcode are concatenated. Ignore barcode index.
        if sorted(read_names) == sorted(['Read 2', 'Barcode index']):
            read_file_titles.append(
                '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP['Read 2']]))
        else:
            read_file_titles.append(
                '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    return read_file_titles


def get_seqfile_readnames(seqfile_item, assay_type: str) -> str:
    """HACK: Use sequence file metadata to infer FASTQ file read types.

    Args:
        seqfile_item (_type_): Return item from the API
        assay_type (str): translated from assay term into atac or rna

    Returns:
        str: e.g., atac_read1, rna_read2
    """
    read_file_titles = []
    curr_read_names = seqfile_item.read_names
    if not curr_read_names:
        raise const.BadDataException(
            'Error: No read names found in the sequence file item.')
    if assay_type == 'atac':
        return get_atac_seqfile_read_titles(read_names=curr_read_names)
    elif assay_type == 'rna':
        return get_rna_seqfile_read_titles(read_names=curr_read_names)


def download_file_via_https(igvf_download_url: str, igvf_api_keys: dict, output_dir: str = './temp') -> str:
    """Download file via href url on the portal

    Args:
        igvf_download_url (str): _description_
        save_file_path (str): _description_

    Returns:
        str: _description_
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    username = igvf_api_keys['public']
    password = igvf_api_keys['secret']
    session = requests.Session()
    session.auth = (username, password)
    response = session.get(igvf_download_url)
    if response.status_code == 200:
        curr_output_file = os.path.join(
            output_dir, igvf_download_url.split('/')[-1])
        with open(curr_output_file, 'wb') as file:
            file.write(response.content)
        return curr_output_file
    else:
        raise const.BadDataException(
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


def find_igvf_acc_in_seqspec(spec: seqspec.Read.Read) -> str | None:
    """Find IGVF accession in a seqspec for fastq files. Trying read_id first, then file_ids, and finally the URL.

    Args:
        spec (seqspec.Read.Read): The seqspec to search after importing a seqspec.

    Raises:
        BadDataException: If no IGVF accession is found.

    Returns:
        str: The found IGVF accession, parsed if has suffixes.
    """
    # If read_id is an IGVF accession, return it
    read_id = spec.read_id
    if const.READ_ID_REGEX.match(read_id):
        return const.READ_ID_REGEX.search(read_id).group(1)
    # If read_id is not an IGVF accession, check the file_ids
    file_specs = spec.files
    for file_spec in file_specs:
        file_id = file_spec.file_id
        # If the read_id is not an IGVF accession, check the file_id
        if const.READ_ID_REGEX.match(file_id):
            return const.READ_ID_REGEX.search(file_id).group(1)
        # If file_id is not an IGVF accession, check the URL
        for item in file_spec.url.split('/'):
            if const.READ_ID_REGEX.match(item):
                return const.READ_ID_REGEX.search(item).group(1)
    return None


def parse_read_ids_from_seqspec_file(seqspec_file_path: str, assay_type: str) -> list[tuple[str, str | None]]:
    """Retrieve read_id values of a specific modality from a seqspec

    Args:
        seqspec_file_path (str): Local seqspec file path
        assay_type (str): rna or atac

    Raises:
        BadDataException: If the read_id does not contain IGVF accession

    Returns:
        list[tuple[str, str | None]]: a list of (read_id, IGVF accession)
    """
    sequence_specs = load_spec(seqspec_file_path).sequence_spec
    read_ids = []
    for spec in sequence_specs:
        if spec.modality == assay_type:
            igvf_accession = find_igvf_acc_in_seqspec(spec=spec)
            read_ids.append((spec.read_id, igvf_accession))
    return read_ids


def get_read_names_from_seqfile(seqfile_obj) -> str:
    """Get the read names from the seqfile object.

    Args:
        seqfile_obj (_type_): IGVF client query return object for sequence file

    Returns:
        str: The read name (e.g., 'Read 1', 'Read 2', 'Barcode index').
    """
    if len(seqfile_obj.read_names) == 1:
        return seqfile_obj.read_names[0]
    # If Read2 and Barcode are concatenated, use only Read 2
    assert sorted(seqfile_obj.read_names) == sorted(
        ['Read 2', 'Barcode index']), 'Malformed seqfile read names.'
    return 'Read 2'


def generate_ordered_read_ids(seqspec_file_path: str, assay_type: str, usage_purpose: str, igvf_api) -> str:
    """Get the read1,read2,barcode files order from the seqspec read_id and portal metadata.

    Args:
        seqspec_file_path (str): seqspec file path
        assay_type (str): rna or atac
        usage_purpose (str): Usage purpose of the read files, 'onlist' or 'index'. 'onlist' will take a joined list, and 'index' will return either barcode index accession if a separate file or read2 if barcode index is combined with read2 or not indicated.
        igvf_api (_type_): igvf client api

    Raises:
        BadDataException: If no IGVF accession is found for a read ID.
        BadDataException: If sequence files have fatal errors in read_names property.

    Returns:
        str: read1_fastaq_accession,read2_fastaq_accession,barcode_fastaq_accession
    """
    # Get read IDs from the seqspec file
    read_ids = parse_read_ids_from_seqspec_file(
        seqspec_file_path=seqspec_file_path, assay_type=assay_type)
    id_by_name = {}
    for (read_id, igvf_accession) in read_ids:
        # If no IGVF accession is found, add to unknown
        if igvf_accession is None:
            raise const.BadDataException(
                f'Error: No IGVF accession found for read ID {read_id}.')
        # Get seq file object if there is IGVF accession
        seqfile_obj = igvf_api.get_by_id(
            f'/sequence-files/{igvf_accession}/').actual_instance
        # If no read names are found, skip this seqfile
        if not seqfile_obj.read_names:
            continue
        read_name = get_read_names_from_seqfile(seqfile_obj=seqfile_obj)
        if read_name in id_by_name:
            raise const.BadDataException(
                f'Error: Multiple sequence files found for the same {read_name}.')
        id_by_name[read_name] = read_id

    if usage_purpose == 'onlist':
        if 'Barcode index' in id_by_name:
            return id_by_name.get('Barcode index')
        return id_by_name.get('Read 2')
    assert usage_purpose == 'index'
    # Sort the dictionary by the order of keys
    ordered_names = [id_by_name.get('Read 1'), id_by_name.get(
        'Read 2'), id_by_name.get('Barcode index')]
    # for use with seqspec index, it takes all read ids joined
    return ','.join(name for name in ordered_names if name is not None)


def _download_file_via_https(seqspec_file_url: str, partial_root_dir: str) -> str:
    """Download file via href url on the portal and save it locally."""
    if not os.path.exists(partial_root_dir):
        os.makedirs(partial_root_dir)
    # TODO: Needs to update this to use API tools
    username = os.getenv('IGVF_API_KEY')
    password = os.getenv('IGVF_SECRET_KEY')
    session = requests.Session()
    session.auth = (username, password)
    response = session.get(seqspec_file_url)
    if response.status_code == 200:
        curr_output_file = os.path.join(
            partial_root_dir, 'temp', seqspec_file_url.split('/')[-1])
        with open(curr_output_file, 'wb') as file:
            file.write(response.content)
        return curr_output_file
    else:
        raise const.BadDataException(
            f'Error: Download failed with status code: {response.status_code}.')


@dataclasses.dataclass
class SeqSpecMetadata:
    """Dataclass to hold seqspec metadata."""
    seqspec_file_path: str
    modality: str
    read_index_string: str
    index_ordered_read_ids: str
    onlist_read_id: str


class SeqSpecMetadataGenerator:
    """Class to generate seqspec metadata for a given seqspec file."""

    def __init__(self, seqspec_file_url: str, igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str):
        self.seqspec_file_url = seqspec_file_url
        self.igvf_api = igvf_api
        self.partial_root_dir = partial_root_dir

    def _get_modality_from_seqspec(self, seqspec_file_path: str) -> str:
        """Get the modality from the seqspec file. Expecting exactly one modality per seqspec yaml."""
        sequence_specs = load_spec(seqspec_file_path).sequence_spec
        modalities = seqspec_info_modalities(sequence_specs)['modalities']
        if len(modalities) != 1:
            raise const.BadDataException(
                f'Error: Expected exactly one modality in seqspec, found {len(modalities)}.')
        return modalities[0]

    def _find_igvf_acc_in_seqspec(spec: seqspec.Read.Read) -> str | None:
        """Find IGVF accession in a seqspec for fastq files. Trying read_id first, then file_ids, and finally the URL.

        Args:
            spec (seqspec.Read.Read): The seqspec to search after importing a seqspec.

        Raises:
            BadDataException: If no IGVF accession is found.

        Returns:
            str: The found IGVF accession, parsed if has suffixes.
        """
        # If read_id is an IGVF accession, return it
        read_id = spec.read_id
        if const.READ_ID_REGEX.match(read_id):
            return const.READ_ID_REGEX.search(read_id).group(1)
        # If read_id is not an IGVF accession, check the file_ids
        file_specs = spec.files
        for file_spec in file_specs:
            file_id = file_spec.file_id
            # If the read_id is not an IGVF accession, check the file_id
            if const.READ_ID_REGEX.match(file_id):
                return const.READ_ID_REGEX.search(file_id).group(1)
            # If file_id is not an IGVF accession, check the URL
            for item in file_spec.url.split('/'):
                if const.READ_ID_REGEX.match(item):
                    return const.READ_ID_REGEX.search(item).group(1)
        return None

    def parse_read_ids_from_seqspec_file(seqspec_file_path: str, assay_type: str) -> list[tuple[str, str | None]]:
        """Retrieve read_id values of a specific modality from a seqspec

        Args:
            seqspec_file_path (str): Local seqspec file path
            assay_type (str): rna or atac

        Raises:
            BadDataException: If the read_id does not contain IGVF accession

        Returns:
            list[tuple[str, str | None]]: a list of (read_id, IGVF accession)
        """
        sequence_specs = load_spec(seqspec_file_path).sequence_spec
        read_ids = []
        for spec in sequence_specs:
            if spec.modality == assay_type:
                igvf_accession = find_igvf_acc_in_seqspec(spec=spec)
                read_ids.append((spec.read_id, igvf_accession))
        return read_ids
