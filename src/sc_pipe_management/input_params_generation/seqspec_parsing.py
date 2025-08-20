import os
import dataclasses
import requests
import subprocess

# Import seqspec modules
from seqspec.utils import load_spec
from seqspec.seqspec_index import run_index
import seqspec

# Import IGVF API client
import igvf_client

# Import custom modules
import sc_pipe_management.input_params_generation.constant as const
import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing


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


def download_file_via_https(seqspec_file_url: str, partial_root_dir: str) -> str:
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


@dataclasses.dataclass(frozen=True)
class SeqSpecMetadata:
    """Dataclass to hold seqspec metadata."""
    seqspec_file_path: str
    modality: str
    ordered_read_ids: list
    onlist_files: list[str]

    def get_index_read_ids(self) -> str:
        """Get the index read IDs based on the usage purpose."""
        return ','.join(self.ordered_read_ids)

    def get_onlist_read_ids(self) -> str:
        """Get the read IDs for onlist usage purpose."""
        return self.ordered_read_ids[-1]


class GetSeqSpecMetadata:
    """Class to get seqspec metadata from a seqspec file."""

    def __init__(self, seqspec_file_path: str, igvf_api: igvf_client.api.igvf_api.IgvfApi):
        self.seqspec_file_path = seqspec_file_path
        self.igvf_api = igvf_api
        self.imported_seqspec = load_spec(self.seqspec_file_path)

    def _get_seqspec_modality(self) -> str:
        modalities = self.imported_seqspec.modalities
        if len(modalities) != 1:
            raise const.BadDataException(
                f'Error: Seqspec file {self.seqspec_file_path.split("/")[-2]} has multiple modalities. Only one modality is supported.')
        return modalities[0]

    def _find_igvf_acc_in_single_read_spec(self, single_read_spec: seqspec.Read.Read) -> str | None:
        """Find IGVF accession in a single !Read spec. Trying read_id first, then file_ids, and finally the URL."""
        # If read_id is an IGVF accession, return it
        read_id = single_read_spec.read_id
        if const.READ_ID_REGEX.match(read_id):
            return const.READ_ID_REGEX.search(read_id).group(1)
        # If read_id is not an IGVF accession, check the file_ids
        file_specs = single_read_spec.files
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

    def _parse_read_ids_from_seqspec_file(self) -> list[tuple[str, str | None]]:
        """Retrieve read_id values of a specific modality from a seqspec."""
        sequence_specs = self.imported_seqspec.sequence_spec
        read_ids = []
        for spec in sequence_specs:
            if spec.modality == self._get_seqspec_modality():
                igvf_accession = self._find_igvf_acc_in_single_read_spec(
                    single_read_spec=spec)
                read_ids.append((spec.read_id, igvf_accession))
        return read_ids

    def _generate_ordered_read_ids(self, seqfiles_metadata: list[portal_parsing.SeqFileMetadata]) -> list[str]:
        """Get the read1,read2,barcode files order from the seqspec read_id and portal metadata."""
        # Get read IDs from the seqspec file
        read_ids = self._parse_read_ids_from_seqspec_file()
        id_by_name = {}
        for (read_id, igvf_accession) in read_ids:
            # If no IGVF accession is found, add to unknown
            if igvf_accession is None:
                raise const.BadDataException(
                    f'Error: No IGVF accession found for read ID {read_id}.')
            for seqfile_metadata in seqfiles_metadata:
                if seqfile_metadata.file_accession == igvf_accession:
                    read_name = seqfile_metadata.get_read_names_from_seqfile()
                    if read_name in id_by_name:
                        raise const.BadDataException(
                            f'Error: Multiple sequence files found for the same {read_name}.')
                    id_by_name[read_name] = read_id
        return [id_by_name[read_name] for read_name in ['Read 1', 'Read 2', 'Barcode index'] if id_by_name.get(read_name) is not None]

    def _get_onlist_files(self) -> list[str]:
        """Get the onlist files from the seqspec file."""
        curr_run_log = subprocess.run(['seqspec', 'file', '-m', self._get_seqspec_modality(), '-s', 'region-type',
                                      '-i', 'barcode', '-f', 'index', '-k', 'url', self.seqspec_file_path], capture_output=True)
        if curr_run_log.returncode != 0:
            raise const.BadDataException(
                f'Error: seqspec file {self.seqspec_file_path} does not have onlist files.')
        onlist_files = curr_run_log.stdout.decode('utf-8').strip().split('\n')
        return onlist_files

    def generate_seqspec_metadata(self, seqfiles_metadata: list[portal_parsing.SeqFileMetadata]) -> SeqSpecMetadata:
        """Generate seqspec metadata from the seqspec file and seqfiles metadata."""
        # Generate ordered read IDs
        ordered_read_ids = self._generate_ordered_read_ids(seqfiles_metadata)
        return SeqSpecMetadata(
            seqspec_file_path=self.seqspec_file_path,
            modality=self._get_seqspec_modality(),
            ordered_read_ids=ordered_read_ids,
            onlist_files=self._get_onlist_files()
        )


@dataclasses.dataclass(frozen=True)
class SeqSpecToolOutput:
    """Dataclass to hold seqspec tool output."""
    read_index_string: str
    onlist_file_path: str


class SeqSpecMetadataGenerator:
    """Class to generate seqspec metadata for a given seqspec file."""

    def __init__(self, seqspec_metadata: SeqSpecMetadata, partial_root_dir: str):
        self.seqspec_metadata = seqspec_metadata
        self.partial_root_dir = partial_root_dir

    def get_seqspec_index(self) -> str:
        """Get the index read IDs based on the usage purpose."""
        read_index_tool = const.ASSAY_TYPE_TO_TOOL_FORMAT.get(
            self.seqspec_metadata.modality, None)
