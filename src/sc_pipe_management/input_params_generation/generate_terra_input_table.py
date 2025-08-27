import pandas as pd
import os
import requests
import dataclasses
import igvf_client


import sc_pipe_management.input_params_generation as const
import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing
import sc_pipe_management.input_params_generation.seqspec_parsing as seqspec_parsing


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

    def _check_emptyinclusion_list(self) -> bool:
        """Check if the final inclusion list created is empty."""
        with open(self.output_barcode_list_file, 'r') as file_obj:
            first_char = file_obj.read(1)
            if not first_char.isalpha():
                return False


class QualityCheckInputData:
    """Class to hold quality check input data."""

    def __init__(self, analysis_set_metadata: portal_parsing.AnalysisSetMetadata, igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str):
        self.analysis_set_metadata = analysis_set_metadata
        self.rna_input_info = self.analysis_set_metadata.rna_input_info
        self.atac_input_info = self.analysis_set_metadata.atac_input_info
        self.igvf_api = igvf_api
        self.partial_root_dir = partial_root_dir

    def _get_seqspec_accession_from_path(self, seqspec_file_path: str) -> str:
        """Get the seqspec accession from the seqspec file path."""
        return const.READ_ID_REGEX.search(seqspec_file_path).group(1)

    def _get_all_seqspec_metadata_per_assay_type(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata]) -> list[seqspec_parsing.SeqSpecMetadata]:
        """Get all seqspec metadata for a given assay type."""
        all_seqspec_metadata = []
        for measet_metadata in measet_metadata_list:
            unique_seqspecs = set()
            for seqfile_item in measet_metadata.seqfile_items:
                if seqfile_item.seqspec_urls:
                    unique_seqspecs.update(seqfile_item.seqspec_urls)
            # Download and parse each unique seqspec file
            for seqspec_url in list(unique_seqspecs):
                # Download the seqspec file
                curr_seqspec_file_path = download_file_via_https(
                    seqspec_file_url=seqspec_url, partial_root_dir=self.partial_root_dir)
                # Get the seqspec metadata
                curr_seqspec_metadata = seqspec_parsing.GetSeqSpecMetadata(
                    seqspec_file_path=curr_seqspec_file_path, igvf_api=self.igvf_api).generate_seqspec_metadata()
                # Save it to the list
                all_seqspec_metadata.append(curr_seqspec_metadata)
        return all_seqspec_metadata

    def _get_all_seqspec_tool_outputs_per_assay_type(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata]) -> list[seqspec_parsing.SeqSpecToolOutputs]:
        """Get all seqspec tool outputs for a given assay type."""
        all_seqspec_tool_outputs = []
        for seqspec_metadata in all_seqspec_metadata:
            curr_seqspec_tool_outputs = seqspec_parsing.SeqSpecToolOutputs(
                seqspec_metadata=seqspec_metadata).generate_seqspec_tool_outputs()
            all_seqspec_tool_outputs.append(curr_seqspec_tool_outputs)
        return all_seqspec_tool_outputs

    def check_if_modality_match(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata], expected_modality: str) -> bool:
        """Check if all seqspec modalities match the expected modality."""
        errors = []
        for seqspec_metadata in all_seqspec_metadata:
            if seqspec_metadata.modality != expected_modality:
                error_seqspec_accession = self._get_seqspec_accession_from_path(
                    seqspec_metadata.seqspec_file_path)
                errors.append(
                    f"Modality mismatch: expected '{expected_modality}', got '{seqspec_metadata.modality}' from seqspec file {error_seqspec_accession} instead.")
        if errors:
            raise const.BadDataException(
                "Modality mismatches found:\n" + "\n".join(errors))
        return True

    def check_if_onlist_match(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata], expected_onlist_files: list[str]) -> bool:
        """Check if all seqspec onlist methods match the expected onlist method."""
        errors = []
        for seqspec_metadata in all_seqspec_metadata:
            unique_onlist_files = sorted(set(seqspec_metadata.onlist_files))
            if unique_onlist_files != sorted(expected_onlist_files):
                error_seqspec_accession = self._get_seqspec_accession_from_path(
                    seqspec_metadata.seqspec_file_path)
                errors.append(
                    f"Onlist files mismatch: expected {sorted(expected_onlist_files)}, got {unique_onlist_files} from seqspec file {error_seqspec_accession} instead.")
        if errors:
            raise const.BadDataException(
                "Onlist method mismatches found:\n" + "\n".join(errors))
        return True

    def check_if_read_index_match(self, all_seqspec_tool_outputs: list[seqspec_parsing.SeqSpecToolOutputs]) -> bool:
        """Check if all seqspec read index strings match the expected read index string."""
        errors = []
        for seqspec_tool_output in all_seqspec_tool_outputs:
            if seqspec_tool_output.read_index_string != expected_read_index:
                error_seqspec_accession = self._get_seqspec_accession_from_path(
                    seqspec_tool_output.seqspec_metadata.seqspec_file_path)
                errors.append(
                    f"Read index string mismatch: expected '{expected_read_index}', got '{seqspec_tool_output.read_index_string}' from seqspec file {error_seqspec_accession} instead.")
        if errors:
            raise const.BadDataException(
                "Read index string mismatches found:\n" + "\n".join(errors))
        return True


@dataclasses.dataclass(frozen=True)
class TerraPipelineParams:
    analysis_set_acc: str = ''
    # Measurement Set accessions
    atac_MeaSetIDs: list[str] = dataclasses.field(default_factory=list)
    rna_MeaSetIDs: list[str] = dataclasses.field(default_factory=list)
    # Sample IDs
    subpool_id: str = ''
    # Taxa (for filtering reference files)
    taxa: str = ''
    # FASTQ file accessions (for accessioning, not used in the pipeline)
    atac_read1_accessions: list[str] = dataclasses.field(default_factory=list)
    atac_read2_accessions: list[str] = dataclasses.field(default_factory=list)
    atac_barcode_accessions: list[str] = dataclasses.field(
        default_factory=list)
    rna_read1_accessions: list[str] = dataclasses.field(default_factory=list)
    rna_read2_accessions: list[str] = dataclasses.field(default_factory=list)
    rna_barcode_accessions: list[str] = dataclasses.field(default_factory=list)
    atac_seqspec_urls: set[str] = dataclasses.field(default_factory=set)
    rna_seqspec_urls: set[str] = dataclasses.field(default_factory=set)
    # FASTQ file URLs and Parse-specific barcode file (for running the pipeline)
    atac_read1: list[str] = dataclasses.field(default_factory=list)
    atac_read2: list[str] = dataclasses.field(default_factory=list)
    atac_barcode: list[str] = dataclasses.field(default_factory=list)
    rna_read1: list[str] = dataclasses.field(default_factory=list)
    rna_read2: list[str] = dataclasses.field(default_factory=list)
    rna_barcode: list[str] = dataclasses.field(default_factory=list)
    barcode_replacement_file: str = ''
    # Seqspec tool outputs (used in the pipeline)
    atac_barcode_inclusion_list: str = ''
    atac_read_format: str = ''
    rna_barcode_inclusion_list: str = ''
    rna_read_format: str = ''
    # Other params
    onlist_mapping: bool = False
    kb_strand: str = 'forward'
    chromap_index: str = ''
    genome_fasta: str = ''
    kb_index: str = ''
    genome_ref: str = ''
    # Error reporting
    possible_errors: str = ''
