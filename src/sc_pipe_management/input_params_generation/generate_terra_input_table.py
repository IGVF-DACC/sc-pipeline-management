import pandas as pd
import os
import requests
import dataclasses
import igvf_client


import sc_pipe_management.input_params_generation as const
import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing
import sc_pipe_management.input_params_generation.seqspec_parsing as seqspec_parsing
import sc_pipe_management.igvf_and_terra_api_tools as api_tools


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


def download_file_via_https(seqspec_file_url: str, igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str) -> str:
    """Download file via href url on the portal and save it locally."""
    download_folder = os.path.join(partial_root_dir, 'temp')
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)
    igvf_api_config = igvf_api.api_client.configuration
    session = requests.Session()
    session.auth = (igvf_api_config.username, igvf_api_config.password)
    response = session.get(seqspec_file_url)
    if response.status_code == 200:
        curr_output_file = os.path.join(
            download_folder, seqspec_file_url.split('/')[-1])
        with open(curr_output_file, 'wb') as file:
            file.write(response.content)
        return curr_output_file
    else:
        raise const.BadDataException(
            f'Error: Download failed with status code: {response.status_code}.')


class QualityCheckInputData:
    """Class method for checking seqspecs against portal metadata."""

    def __init__(self, analysis_set_metadata: portal_parsing.InputAnalysisSetMetadata, igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str):
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

    def _get_all_seqspec_tool_outputs_per_assay_type(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata]) -> list[seqspec_parsing.SeqSpecToolOutput]:
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

    def check_if_onlist_method_match(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata], expected_onlist_method: str) -> bool:
        """Check if all seqspec onlist methods match the expected onlist method."""
        errors = []
        for seqspec_metadata in all_seqspec_metadata:
            # Check 1: if expected is 'no combination', then there should be only one onlist file
            error_condition_1 = (len(seqspec_metadata.onlist_files) > 1) and (
                expected_onlist_method == 'no combination')
            # Check 2: if expected is combinatorial, then there should be more than one onlist file
            error_condition_2 = (len(seqspec_metadata.onlist_files) == 1) and (
                expected_onlist_method != 'no combination')
            # Collect errors
            if error_condition_1 or error_condition_2:
                error_seqspec_accession = self._get_seqspec_accession_from_path(
                    seqspec_metadata.seqspec_file_path)
                errors.append(
                    f"Onlist method mismatch: expected '{expected_onlist_method}', got {len(seqspec_metadata.onlist_files)} from seqspec file {error_seqspec_accession} instead.")
        if errors:
            raise const.BadDataException(
                "Onlist method mismatches found:\n" + "\n".join(errors))
        return True

    def check_if_read_index_match(self, all_seqspec_tool_outputs: list[seqspec_parsing.SeqSpecToolOutput], assay_type: str) -> bool:
        """Check if all seqspec read index strings match the first read index string."""
        mismatch_cnt = 0
        # Use the first read index string as the reference
        reference_read_index = all_seqspec_tool_outputs[0].read_index
        for seqspec_tool_output in all_seqspec_tool_outputs:
            if seqspec_tool_output.read_index != reference_read_index:
                mismatch_cnt += 1
        if mismatch_cnt > 0:
            raise const.BadDataException(
                f"{mismatch_cnt} {assay_type} read index string mismatches found.")
        return True

    def check_if_empty_inclusion_list(self, all_seqspec_tool_outputs: list[seqspec_parsing.SeqSpecToolOutput], assay_type: str) -> bool:
        """Check if any of the final inclusion lists created is empty."""
        empty_cnt = 0
        for seqspec_tool_output in all_seqspec_tool_outputs:
            with open(seqspec_tool_output.final_barcode_file, 'r') as file_obj:
                first_char = file_obj.read(1)
                if not first_char.isalpha():
                    empty_cnt += 1
        if empty_cnt > 0:
            raise const.BadDataException(
                f"{empty_cnt} {assay_type} final inclusion lists are empty.")
        return True

    def _quality_check_one_assay_type(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata], assay_type: str) -> str | None:
        """Quality check the input analysis set metadata."""
        possible_errors = ''
        try:
            rna_all_seqspec_metadata = self._get_all_seqspec_metadata_per_assay_type(
                measet_metadata_list=self.rna_input_info)
            # Assuming all RNA measets have the same onlist files (or else the portal will audit)
            expected_onlist_files = self.rna_input_info[0].onlist_files
            expected_onlist_method = self.rna_input_info[0].onlist_method
            # Check modality (fixed to rna)
            self.check_if_modality_match(
                all_seqspec_metadata=rna_all_seqspec_metadata, expected_modality=assay_type.lower())
            # Check onlist files
            self.check_if_onlist_match(
                all_seqspec_metadata=rna_all_seqspec_metadata, expected_onlist_files=expected_onlist_files)
            # Check onlist method
            self.check_if_onlist_method_match(
                all_seqspec_metadata=rna_all_seqspec_metadata, expected_onlist_method=expected_onlist_method)
            # Get all seqspec tool outputs
            rna_all_seqspec_tool_outputs = self._get_all_seqspec_tool_outputs_per_assay_type(
                all_seqspec_metadata=rna_all_seqspec_metadata)
            # Check read index strings
            self.check_if_read_index_match(
                all_seqspec_tool_outputs=rna_all_seqspec_tool_outputs, assay_type=assay_type.upper())
            # Check empty inclusion lists
            self.check_if_empty_inclusion_list(
                all_seqspec_tool_outputs=rna_all_seqspec_tool_outputs, assay_type=assay_type.upper())
        except const.BadDataException as e:
            possible_errors += str(e)
        if possible_errors:
            return possible_errors
        return None

    def quality_check_input_data(self) -> str:
        """Quality check the input analysis set metadata."""
        all_possible_errors = ''
        if self.rna_input_info:
            rna_possible_errors = self._quality_check_one_assay_type(
                measet_metadata_list=self.rna_input_info, assay_type='rna')
            if rna_possible_errors:
                all_possible_errors += rna_possible_errors
        if self.atac_input_info:
            atac_possible_errors = self._quality_check_one_assay_type(
                measet_metadata_list=self.atac_input_info, assay_type='atac')
            if atac_possible_errors:
                all_possible_errors += atac_possible_errors
        if all_possible_errors:
            return all_possible_errors
        # This is specifically for Terra, which uses "None" to indicate no value for string fields
        return "None"


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
