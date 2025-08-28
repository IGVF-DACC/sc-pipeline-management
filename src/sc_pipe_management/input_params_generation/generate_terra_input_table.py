import pandas as pd
import os
import requests
import dataclasses
import igvf_client


import sc_pipe_management.input_params_generation.constant as const
import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing
import sc_pipe_management.input_params_generation.seqspec_parsing as seqspec_parsing


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


def get_seqspec_accession_from_path(seqspec_file_path: str) -> str:
    """Get the seqspec accession from the seqspec file path."""
    return const.READ_ID_REGEX.search(seqspec_file_path).group(1)


class QCandParseSeqspecs:
    """Class method for checking seqspecs against portal metadata."""

    def __init__(self, analysis_set_metadata: portal_parsing.InputAnalysisSetMetadata, igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str):
        self.analysis_set_metadata = analysis_set_metadata
        self.analysis_set_acc = self.analysis_set_metadata.analysis_set_acc
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
            for seqfile_item in measet_metadata.seqfiles:
                if seqfile_item.seqspec_urls:
                    unique_seqspecs.update(seqfile_item.seqspec_urls)
            # Download and parse each unique seqspec file
            for seqspec_url in list(unique_seqspecs):
                # Download the seqspec file
                curr_seqspec_file_path = download_file_via_https(
                    seqspec_file_url=seqspec_url, igvf_api=self.igvf_api, partial_root_dir=self.partial_root_dir)
                # Get the seqspec metadata
                curr_seqspec_metadata = seqspec_parsing.GetSeqSpecMetadata(
                    seqspec_file_path=curr_seqspec_file_path, igvf_api=self.igvf_api).generate_seqspec_metadata(seqfiles_metadata=measet_metadata.seqfiles)
                # Save it to the list
                all_seqspec_metadata.append(curr_seqspec_metadata)
        return all_seqspec_metadata

    def _get_all_seqspec_tool_outputs_per_assay_type(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata], onlist_method: str) -> list[seqspec_parsing.SeqSpecToolOutput]:
        """Get all seqspec tool outputs for a given assay type."""
        all_seqspec_tool_outputs = []
        for seqspec_metadata in all_seqspec_metadata:
            curr_final_barcode_file_name = '_'.join([self.analysis_set_acc, seqspec_metadata.modality,
                                                     get_seqspec_accession_from_path(
                                                         seqspec_file_path=seqspec_metadata.seqspec_file_path),
                                                     'final_inclusion_list.txt'])
            curr_seqspec_tool_outputs = seqspec_parsing.GetSeqSpecToolOutput(
                seqspec_metadata=seqspec_metadata,
                onlist_method=onlist_method,
                output_barcode_list_file=os.path.join(self.partial_root_dir, curr_final_barcode_file_name)).generate_seqspec_tool_output()
            all_seqspec_tool_outputs.append(curr_seqspec_tool_outputs)
        return all_seqspec_tool_outputs

    def check_if_modality_match(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata], expected_modality: str) -> bool:
        """Check if all seqspec modalities match the expected modality."""
        errors = []
        for seqspec_metadata in all_seqspec_metadata:
            if seqspec_metadata.modality != expected_modality:
                error_seqspec_accession = get_seqspec_accession_from_path(
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
            # Sort and parse IGVF accessions from the onlist files in test files
            unique_onlist_file_accs = [get_seqspec_accession_from_path(
                entry) for entry in sorted(set(seqspec_metadata.onlist_files))]
            # Sort and parse IGVF accessions from the expected onlist files
            expected_onlist_file_accs = [get_seqspec_accession_from_path(
                entry) for entry in sorted(set(expected_onlist_files))]
            # Compare the two lists
            if unique_onlist_file_accs != expected_onlist_file_accs:
                error_seqspec_accession = get_seqspec_accession_from_path(
                    seqspec_metadata.seqspec_file_path)
                errors.append(
                    f"Onlist files mismatch: expected {expected_onlist_file_accs}, got {unique_onlist_file_accs} from seqspec file {error_seqspec_accession} instead.")
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
                error_seqspec_accession = get_seqspec_accession_from_path(
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
        reference_read_index = all_seqspec_tool_outputs[0].read_index_string
        for seqspec_tool_output in all_seqspec_tool_outputs:
            if seqspec_tool_output.read_index_string != reference_read_index:
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
            all_seqspec_metadata = self._get_all_seqspec_metadata_per_assay_type(
                measet_metadata_list=measet_metadata_list)
            # Assuming all RNA measets have the same onlist files (or else the portal will audit)
            expected_onlist_files = measet_metadata_list[0].onlist_files
            expected_onlist_method = measet_metadata_list[0].onlist_method
            # Check modality (fixed to rna)
            self.check_if_modality_match(
                all_seqspec_metadata=all_seqspec_metadata, expected_modality=assay_type.lower())
            # Check onlist files
            self.check_if_onlist_match(
                all_seqspec_metadata=all_seqspec_metadata, expected_onlist_files=expected_onlist_files)
            # Check onlist method
            self.check_if_onlist_method_match(
                all_seqspec_metadata=all_seqspec_metadata, expected_onlist_method=expected_onlist_method)
            # Get all seqspec tool outputs (it assumes that all measurement sets should have the same onlist method)
            # Portal checks for inconsistent onlist methods across measurement sets of the same assay type
            all_seqspec_tool_outputs = self._get_all_seqspec_tool_outputs_per_assay_type(
                all_seqspec_metadata=all_seqspec_metadata, onlist_method=expected_onlist_method)
            # Check read index strings
            self.check_if_read_index_match(
                all_seqspec_tool_outputs=all_seqspec_tool_outputs, assay_type=assay_type.upper())
            # Check empty inclusion lists
            self.check_if_empty_inclusion_list(
                all_seqspec_tool_outputs=all_seqspec_tool_outputs, assay_type=assay_type.upper())
        except const.BadDataException as e:
            possible_errors += str(e)
        if possible_errors:
            # This syntax is to satisfy Terra input table format
            return seqspec_parsing.SeqSpecToolOutput(read_index_string="None", final_barcode_file="None", errors=possible_errors)
        return all_seqspec_tool_outputs[0]

    def quality_check_input_data(self) -> dict:
        """Quality check the input analysis set metadata."""
        seqspecs_for_input_params = {'rna': None, 'atac': None}
        if self.rna_input_info:
            seqspecs_for_input_params['rna'] = self._quality_check_one_assay_type(
                measet_metadata_list=self.rna_input_info, assay_type='rna')
        if self.atac_input_info:
            seqspecs_for_input_params['atac'] = self._quality_check_one_assay_type(
                measet_metadata_list=self.atac_input_info, assay_type='atac')
        return seqspecs_for_input_params


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
