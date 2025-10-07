import pandas as pd
import os
import requests
import dataclasses
import igvf_client
import datetime
import logging


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
    return const.IGVF_ACCESSION_REGEX.search(seqspec_file_path).group(1)


class QCandParseSeqspecs:
    """Class method for checking seqspecs against portal metadata."""

    def __init__(self, analysis_set_metadata: portal_parsing.InputAnalysisSetMetadata, igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str):
        self.analysis_set_metadata = analysis_set_metadata
        self.analysis_set_acc = self.analysis_set_metadata.analysis_set_acc
        self.rna_input_info = self.analysis_set_metadata.rna_input_info
        self.atac_input_info = self.analysis_set_metadata.atac_input_info
        self.igvf_api = igvf_api
        self.partial_root_dir = partial_root_dir

    def _validate_final_barcode_file(self, seqspec_tool_output: seqspec_parsing.SeqSpecToolOutput, assay_type: str):
        """Check if any of the final inclusion lists created is empty."""
        with open(seqspec_tool_output.final_barcode_file, 'r') as file_obj:
            first_char = file_obj.read(1)
            if not first_char.isalpha():
                raise const.BadDataException(
                    f"The {assay_type} final inclusion lists are empty.")

    def _get_all_seqspec_metadata_per_assay_type(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata], assay_type: str) -> list[seqspec_parsing.SeqSpecMetadata]:
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
                    seqspec_file_path=curr_seqspec_file_path, igvf_api=self.igvf_api).generate_seqspec_metadata(seqfiles_metadata=measet_metadata.seqfiles, assay_type=assay_type)
                # Save it to the list
                all_seqspec_metadata.append(curr_seqspec_metadata)
        return all_seqspec_metadata

    def _get_single_seqspec_tool_output_per_assay_type(self, seqspec_metadata: seqspec_parsing.SeqSpecMetadata, onlist_method: str, assay_type: str) -> seqspec_parsing.SeqSpecToolOutput:
        """Get one seqspec tool output for a given assay type."""
        final_barcode_dir = os.path.join(
            self.partial_root_dir, 'final_barcode_onlist')
        if not os.path.exists(final_barcode_dir):
            os.makedirs(final_barcode_dir)
        curr_final_barcode_file_name = '_'.join(
            [self.analysis_set_acc, seqspec_metadata.modality, 'final_barcode_inclusion_list.txt'])
        curr_seqspec_tool_output = seqspec_parsing.GetSeqSpecToolOutput(
            seqspec_metadata=seqspec_metadata,
            onlist_method=onlist_method,
            output_barcode_list_file=os.path.join(final_barcode_dir, curr_final_barcode_file_name)).generate_seqspec_tool_output()
        # Validate the final barcode file is not empty
        self._validate_final_barcode_file(
            seqspec_tool_output=curr_seqspec_tool_output, assay_type=assay_type)
        return curr_seqspec_tool_output

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

    def check_if_read_index_match(self, all_seqspec_metadata: list[seqspec_parsing.SeqSpecMetadata], onlist_method: str, assay_type: str) -> bool:
        """Check if all seqspec read index strings match the first read index string."""
        mismatch_cnt = 0
        # First, get all read strings
        all_read_index_strings = []
        for seqspec_metadata in all_seqspec_metadata:
            curr_seqspec_tool_mthd = seqspec_parsing.GetSeqSpecToolOutput(
                seqspec_metadata=seqspec_metadata,
                onlist_method=onlist_method,
                output_barcode_list_file='place_holder.txt')
            all_read_index_strings.append(
                curr_seqspec_tool_mthd._get_seqspec_index())
        # Then, check if they all match the first one
        for read_string in all_read_index_strings[1:]:
            if read_string != all_read_index_strings[0]:
                mismatch_cnt += 1
        if mismatch_cnt > 0:
            raise const.BadDataException(
                f"{mismatch_cnt} {assay_type} read index string mismatches found.")
        return True

    def quality_check_one_assay_type(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata], assay_type: str) -> str | None:
        """Quality check the input analysis set metadata."""
        possible_errors = ''
        try:
            all_seqspec_metadata = self._get_all_seqspec_metadata_per_assay_type(
                measet_metadata_list=measet_metadata_list, assay_type=assay_type)
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
            # Check read index strings
            self.check_if_read_index_match(
                all_seqspec_metadata=all_seqspec_metadata, onlist_method=expected_onlist_method, assay_type=assay_type.upper())
            # If all good, use the first seqspec metadata to generate the seqspec tool output
            return self._get_single_seqspec_tool_output_per_assay_type(seqspec_metadata=all_seqspec_metadata[0], onlist_method=expected_onlist_method, assay_type=assay_type)
        except const.BadDataException as e:
            possible_errors += str(e)
        if possible_errors:
            # This syntax is to satisfy Terra input table format
            return seqspec_parsing.SeqSpecToolOutput(read_index_string="None", final_barcode_file="None", errors=possible_errors)


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
    atac_seqspec_urls: set[str] = dataclasses.field(default_factory=list)
    rna_seqspec_urls: set[str] = dataclasses.field(default_factory=list)
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


class GenerateTerraInputParams:
    """Class method for generating Terra input table."""

    def __init__(self, analysis_set_metadata: portal_parsing.InputAnalysisSetMetadata, igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str):
        # Analysis set metadata
        self.analysis_set_metadata = analysis_set_metadata
        self.analysis_set_acc = self.analysis_set_metadata.analysis_set_acc
        # Measurement sets metadata
        self.rna_input_info = self.analysis_set_metadata.rna_input_info
        self.atac_input_info = self.analysis_set_metadata.atac_input_info
        # Sample Info
        self.samples = self.analysis_set_metadata.sample_info
        # IGVF API client
        self.igvf_api = igvf_api
        # Partial root directory for saving intermediate files
        self.partial_root_dir = partial_root_dir

    def _check_and_get_seqspecs_for_input_params(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata], assay_type: str) -> seqspec_parsing.SeqSpecToolOutput:
        """Quality check the input analysis set metadata and get seqspec tool outputs."""
        seqspec_qc_generator_mthd = QCandParseSeqspecs(analysis_set_metadata=self.analysis_set_metadata,
                                                       igvf_api=self.igvf_api,
                                                       partial_root_dir=self.partial_root_dir)
        return seqspec_qc_generator_mthd.quality_check_one_assay_type(measet_metadata_list=measet_metadata_list, assay_type=assay_type)

    def _get_input_measet_accsions(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata]) -> list[str]:
        """Get all measurement set accessions for a given assay type."""
        return [measet_metadata.measet_acc for measet_metadata in measet_metadata_list]

    def _get_input_seqspec_urls(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata]) -> list[str]:
        """Get all seqspec URLs for a given assay type."""
        unique_seqspec_urls = set()
        for measet_metadata in measet_metadata_list:
            for seqfile_item in measet_metadata.seqfiles:
                if seqfile_item.seqspec_urls:
                    unique_seqspec_urls.update(seqfile_item.seqspec_urls)
        return sorted(unique_seqspec_urls)

    def _get_barcode_replacement_file(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata]) -> str:
        """Get the barcode replacement file for the analysis set."""
        unique_replacement_files = set()
        for measet_metadata in measet_metadata_list:
            if measet_metadata.barcode_replacement_file:
                unique_replacement_files.add(
                    measet_metadata.barcode_replacement_file)
        if len(unique_replacement_files) > 1:
            raise const.BadDataException(
                f"Multiple barcode replacement files found for analysis set {self.analysis_set_acc}.")
        elif len(unique_replacement_files) == 1:
            return list(unique_replacement_files)[0]
        else:
            # Terra takes this an optional File type. Because it is optional, return empty string if not present
            # Terra will read it in as NaN/Null
            return ""

    def _generate_fastq_urls_per_assay(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata], assay_type: str) -> dict:
        """Generate all FASTQ file info for a given assay type."""
        assay_type = assay_type.lower()
        # Get all sequence file metadata from all measurement sets of a specific assay type
        all_sequence_file_metadata = []
        for measet_metadata in measet_metadata_list:
            # for seqfile_metadata in measet_metadata.seqfiles:
            # curr_pipeline_read_type = seqfile_metadata.get_read_names_from_seqfile()
            all_sequence_file_metadata.extend(measet_metadata.seqfiles)
        # Sort them by file set, read type, sequencing run, lane, and flowcell ID
        sorted_all_sequence_file_metadata = sorted(
            all_sequence_file_metadata,
            key=lambda x: (
                x.file_set,
                x.illumina_read_type,
                x.sequencing_run,
                x.lane if x.lane is not None else 99,
                x.flowcell_id if x.flowcell_id is not None else 99
            )
        )

        # Output the FASTQ file URLs
        seqfile_urls_and_accs_by_reads = {f'{assay_type}_read1': [],
                                          f'{assay_type}_read2': [],
                                          f'{assay_type}_barcode': []
                                          }
        for seqfile_metadata in sorted_all_sequence_file_metadata:
            curr_pipeline_read_types = seqfile_metadata.get_read_names_from_seqfile(
                assay_type=assay_type)
            for curr_pipeline_read_type in curr_pipeline_read_types:
                curr_dict_key = '_'.join(
                    [assay_type, curr_pipeline_read_type])
                seqfile_urls_and_accs_by_reads.setdefault(curr_dict_key, []).append(
                    seqfile_metadata.file_url)
        # Add on accessions
        for read_type in list(seqfile_urls_and_accs_by_reads.keys()):
            urls = seqfile_urls_and_accs_by_reads[read_type]
            seqfile_urls_and_accs_by_reads[f'{read_type}_accessions'] = [
                get_seqspec_accession_from_path(url) for url in urls
            ]
        return seqfile_urls_and_accs_by_reads

    def _get_all_input_params_per_assay_type(self, measet_metadata_list: list[portal_parsing.MeasurementSetMetadata], assay_type: str) -> dict[str, list | str]:
        """Get all input params for a given assay type."""
        # Get all measurement set accessions
        measet_accessions = self._get_input_measet_accsions(
            measet_metadata_list=measet_metadata_list)
        # Get all seqspec URLs
        seqspec_urls = self._get_input_seqspec_urls(
            measet_metadata_list=measet_metadata_list)
        # Generate FASTQ file URLs and accessions
        seqfile_urls_and_accs_by_reads = self._generate_fastq_urls_per_assay(
            measet_metadata_list=measet_metadata_list, assay_type=assay_type)
        # Seqspec tool outputs
        seqspec_tool_outputs = self._check_and_get_seqspecs_for_input_params(
            measet_metadata_list=measet_metadata_list, assay_type=assay_type)
        # Build the return dict
        partial_params_dict = dict({
            f'{assay_type}_MeaSetIDs': measet_accessions,
            f'{assay_type}_read1_accessions': seqfile_urls_and_accs_by_reads.get(f'{assay_type}_read1_accessions', []),
            f'{assay_type}_read2_accessions': seqfile_urls_and_accs_by_reads.get(f'{assay_type}_read2_accessions', []),
            f'{assay_type}_barcode_accessions': seqfile_urls_and_accs_by_reads.get(f'{assay_type}_barcode_accessions', []),
            f'{assay_type}_seqspec_urls': seqspec_urls,
            f'{assay_type}_read1': seqfile_urls_and_accs_by_reads.get(f'{assay_type}_read1', []),
            f'{assay_type}_read2': seqfile_urls_and_accs_by_reads.get(f'{assay_type}_read2', []),
            f'{assay_type}_barcode': seqfile_urls_and_accs_by_reads.get(f'{assay_type}_barcode', []),
            f'{assay_type}_barcode_inclusion_list': seqspec_tool_outputs.final_barcode_file,
            f'{assay_type}_read_format': seqspec_tool_outputs.read_index_string
        })
        if assay_type == 'rna':
            # Barcode replacement file if RNAseq
            partial_params_dict['barcode_replacement_file'] = self._get_barcode_replacement_file(
                measet_metadata_list=measet_metadata_list)
        if seqspec_tool_outputs.errors:
            partial_params_dict['possible_errors'] = seqspec_tool_outputs.errors
        return partial_params_dict

    def _get_taxa(self) -> str:
        """Get the taxa for the analysis set."""
        unique_taxa = set()
        for sample in self.samples:
            if sample.taxa:
                unique_taxa.add(sample.taxa)
        if len(unique_taxa) > 1:
            raise const.BadDataException(
                f"Multiple taxa found for analysis set {self.analysis_set_acc}.")
        elif len(unique_taxa) == 1:
            if unique_taxa != {'Mixed species'}:
                return list(unique_taxa)[0]
            else:
                raise const.BadDataException(
                    f"Mixed species taxa found for analysis set {self.analysis_set_acc}.")
        else:
            raise const.BadDataException(
                f"No taxa found for analysis set {self.analysis_set_acc}.")

    def _get_reference_files(self, taxa: str) -> const.RunReferenceFiles:
        """Get the reference files for the analysis set."""
        if taxa not in const.TAXA_TO_GENOME_REF_FILES:
            raise const.BadDataException(
                f"Taxa '{taxa}' not recognized for analysis set {self.analysis_set_acc}.")
        return const.TAXA_TO_GENOME_REF_FILES[taxa]

    def _get_subpool_id(self) -> str:
        """Get the subpool ID for the analysis set."""
        unique_subpool_ids = set()
        for sample in self.samples:
            if sample.subpool_id:
                unique_subpool_ids.add(sample.subpool_id)
        if len(unique_subpool_ids) >= 1:
            return '-'.join(sorted(unique_subpool_ids))
        else:
            raise const.BadDataException(
                f"No subpool ID found for analysis set {self.analysis_set_acc}.")

    def _get_onlist_mapping_bool(self) -> bool:
        """Determine if onlist mapping is needed for the analysis set."""
        # If any of the measurement sets has combinatorial onlist method, then we need onlist mapping
        if any(assay_title in const.ASSAYS_NEED_ONLIST_MAPPING for assay_title in self.analysis_set_metadata.preferred_assay_titles):
            return True
        return False

    def generate_terra_input_params(self) -> TerraPipelineParams:
        possible_errors = []
        params = {}
        params['kb_strand'] = 'forward'
        params['analysis_set_acc'] = self.analysis_set_acc

        try:
            params['onlist_mapping'] = self._get_onlist_mapping_bool()
        except const.BadDataException as e:
            possible_errors.append(f"onlist_mapping: {e}")

        try:
            params['taxa'] = self._get_taxa()
        except const.BadDataException as e:
            possible_errors.append(f"taxa: {e}")
            params['taxa'] = ''

        try:
            params['subpool_id'] = self._get_subpool_id()
        except const.BadDataException as e:
            possible_errors.append(f"subpool_id: {e}")
            params['subpool_id'] = ''

        try:
            reference_files = self._get_reference_files(
                taxa=params.get('taxa', ''))
            params.update(dataclasses.asdict(reference_files))
        except const.BadDataException as e:
            possible_errors.append(f"reference_files: {e}")

        if self.analysis_set_metadata.rna_input_info:
            try:
                rna_partial = self._get_all_input_params_per_assay_type(
                    measet_metadata_list=self.analysis_set_metadata.rna_input_info, assay_type='rna')
                params.update(rna_partial)
                if 'possible_errors' in rna_partial and rna_partial['possible_errors']:
                    possible_errors.append(
                        f"rna: {rna_partial['possible_errors']}")
            except const.BadDataException as e:
                possible_errors.append(f"rna: {e}")

        if self.analysis_set_metadata.atac_input_info:
            try:
                atac_partial = self._get_all_input_params_per_assay_type(
                    measet_metadata_list=self.analysis_set_metadata.atac_input_info, assay_type='atac')
                params.update(atac_partial)
                if 'possible_errors' in atac_partial and atac_partial['possible_errors']:
                    possible_errors.append(
                        f"atac: {atac_partial['possible_errors']}")
            except const.BadDataException as e:
                possible_errors.append(f"atac: {e}")

        # Fill in any missing fields with their default values
        for field in dataclasses.fields(TerraPipelineParams):
            if field.name not in params:
                params[field.name] = field.default if field.default is not dataclasses.MISSING else field.default_factory(
                ) if field.default_factory is not dataclasses.MISSING else None

        params['possible_errors'] = "\n".join(
            possible_errors) if possible_errors else "None"
        terra_input_params = TerraPipelineParams(**params)
        return terra_input_params


class ConvertParamsToTerraTable:
    """Class method for converting Terra input params to a Terra input table."""

    def __init__(self, all_terra_input_params: list[TerraPipelineParams], terra_etype: str, local_barcode_file_dir: str, gs_barcode_list_bucket: str):
        self.all_terra_input_params = all_terra_input_params
        self.terra_etype = terra_etype
        self.local_barcode_file_dir = local_barcode_file_dir
        self.gs_barcode_list_bucket = gs_barcode_list_bucket

    def _terra_str_formatter(self, input_strs: list) -> str:
        """Convert Python default array of string output to Terra format (in [] with double quotes)."""
        return '[' + ','.join([f'"{v}"' for v in input_strs]) + ']'

    def _reformat_arrays_to_terra_format(self):
        """Reformat all read file URLs arrays to Terra format."""
        array_attrs = [
            'atac_read1', 'atac_read2', 'atac_barcode',
            'rna_read1', 'rna_read2', 'rna_barcode'
        ]
        # Barcode replacement file is an optional value, it must be Null/NaN in Terra if not present
        string_attrs = [
            'atac_read_format', 'rna_read_format',
            'atac_barcode_inclusion_list', 'rna_barcode_inclusion_list'
        ]
        all_terra_input_params_formatted = []
        for terra_input_params in self.all_terra_input_params:
            params_dict = dataclasses.asdict(terra_input_params)
            # Format arrays
            for attr in array_attrs:
                val = params_dict[attr]
                if isinstance(val, list) and val:
                    params_dict[attr] = self._terra_str_formatter(val)
                elif isinstance(val, list) and not val:
                    params_dict[attr] = "[]"
            # Format strings
            for attr in string_attrs:
                val = params_dict[attr]
                if isinstance(val, str) and val == '':
                    params_dict[attr] = "None"
            all_terra_input_params_formatted.append(
                TerraPipelineParams(**params_dict))
        return all_terra_input_params_formatted

    def _mod_input_table_for_terra(self, pipeline_input_table: pd.DataFrame) -> pd.DataFrame:
        """Modify the pipeline input table to have Terra table ID in column name and replace local file path with gs:// bucket path."""
        pipeline_input_table.index = pipeline_input_table['analysis_set_acc']
        pipeline_input_table.index.name = f'entity:{self.terra_etype}_id'
        pipeline_input_table_with_gspaths = pipeline_input_table.replace(
            self.local_barcode_file_dir.rstrip('/'), self.gs_barcode_list_bucket.rstrip('/'), regex=True)
        return pipeline_input_table_with_gspaths

    def generate_terra_input_table(self) -> pd.DataFrame:
        """Generate the Terra input table as a pandas DataFrame."""
        formatted_terra_input_params = self._reformat_arrays_to_terra_format()
        params_dicts = [dataclasses.asdict(params)
                        for params in formatted_terra_input_params]
        terra_input_df = pd.DataFrame(params_dicts)
        return self._mod_input_table_for_terra(pipeline_input_table=terra_input_df)


class CompleteTerraForming:
    """Class method for generating and saving the complete Terra input table."""

    def __init__(self, analysis_set_accessions: list[str], igvf_api: igvf_client.api.igvf_api.IgvfApi, partial_root_dir: str, terra_etype: str, local_barcode_file_dir: str, gs_barcode_list_bucket: str):
        self.analysis_set_accessions = analysis_set_accessions
        self.igvf_api = igvf_api
        self.partial_root_dir = partial_root_dir
        self.terra_etype = terra_etype
        self.local_barcode_file_dir = local_barcode_file_dir
        self.gs_barcode_list_bucket = gs_barcode_list_bucket

    def generate_complete_terra_input_table(self) -> pd.DataFrame:
        """Generate the complete Terra input table for all analysis sets."""
        all_terra_input_params = []
        for analysis_set_acc in self.analysis_set_accessions:
            logging.info(
                f'>>>>>>>>> Processing analysis set {analysis_set_acc}...')
            # Get the analysis set metadata
            analysis_set_metadata = portal_parsing.GetAnalysisSetMetadata(
                igvf_api=self.igvf_api, analysis_set_accession=analysis_set_acc).get_input_analysis_set_metadata()
            # Generate Terra input params for the analysis set
            terra_input_params = GenerateTerraInputParams(analysis_set_metadata=analysis_set_metadata,
                                                          igvf_api=self.igvf_api,
                                                          partial_root_dir=self.partial_root_dir).generate_terra_input_params()
            all_terra_input_params.append(terra_input_params)
        # Convert all Terra input params to a Terra input table
        terra_input_table = ConvertParamsToTerraTable(all_terra_input_params=all_terra_input_params,
                                                      terra_etype=self.terra_etype,
                                                      local_barcode_file_dir=self.local_barcode_file_dir,
                                                      gs_barcode_list_bucket=self.gs_barcode_list_bucket).generate_terra_input_table()
        return terra_input_table


def save_pipeline_input_table(terra_data_table: pd.DataFrame, partial_root_dir: str) -> pd.DataFrame:
    """Save the pipeline input table to a TSV file.

    Args:
        pipeline_input_table (pd.DataFrame): The pipeline input table
        output_dir (str): The output directory
    """
    curr_datetime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    if os.path.exists(partial_root_dir) is False:
        os.makedirs(partial_root_dir)
    # Save the pipeline input table to a TSV file
    terra_data_table.to_csv(os.path.join(
        partial_root_dir, f'single-cell_uniform_pipeline_input_table_{curr_datetime}_withGSpath.tsv'), sep='\t')
