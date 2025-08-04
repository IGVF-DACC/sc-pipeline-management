"""
Utility functions for Terra to Portal posting operations.
Functions for parsing, string generation, file handling, and other utilities
that don't directly call .post or .patch operations.
"""

import pandas as pd
import dataclasses
import re
import firecloud.api as fapi
import io

fapi._set_session()

# IGVF file url parsing regex for accession
IGVF_URL_PATH_REGEX = re.compile(
    r'\/.*-files\/(IGVF[A-Z0-9]+|TST[A-Z0-9]+)\/@@download')

# GS file path regex
GS_FILE_PATH_REGEX = re.compile(
    r'gs://([a-z0-9\-]+)/submissions/final-outputs/([a-z0-9\-]+)/single_cell_pipeline/([a-z0-9\-]+)/[a-z\-]+/[a-z\_]+/([a-z0-9\-]+)/.*')


@dataclasses.dataclass(frozen=True)
class TerraJobUUIDs:
    """Data class to hold pipeline output UUIDs."""
    gcloud_bucket: str
    submission_id: str
    workflow_id: str
    subworkflow_id: str

    def aliases(self) -> str:
        """Generate a string of the form submissionID_workflowID_subworkflowID as file."""
        return f'{self.gcloud_bucket}_{self.submission_id}_{self.workflow_id}_{self.subworkflow_id}'

    def input_param_aliases(self) -> str:
        """Generate a string of the form submissionID_workflowID."""
        return f'{self.submission_id}_{self.workflow_id}'


@dataclasses.dataclass(frozen=True)
class InputFileAccs:
    """Data class to hold derived from files."""
    sequence_files: list[str]
    seqspec_files: list[str]
    rna_barcode_replacement: list[str]
    reference_files: list[str]


@dataclasses.dataclass(frozen=True)
class InputFileInfo:
    sequence_file_headers: list[str]
    seqspec_file_headers: str
    reference_file_headers: list[str]
    barcode_replacement_file_header: str | None


@dataclasses.dataclass(frozen=True)
class TerraPipelineDataTable:
    terra_datatable: pd.DataFrame


# Reference file and derived from file info by assay type
INPUT_FILE_HEADERS_BY_ASSAY_TYPE = {
    'atac': InputFileInfo(
        sequence_file_headers=['atac_read1_accessions',
                               'atac_read2_accessions', 'atac_barcode_accessions'],
        seqspec_file_header='atac_seqspec_urls',
        reference_file_headers=['chromap_index', 'genome_fasta'],
        barcode_replacement_file_header=None
    ),
    'rna': InputFileInfo(
        sequence_file_headers=['rna_read1_accessions',
                               'rna_read2_accessions', 'rna_barcode_accessions'],
        seqspec_file_header='rna_seqspec_urls',
        reference_file_headers=['kb_index'],
        barcode_replacement_file_header='barcode_replacement_file'
    )
}


def _parse_terra_str_list(terra_str_lists: list[str]) -> list[str]:
    """Parse a string list from Terra data table to a Python list."""
    parse_list_of_strs = []
    for value in terra_str_lists:
        if str(value) == 'nan':
            continue
        if str(value) == '[]':
            continue
        if str(value) == '[""]':
            continue
        if str(value).startswith("[") and str(value).endswith("]"):
            parse_list_of_strs.extend([ele.strip("'")
                                      for ele in value[1:-1].split(', ')])
        else:
            parse_list_of_strs.append(value)
    return parse_list_of_strs


def _get_multiome_types(terra_data_record: pd.Series) -> list:
    """Check based on input measurement set IDs to see if the pipeline is RNAseq-only, ATACseq-only, or multiome-seq."""
    multiome_map = {'atac_MeaSetIDs': 'ATACseq', 'rna_MeaSetIDs': 'RNAseq'}
    pipeline_output_content = []
    for key, value in multiome_map.items():
        if len(_parse_terra_str_list(terra_str_lists=[terra_data_record[key]])) > 0:
            pipeline_output_content.append(value)
    return pipeline_output_content


class TerraOutputMetadata:
    def __init__(self, terra_data_record: pd.Series, igvf_client_api):
        self.terra_data_record = terra_data_record
        self.igvf_api = igvf_client_api
        self.anaset_accession = self.terra_data_record['analysis_set_acc']
        self.taxa = self.terra_data_record['taxa']
        self.multiome_type = _get_multiome_types(
            terra_data_record=self.terra_data_record)

    def _get_gs_path_for_terra_output_cols(self) -> str:
        """Get GS path from terra data record for submission and workflow IDs parsing, trying RNA first, then ATAC."""
        # Possible columns to check for GS path (RNA first, then ATAC)
        possible_columns = ['rna_kb_h5ad', 'atac_bam']
        for col in possible_columns:
            # If the column is empty, it returns np.float64(nan)
            if (col in list(self.terra_data_record.index)) and (pd.notna(self.terra_data_record[col])):
                return self.terra_data_record[col]
        raise ValueError(
            'No valid GS path found in the Terra data record for workflow UUID parsing.')

    def _parse_workflow_uuids_from_gs_path(self, gs_path_regex: re.Pattern = GS_FILE_PATH_REGEX) -> TerraJobUUIDs:
        """Parse workflow UUID from a gs path."""
        usable_gs_path = self._get_gs_path_for_terra_output_cols()
        matches = re.search(gs_path_regex, usable_gs_path)
        if not matches:
            raise ValueError(
                'Unable to parse workflow UUIDs from file GCP path.')
        uuids = matches.groups()
        return TerraJobUUIDs(gcloud_bucket=uuids[0],
                             submission_id=uuids[1],
                             workflow_id=uuids[2],
                             subworkflow_id=uuids[3])

    def _get_input_file_accs_from_table(self, assay_type: str) -> list:
        """Get a list of sequence file accessions (for derived_from) from Terra table."""
        input_file_info = INPUT_FILE_HEADERS_BY_ASSAY_TYPE.get(assay_type)
        # Build the derived from file info
        # Get sequence file accessions from the Terra table
        seqfile_accessions = self._parse_terra_str_list(
            terra_str_lists=[self.terra_data_record[seqfile_col] for seqfile_col in input_file_info.sequence_file_headers])
        # Get sequence specification file accessions from the file URLs in the Terra table
        seqspec_accessions = self._parse_igvf_accessions_from_urls(
            igvf_file_urls=[self.terra_data_record[input_file_info.seqspec_file_header]])
        # Get RNA barcode replacement accessions from the file URLs in the Terra table (can be empty for non-Parse)
        barcode_replacement_accession = self._parse_igvf_accessions_from_urls(
            igvf_file_urls=[self.terra_data_record.get(input_file_info.barcode_replacement_file_header)])
        # Get reference files from the file URLs in the Terra table
        reference_files = self._parse_igvf_accessions_from_urls(
            igvf_file_urls=[self.terra_data_record[ref_file_header] for ref_file_header in input_file_info.reference_file_headers])
        # Build the InputFileAccs dataclass
        return InputFileAccs(sequence_files=seqfile_accessions,
                             seqspec_files=seqspec_accessions,
                             rna_barcode_replacement=barcode_replacement_accession,
                             reference_files=reference_files)
