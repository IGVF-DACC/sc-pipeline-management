import urllib.parse
import urllib
import dataclasses
import igvf_client


import sc_pipe_management.input_params_generation.constant as const


def construct_full_href_url(igvf_href: str) -> str:
    """_summary_

    Args:
        igvf_href (str): _description_

    Returns:
        str: _description_
    """
    return urllib.parse.urljoin(base=const.BASE_IGVF_PORTAL_URL, url=igvf_href)


@dataclasses.dataclass(frozen=True)
class SeqFileMetadata:
    """Dataclass to hold sequence file metadata."""
    file_accession: str
    file_set: str
    illumina_read_type: str
    sequencing_run: int
    lane: int | None
    flowcell_id: str | None
    file_url: str
    seqspec_urls: list[str]
    read_names: list[str]

    def get_read_names_from_seqfile(self, assay_type: str) -> list[str]:
        """Get the read names from the seqfile object. For ATACseq, if Read2 and Barcode index are concatenated, return both.
        For RNAseq, if Read2 and Barcode index are concatenated, return Read2 only."""
        if len(self.read_names) == 1:
            if self.read_names[0] != 'Barcode index':
                # return read1 or read2
                return [self.read_names[0].lower().replace(' ', '')]
            else:
                return ['barcode']
        elif self.read_names == ['Read 2', 'Barcode index']:
            # If Read2 and Barcode index are concatenated, RNAseq only needs read2
            if assay_type == 'rna':
                return ['read2']
            # If ATACseq, return the same fastq file for both read2 and barcode
            elif assay_type == 'atac':
                return ['read2', 'barcode']
        else:
            raise const.BadDataException(
                f'Malformed seqfile read names: {self.read_names}. Expected only one read name or Read 2 and Barcode index.')


@dataclasses.dataclass(frozen=True)
class MeasurementSetMetadata:
    """Dataclass to hold measurement set metadata."""
    measet_acc: str
    assay_type: str
    onlist_mapping: bool
    seqfiles: list[SeqFileMetadata]
    onlist_method: str
    onlist_files: list[str]
    barcode_replacement_file: str | None


@dataclasses.dataclass(frozen=True)
class SampleMetadata:
    """Dataclass to hold sample metadata."""
    taxa: str
    subpool_id: str


@dataclasses.dataclass(frozen=True)
class InputAnalysisSetMetadata:
    analysis_set_acc: str
    # Input file sets
    rna_input_info: list[MeasurementSetMetadata] | None
    atac_input_info: list[MeasurementSetMetadata] | None
    # Sample metadata
    sample_info: list[SampleMetadata]
    # Assay titles
    preferred_assay_titles: list[str]


class GetSeqFileMetadata:
    """Class to get sequence file metadata from a list of sequence file objects."""

    def __init__(self, seqfile_id: str, igvf_api: igvf_client.api.igvf_api.IgvfApi):
        """Initialize with a list of sequence file objects."""
        self.seqfile_id = seqfile_id
        self.igvf_api = igvf_api
        self.seqfile_obj = self.igvf_api.get_by_id(seqfile_id).actual_instance

    def _get_seqspec_urls(self) -> list[str]:
        """Get the seqspec URLs from the sequence file object."""
        if not self.seqfile_obj.seqspecs:
            return []
        if self.seqfile_obj.seqspecs:
            seqspec_urls = set()
            for seqspec_id in self.seqfile_obj.seqspecs:
                curr_seqspec_item = self.igvf_api.get_by_id(
                    seqspec_id).actual_instance
                if curr_seqspec_item.status not in const.FILE_DEPRECATED_STATUSES:
                    seqspec_urls.add(construct_full_href_url(
                        igvf_href=curr_seqspec_item.href))
        return sorted(seqspec_urls)

    def get_seqfile_metadata(self) -> SeqFileMetadata:
        """Get the sequence file metadata."""
        return SeqFileMetadata(
            file_accession=self.seqfile_obj.accession,
            file_set=self.seqfile_obj.file_set,
            illumina_read_type=self.seqfile_obj.illumina_read_type,
            sequencing_run=self.seqfile_obj.sequencing_run,
            lane=self.seqfile_obj.lane,
            flowcell_id=self.seqfile_obj.flowcell_id,
            file_url=construct_full_href_url(igvf_href=self.seqfile_obj.href),
            seqspec_urls=self._get_seqspec_urls(),
            read_names=self.seqfile_obj.read_names
        )


class GetMeasurementSetMetadata:
    """Class to get measurement set metadata from a list of measurement set objects."""

    def __init__(self, measet_id: str, igvf_api: igvf_client.api.igvf_api.IgvfApi):
        """Initialize with a measurement set ID and IGVF API client."""
        self.measet_id = measet_id
        self.igvf_api = igvf_api
        # Get the measurement set object
        self.measet_obj = self.igvf_api.get_by_id(
            measet_id).actual_instance

    def _check_needs_onlist_mapping(self) -> bool:
        """Get the onlist mapping status of the measurement set items."""
        preferred_assay_titles = self.measet_obj.preferred_assay_titles
        return any(assay_title in const.ASSAYS_NEED_ONLIST_MAPPING for assay_title in preferred_assay_titles)

    def _get_seqfile_metadata(self) -> list[SeqFileMetadata]:
        """Get sequence file metadata from the measurement set object."""
        seqfile_metadata_list = []
        for seqfile_id in self.measet_obj.files:
            if seqfile_id.startswith('/sequence-files/'):
                curr_seqfile_metadata = GetSeqFileMetadata(
                    seqfile_id=seqfile_id, igvf_api=self.igvf_api).get_seqfile_metadata()
                if curr_seqfile_metadata.read_names is None:
                    continue
                seqfile_metadata_list.append(curr_seqfile_metadata)
        return seqfile_metadata_list

    def _get_barcode_replacement_file_url(self) -> str | None:
        """Get the barcode replacement file URL if exists."""
        if self.measet_obj.barcode_replacement_file:
            brf_obj = self.igvf_api.get_by_id(
                self.measet_obj.barcode_replacement_file).actual_instance
            return construct_full_href_url(
                igvf_href=brf_obj.href)
        return None

    def get_measurement_set_metadata(self) -> MeasurementSetMetadata:
        """Get the measurement set metadata."""
        return MeasurementSetMetadata(
            measet_acc=self.measet_obj.accession,
            assay_type=const.ASSAY_NAMES_CONVERSION_REF[self.measet_obj.assay_term],
            onlist_mapping=self._check_needs_onlist_mapping(),
            seqfiles=self._get_seqfile_metadata(),
            # The following fields are checked by portal audits
            onlist_method=self.measet_obj.onlist_method,
            onlist_files=self.measet_obj.onlist_files,
            barcode_replacement_file=self._get_barcode_replacement_file_url()
        )


class GetAnalysisSetMetadata:
    def __init__(self, analysis_set_accession: str, igvf_api: igvf_client.api.igvf_api.IgvfApi):
        """Initialize with analysis set accession and IGVF API client."""
        self.igvf_api = igvf_api
        # Get the analysis set accession and the analysis set object
        self.analysis_set_accession = analysis_set_accession
        self.analysis_set_obj = self.igvf_api.get_by_id(
            f'/analysis-sets/{analysis_set_accession}').actual_instance
        # Get the input measurement set metadata

    def _get_input_samples_info(self):
        """Get taxa and subpool ID from the analysis set samples."""
        sample_list = []
        for sample_id in self.analysis_set_obj.samples:
            # Get the sample object
            sample_obj = self.igvf_api.get_by_id(sample_id).actual_instance
            # Get taxa and subpool ID
            sample_list.append(SampleMetadata(
                taxa=sample_obj.taxa,
                subpool_id=sample_obj.accession))
        # Output the sample metadata
        return sample_list

    def _get_input_measet_metadata_by_assay_type(self) -> dict[str, MeasurementSetMetadata]:
        """Get input measurement set metadata by assay type."""
        measet_metadata_by_assay_type = {'rna': [], 'atac': []}
        for fileset_id in self.analysis_set_obj.input_file_sets:
            if not fileset_id.startswith('/measurement-sets/'):
                continue
            # Get the measurement set object
            measet_metadata_mthd = GetMeasurementSetMetadata(
                measet_id=fileset_id, igvf_api=self.igvf_api)
            curr_measet_metadata = measet_metadata_mthd.get_measurement_set_metadata()
            curr_assay_type = curr_measet_metadata.assay_type
            measet_metadata_by_assay_type[curr_assay_type].append(
                curr_measet_metadata)
        # If no measurement sets are found for a specific assay type, return an empty list
        return measet_metadata_by_assay_type

    def get_input_analysis_set_metadata(self) -> InputAnalysisSetMetadata:
        """Get the input analysis set metadata."""
        measet_metadata_by_assay_type = self._get_input_measet_metadata_by_assay_type()
        return InputAnalysisSetMetadata(
            analysis_set_acc=self.analysis_set_accession,
            rna_input_info=measet_metadata_by_assay_type['rna'],
            atac_input_info=measet_metadata_by_assay_type['atac'],
            sample_info=self._get_input_samples_info(),
            preferred_assay_titles=self.analysis_set_obj.preferred_assay_titles
        )
