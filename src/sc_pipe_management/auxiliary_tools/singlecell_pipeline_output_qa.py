from igvf_client import Configuration
from igvf_client import ApiClient
from igvf_client import IgvfApi
import os
import dataclasses
import re


@dataclasses.dataclass(frozen=True)
class QCMetricsMetadata:
    """Dataclass to hold QC metrics metadata."""
    qc_metrics_type: str
    qc_metrics_desc: str
    qc_metrics_attachments: list[str] | None


@dataclasses.dataclass(frozen=True)
class FileMetadata:
    """Dataclass to hold file metadata."""
    file_object:
    file_type: str
    content_type: str
    file_format: str
    description: str
    file_format_type: str | None
    analysis_step_version: str
    file_format_specifications: list[str] | None


@dataclasses.dataclass(frozen=True)
class PipelineDocument:
    description: str
    document_type: str


@dataclasses.dataclass(frozen=True)
class AnalysisSetMetadata:
    """Dataclass to hold analysis set metadata."""
    accession: str
    has_rnaseq: bool
    has_atacseq: bool
    documents: list[str]
    files: {str, list[FileMetadata]}
    uniform_pipeline_status: str | None


@dataclasses.dataclass(frozen=True)
class PerObjectQAResults:
    """Dataclass to hold per file QA results."""
    file_metadata_check: bool | str | None  # True if pass, error message if fail
    # True if pass, error message if fail, None if no QC expected
    qc_metrics_check: bool | str | None
    # True if pass, error message if fail
    submitted_file_name_check: bool | str | None
    document_info_check: bool | str | None  # True if pass, error message if fail


@dataclasses.dataclass(frozen=True)
class QACheckResults:
    """Dataclass to hold QA check results."""
    analysis_set_accession: str
    h5ad_check_results: PerObjectQAResults
    tar_check_results: PerObjectQAResults
    bam_check_results: PerObjectQAResults
    bed_check_results: PerObjectQAResults
    documents_check_results: PerObjectQAResults


FILE_TYPE_REGEX = re.compile(r'\/([a-z]*)\-files\/(IGVFFI[0-9]+)')

DEPRECATED_STATUSES = ['deleted', 'revoke', 'archived', 'replaced']

PROPER_DOCUMENT_INFO = PipelineDocument(description='Single-cell RNAseq and ATACseq pipeline output',
                                        document_type='pipeline parameters')

PROPER_QC_METRICS = {
    'bam': QCMetricsMetadata(
        qc_metrics_type='single-cell-atac-seq-quality-metrics',
        qc_metrics_desc='ATACseq chromap alignment QC metric',
        qc_metrics_attachments=['atac_bam_summary_stats']
    ),
    'bed': QCMetricsMetadata(
        qc_metrics_type='single-cell-atac-seq-quality-metrics',
        qc_metrics_desc='ATACseq chromap fragments QC metric',
        qc_metrics_attachments=[
            'atac_fragment_summary_stats', 'atac_fragments_alignment_stats']
    ),
    'gene_counts': QCMetricsMetadata(
        qc_metrics_type='single-cell-rna-seq-quality-metrics',
        qc_metrics_desc='RNAseq Kallisto Bustools QC metric',
        qc_metrics_attachments=['rnaseq_kb_info']
    ),
}

PROPER_ATAC_OUTPUT_FILES = {
    'bam': {
        'bam': FileMetadata(
            file_type='alignment',
            content_type='alignments',
            file_format='bam',
            file_format_type=None,
            description='Raw Aligned bam file from Chromap',
            analysis_step_version='igvf:single-cell-uniform-pipeline-chromap-atacseq-step-v1.1.0'
        ),
        'bai': FileMetadata(
            file_type='index',
            content_type='index',
            file_format='bai',
            file_format_type=None,
            description='Raw Aligned bam file index from Chromap',
            analysis_step_version='igvf:single-cell-uniform-pipeline-chromap-atacseq-step-v1.1.0'
        )
    },
    'bed': {
        'bed': FileMetadata(
            file_type='tabular',
            content_type='fragments',
            file_format='bed',
            file_format_type='bed3+',
            description='Raw Fragment file from Chromap',
            analysis_step_version='igvf:single-cell-uniform-pipeline-chromap-atacseq-step-v1.1.0',
            file_format_specifications=[
                'buenrostro-bernstein:igvf-sc-pipeline-atac-fragments-hdf5-specification'
            ]
        ),
        'tbi': FileMetadata(
            file_type='index',
            content_type='index',
            file_format='tbi',
            file_format_type=None,
            description='ATACseq chromap fragments QC metric',
            analysis_step_version='igvf:single-cell-uniform-pipeline-chromap-atacseq-step-v1.1.0',
        )
    }
}

PROPER_RNA_MATRIX_FILES = {
    'tar': FileMetadata(
        file_type='matrix',
        content_type='kallisto single cell RNAseq output',
        file_format='tar',
        file_format_type=None,
        description='Raw Tarball containing all the matrices, logs, and bus files generated from kb',
        analysis_step_version='igvf:single-cell-uniform-pipeline-kallisto-rnaseq-step-v1.1.0',
        file_format_specifications=sorted([
            'buenrostro-bernstein:igvf-sc-pipeline-matrix-tar-specification',
            'igvf:igvf-sc-pipeline-rna-tar-mtx-per-file-specification'
        ])
    ),
    'h5ad': FileMetadata(
        file_type='matrix',
        content_type='sparse gene count matrix',
        file_format='h5ad',
        file_format_type=None,
        description='Raw h5ad containing four separated count matrices: Spliced, Unspliced, Ambiguous, and Total',
        analysis_step_version='igvf:single-cell-uniform-pipeline-kallisto-rnaseq-step-v1.1.0',
        file_format_specifications=[
            'buenrostro-bernstein:igvf-sc-pipeline-matrix-h5-specification'
        ]
    )
}


def get_igvf_client_auth(igvf_api_keys: dict, igvf_endpoint: str = 'sandbox'):
    """Set up IGVF data portal access and set up IGVF python client api.

    Args:
        igvf_api_keys (dict): A dictionary containing the public and secret API keys.
        igvf_endpoint (str, optional): The IGVF site to use. Defaults to 'sandbox'.

    Returns:
        IgvfApi: An instance of the IgvfApi client.
    """
    config = Configuration(
        access_key=os.environ.get('IGVF_API_KEY'),
        secret_access_key=os.environ.get('IGVF_SECRET_KEY'),
        host='prod',
    )
    client = ApiClient(config)
    return IgvfApi(client)


class ParseAnalysisSetContent:
    """Class to handle the QA of single-cell pipeline outputs."""

    def __init__(self, igvf_client_api: IgvfApi, analysis_set_accession: str):
        self.igvf_client_api = igvf_client_api
        self.analysis_set_accession = analysis_set_accession
        self.analysis_object = self.igvf_client_api.get_by_id(
            f'/analysis_sets/{analysis_set_accession}/').actual_instance

    def _check_assay_types(self, assay_type: str) -> bool:
        """Check the assay types in the analysis set."""
        return any(assay for assay in self.analysis_object.assay_titles if assay_type in assay)

    def _get_main_file_metadata(self, file_id: str) -> FileMetadata:
        """Get file metadata from the provided data, specifically on the main data files (bam, h5ad, and tar)."""
        if file_id.startswith('/index-files/'):
            return None
        file_obj = self.igvf_client_api.get_by_id(file_id).actual_instance
        if file_obj.status in DEPRECATED_STATUSES:
            return None
        return FileMetadata(
            file_type=FILE_TYPE_REGEX.search(file_obj.type).group(0),
            content_type=file_obj.content_type,
            file_format=file_obj.file_format,
            file_format_type=file_obj.file_format_type,
            description=file_obj.description,  # Add this
            analysis_step_version=file_obj.analysis_step_version,  # Add this
            file_format_specifications=sorted(getattr(
                file_obj, 'file_format_specifications', None))  # Add this
        )

    def get_anaset_metadata(self) -> AnalysisSetMetadata:
        """Get analysis set metadata from the provided data."""
        analysis_set_main_files = {}
        for file_id in self.analysis_object.files:
            curr_file_metadata = self._get_main_file_metadata(file_id)
            # Skip files that are not main files or have been deprecated
            if curr_file_metadata is None:
                continue
            analysis_set_main_files.setdefault(
                curr_file_metadata.file_format, []).append(file_id)

        return AnalysisSetMetadata(
            accession=self.analysis_set_accession,
            has_rnaseq=self._check_assay_types('RNA'),
            has_atacseq=self._check_assay_types('ATAC'),
            documents=self.analysis_object.documents,
            files=analysis_set_main_files,
            uniform_pipeline_status=self.analysis_object.uniform_pipeline_status
        )


class SingleCellPipelineOutputQA:
    """Class to handle the QA of single-cell pipeline outputs."""

    def __init__(self, igvf_client_api: IgvfApi, analysis_set_metadata: AnalysisSetMetadata):
        self.igvf_client_api = igvf_client_api
        self.analysis_set_metadata = analysis_set_metadata
        self.analysis_object = self.analysis_set_metadata.analysis_object

    def _has_correct_submitted_file_name(self, file_id: str) -> bool:
        """Check if the submitted file name matches the expected format."""
        file_obj = self.igvf_client_api.get_by_id(file_id).actual_instance
        filename = file_obj.submitted_file_name.split('/')[-1]
        if not filename.startswith(self.analysis_set_metadata.accession):
            raise ValueError(
                f"File {file_obj.accession} has incorrect submitted_file_name."
            )
        return True

    def _has_correct_file_metadata(self, file_metadata: FileMetadata, qa_standards: dict) -> bool:
        """Check if the file format matches the expected format."""
        # Check file format
        if file_metadata.file_format not in qa_standards:
            raise ValueError(
                f"File {file_metadata.content_type} is not a valid format for the analysis set."
            )
        # Check file metadata
        expected_metadata = qa_standards[file_metadata.file_format]
        for field in expected_metadata.__dataclass_fields__:
            actual_value = getattr(file_metadata, field)
            expected_value = getattr(expected_metadata, field)

            if actual_value != expected_value:
                raise ValueError(
                    f"File {field} mismatch: expected '{expected_value}', got '{actual_value}'"
                )

        return True

    def _has_correct_qc_metrics(self, file_id: str, file_metadata: FileMetadata, qc_metadata: QCMetricsMetadata) -> bool:
        # Get file obj
        file_obj = self.igvf_client_api.get_by_id(file_id).actual_instance
        # Check if the file has QC metrics
        file_qc_metrics = file_obj.quality_metrics
        if not file_qc_metrics:
            raise ValueError(
                f"{file_metadata.content_type} file does not have QC metrics."
            )
        # Get all active QC metrics objects
        all_file_qc_objects = [self.igvf_client_api.get_by_id(
            metric_id).actual_instance for metric_id in file_qc_metrics]
        active_qc_metrics = [
            qc_obj for qc_obj in all_file_qc_objects
            if qc_obj.status not in DEPRECATED_STATUSES
        ]
        if len(active_qc_metrics) != 1:
            raise ValueError(
                f"{file_metadata.content_type} has an incorrect number of pipeline parameters."
            )
        # Get the actual QC metrics object
        curr_qc_metrics_obj = active_qc_metrics[0]
        if curr_qc_metrics_obj.type != qc_metadata.qc_metrics_type:
            raise ValueError(
                f"{file_metadata.content_type} file has incorrect QC metrics type."
            )
        if curr_qc_metrics_obj.description != qc_metadata.qc_metrics_desc:
            raise ValueError(
                f"{file_metadata.content_type} has incorrect QC metrics description."
            )
        # Check if the QC metrics object has the required attachments
        missing_attachments = []
        for attachment in qc_metadata.qc_metrics_attachments:
            if getattr(curr_qc_metrics_obj, attachment) is None:
                missing_attachments.append(attachment)
                continue
            if not curr_qc_metrics_obj.attachment:
                missing_attachments.append(attachment)
        if missing_attachments:
            raise ValueError(
                f"{file_metadata.content_type} file does not have attachment {','.join(attachment)}."
            )
        return True

    def _has_correct_document_info(self) -> bool:
        if not self.analysis_set_metadata.documents:
            raise ValueError("Analysis set does not have any documents.")

        pipeline_documents = []
        for document_id in self.analysis_set_metadata.documents:
            document_obj = self.igvf_client_api.get_by_id(
                document_id).actual_instance
            if document_obj.status not in DEPRECATED_STATUSES:
                pipeline_documents.append(
                    PipelineDocument(
                        description=document_obj.description,
                        document_type=document_obj.document_type
                    )
                )

        if not pipeline_documents:
            raise ValueError(
                "Analysis set does not have any active pipeline documents.")

        # Check for correct document count
        matching_docs = [
            doc for doc in pipeline_documents if doc == PROPER_DOCUMENT_INFO]
        if len(matching_docs) == 0:
            raise ValueError(
                "Analysis set does not have the correct pipeline parameters document.")
        elif len(matching_docs) > 1:
            raise ValueError(
                "Analysis set has multiple active pipeline parameters documents.")

        return True

    def _check_single_rna_main_file(self, file_id: str, file_metadata: FileMetadata) -> PerObjectQAResults:
        # Check file metadata
        try:
            self._has_correct_file_metadata(
                file_metadata, PROPER_RNA_MATRIX_FILES)
            metadata_result = True
        except ValueError as e:
            metadata_result = str(e)

        # Check submitted file name
        try:
            self._has_correct_submitted_file_name(file_metadata.file_type)
            file_name_result = True
        except ValueError as e:
            file_name_result = str(e)

        # Check QC metrics
        try:
            self._has_correct_qc_metrics(file_id=file_id,
                                         file_metadata=file_metadata,
                                         qc_metadata=PROPER_QC_METRICS['gene_counts'])
            qc_result = True
        except ValueError as e:
            qc_result = str(e)

        return PerObjectQAResults(
            file_metadata_check=metadata_result,
            qc_metrics_check=qc_result,
            submitted_file_name_check=file_name_result,
            document_info_check=None
        )

    def _check_single_atac_main_file(self, file_metadata: FileMetadata) -> PerObjectQAResults:
        """Check a single main file (bam, h5ad, or tar) for correctness."""
        # Check file metadata
        try:
            self._has_correct_file_metadata(
                file_metadata, PROPER_ATAC_OUTPUT_FILES[file_metadata.file_type])
            metadata_result = True
        except ValueError as e:
            metadata_result = str(e)

        # Check submitted file name
        try:
            self._has_correct_submitted_file_name(file_metadata.file_type)
            file_name_result = True
        except ValueError as e:
            file_name_result = str(e)

        # Check QC metrics
        try:
            self._has_correct_qc_metrics(
                file_metadata, PROPER_QC_METRICS[file_metadata.file_type])
            qc_result = True
        except ValueError as e:
            qc_result = str(e)

        return PerObjectQAResults(
            file_metadata_check=metadata_result,
            qc_metrics_check=qc_result,
            submitted_file_name_check=file_name_result,
            document_info_check=None
        )

    def check_one_analysis_set(self):
        if self.analysis_set_metadata.has_rnaseq:
            for file_format in PROPER_RNA_MATRIX_FILES:
                curr_rna_file_metadata = self.analysis_set_metadata.files.get(
                    file_format)
                if not curr_rna_file_metadata:
                    raise ValueError(
                        f"Analysis set does not have the required RNAseq file: {file_format}.")
                if len(curr_rna_file_metadata) != 1:
                    raise ValueError(
                        f"Analysis set has multiple {file_format} files, expected only one.")
