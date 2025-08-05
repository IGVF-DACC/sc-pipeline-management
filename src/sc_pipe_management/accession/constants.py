"""
Constants for Terra to Portal posting operations.
"""

import dataclasses

# The lab and award to post data under
OUTPUT_SUBMITTER_INFO = {
    'lab': '/labs/igvf-dacc-processing-pipeline/',
    'award': '/awards/HG012012/'
}


# Convert Terra table's info to portal enum
GENOME_ASSEMBLY_INFO = {'Homo sapiens': 'GRCh38',
                        'Mus musculus': 'GRCm39'
                        }


# AnalysisStepVersion
ANALYSIS_STEP_VERSION = {
    'atacseq': 'igvf:single-cell-uniform-pipeline-chromap-atacseq-step-v1.1.0',
    'rnaseq': 'igvf:single-cell-uniform-pipeline-kallisto-bustools-rnaseq-step-v1.1.0'
}


@dataclasses.dataclass(frozen=True)
class FileObjMetadata:
    """Payload metadata for file objects that will be accessioned to the IGVF portal."""
    file_format: str    # file format enum used for the output file
    description: str    # description for the output file
    content_type: str   # content type enum used for the output file
    # file format specification documents' aliases for the output file
    file_specifications: list
    analysis_step_version: str  # ASV used for for the output file
    assay_type: str

# The following is a mapping of file types to their metadata.
# The dict key is the name of the Terra data table column that hosts the GCP path to the file.


# Metadata for alignment files, following the IGVF portal schema.
ALIGNMENT_FILETYPES = {
    'atac_bam': FileObjMetadata(file_format='bam',
                                description='Raw Aligned bam file from Chromap',
                                content_type='alignments',
                                file_specifications=[
                                    'buenrostro-bernstein:igvf-sc-pipeline-alignment-bam-specification'],
                                analysis_step_version=ANALYSIS_STEP_VERSION['atacseq'],
                                assay_type='atac'
                                )
}

# Metadata for index files, following the IGVF portal schema
INDEX_FILETYPES = {
    'atac_fragments_index': FileObjMetadata(file_format='tbi',
                                            description='Raw Fragment file index from Chromap',
                                            content_type='index',
                                            file_specifications=[
                                                'buenrostro-bernstein:igvf-sc-pipeline-fragment-file-specification'],
                                            analysis_step_version=ANALYSIS_STEP_VERSION['atacseq'],
                                            assay_type='atac'
                                            ),
    'atac_bam_index': FileObjMetadata(file_format='bai',
                                      description='Raw Aligned bam file index from Chromap',
                                      content_type='index',
                                      file_specifications=[
                                          'buenrostro-bernstein:igvf-sc-pipeline-alignment-bam-index-specification'],
                                      analysis_step_version=ANALYSIS_STEP_VERSION['atacseq'],
                                      assay_type='atac'
                                      )
}

# Metadata for tabular files, following the IGVF portal schema
TABULAR_FILETYPES = {
    'atac_fragments': FileObjMetadata(file_format='bed',
                                      description='Raw Fragment file from Chromap',
                                      content_type='fragments',
                                      file_format_specifications=[
                                          'buenrostro-bernstein:igvf-single-cell-pipeline-fragment-file-specification'],
                                      analysis_step_version=ANALYSIS_STEP_VERSION['atacseq'],
                                      assay_type='atac'
                                      )
}

# Metadata for matrix files, following the IGVF portal schema
MATRIX_FILETYPES = {
    'rna_kb_output_folder_tar_gz': FileObjMetadata(
        file_format='tar',
        description='Raw Tarball containing all the matrices, logs, and bus files generated from kb',
        content_type='kallisto single cell RNAseq output',
        file_format_specifications=['buenrostro-bernstein:igvf-sc-pipeline-matrix-tar-specification',
                                    'igvf:igvf-sc-pipeline-rna-tar-mtx-per-file-specification'],
        analysis_step_version=ANALYSIS_STEP_VERSION['rnaseq'],
        assay_type='rna'
    ),
    'rna_kb_h5ad': FileObjMetadata(
        file_format='h5ad',
        description='Raw h5ad containing four separated count matrices: Spliced, Unspliced, Ambiguous, and Total',
        content_type='sparse gene count matrix',
        file_format_specifications=[
            'buenrostro-bernstein:igvf-sc-pipeline-matrix-h5-specification'],
        analysis_step_version=ANALYSIS_STEP_VERSION['rnaseq'],
        assay_type='rna'
    )
}


# Dataclass for QC info map that will be used to generate QC metrics payloads
@dataclasses.dataclass(frozen=True)
class QCInfoMap:
    """Dataclass to hold QC information for posting to the IGVF portal."""
    metadata: list  # The output files will be parsed into a JSON schema payload
    object_type: str    # The object type for the QC metrics
    attachment: dict    # The output files that will be uploaded as attachments
    metadata_map: dict  # The metadata map to convert the QC JSON files to the portal schema
    description: str    # Description of the QC metrics
    analysis_step_version: str  # The analysis step version for the QC metrics


# QC objects using QCInfoMap dataclass
TERRA_QC_OUTPUTS = {
    'rnaseq': {
        'gene_count': QCInfoMap(
            metadata=['rna_kb_library_qc_metrics_json',
                      'rna_kb_run_info_json'],
            description='RNAseq Kallisto Bustools QC metric',
            attachment={'rnaseq_kb_info': 'rna_kb_parameters_json'},
            object_type='single_cell_rna_seq_quality_metric',
            metadata_map={
                'numRecords': 'n_records',   # starting here is inspect.json
                'numReads': 'n_reads',
                'numBarcodes': 'n_barcodes',
                'medianReadsPerBarcode': 'median_reads_per_barcode',
                'meanReadsPerBarcode': 'mean_reads_per_barcode',
                'numUMIs': 'total_umis',
                'numBarcodeUMIs': 'n_barcode_umis',
                'medianUMIsPerBarcode': 'median_umis_per_barcode',
                'meanUMIsPerBarcode': 'mean_umis_per_barcode',
                'gtRecords': 'gt_records',
                'numBarcodesOnOnlist': 'num_barcodes_on_onlist',
                'percentageBarcodesOnOnlist': 'percentage_barcodes_on_onlist',
                'numReadsOnOnlist': 'num_reads_on_onlist',
                'percentageReadsOnOnlist': 'percentage_reads_on_onlist',
                'n_targets': 'n_targets',    # starting here is run_info.json
                'n_bootstraps': 'n_bootstraps',
                'n_processed': 'n_processed',
                'n_pseudoaligned': 'n_pseudoaligned',
                'n_unique': 'n_unique',
                'p_pseudoaligned': 'p_pseudoaligned',
                'p_unique': 'p_unique',
                'index_version': 'index_version',
                'k-mer length': 'kmer_length'
            },
            analysis_step_version=ANALYSIS_STEP_VERSION['rnaseq']
        )
    },
    'atacseq': {
        'alignment': QCInfoMap(
            metadata=None,
            description='ATACseq chromap alignment QC metric',
            attachment={'atac_bam_summary_stats': 'atac_bam_summary_stats'},
            object_type='single_cell_atac_seq_quality_metric',
            metadata_map={},  # No metadata mapping for alignment QC
            analysis_step_version=ANALYSIS_STEP_VERSION['atacseq']
        ),
        'fragment': QCInfoMap(
            metadata=['atac_fragments_metrics'],
            description='ATACseq chromap fragments QC metric',
            attachment={
                'atac_fragment_summary_stats': 'atac_fragments_barcode_summary',
                'atac_fragments_alignment_stats': 'atac_fragments_alignment_stats',
            },
            object_type='single_cell_atac_seq_quality_metric',
            metadata_map={
                'Number_of_reads': 'n_reads',
                'Number_of_mapped_reads': 'n_mapped_reads',
                'Number_of_uniquely_mapped_reads': 'n_uniquely_mapped_reads',
                'Number_of_reads_have_multi-mappings': 'n_reads_with_multi_mappings',
                'Number_of_candidates': 'n_candidates',
                'Number_of_mappings': 'n_mappings',
                'Number_of_uni-mappings': 'n_uni_mappings',
                'Number_of_multi-mappings': 'n_multi_mappings',
                'Number_of_barcodes_in_whitelist': 'n_barcodes_on_onlist',
                'Number_of_corrected_barcodes': 'n_corrected_barcodes',
                'Number_of_output_mappings_(passed_filters)': 'n_output_mappings',
                'uni-mappings': 'uni_mappings',
                'multi-mappings': 'multi_mappings',
                'total': 'total',
                'percentage_duplicates': 'pct_duplicates',
            },
            analysis_step_version=ANALYSIS_STEP_VERSION['atacseq']
        )
    }
}
