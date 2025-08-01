"""
Constants for Terra to Portal posting operations.
"""

import re

# Hard code lab and award for sandbox and production
POSTING_LAB = '/labs/igvf-dacc-processing-pipeline/'
POSTING_AWARD = '/awards/HG012012/'

# Hard code for Terra output table column types
TERRA_OUTPUT_TABLE_COLUMN_TYPES = {
    # ATACseq data
    'alignment_file': {
        'atac_bam': {
            'file_format': 'bam',
            'description': 'Raw Aligned bam file from Chromap',
            'content_type': 'alignments'
        }
    },
    'index_file': {
        'atac_fragments_index': {
            'file_format': 'tbi',
            'description': 'Raw Fragment file index from Chromap',
            'content_type': 'index'
        },
        'atac_bam_index': {
            'file_format': 'bai',
            'description': 'Raw Aligned bam file index from Chromap',
            'content_type': 'index'
        },
    },
    'tabular_file': {
        'atac_fragments': {
            'file_format': 'bed',
            'description': 'Raw Fragment file from Chromap',
            'content_type': 'fragments',
            'file_format_specifications': ['buenrostro-bernstein:igvf-single-cell-pipeline-fragment-file-specification']
        }
    },
    # RNAseq data
    'matrix_file': {
        'rna_kb_output_folder_tar_gz': {
            'file_format': 'tar',
            'description': 'Raw Tarball containing all the matrices, logs, and bus files generated from kb',
            'content_type': 'kallisto single cell RNAseq output',
            'file_format_specifications': ['buenrostro-bernstein:igvf-sc-pipeline-matrix-tar-specification',
                                           'igvf:igvf-sc-pipeline-rna-tar-mtx-per-file-specification']
        },
        'rna_kb_h5ad': {
            'file_format': 'h5ad',
            'description': 'Raw h5ad containing four separated count matrices: Spliced, Unspliced, Ambiguous, and Total',
            'content_type': 'sparse gene count matrix',
            'file_format_specifications': ['buenrostro-bernstein:igvf-sc-pipeline-matrix-h5-specification']
        },
    },
}

# Hard code for columns with IGVF accessions
ACCESSION_HEADERS_BY_ASSAY_TYPES = {'atac': ['atac_read1_accessions', 'atac_read2_accessions', 'atac_barcode_accessions'],
                                    'rna': ['rna_read1_accessions', 'rna_read2_accessions', 'rna_barcode_accessions']
                                    }

# NOTE: These may end up not being needed if genome tsv reference and assembly column is available in input
# Convert Terra table's info to portal enum
GENOME_ASSEMBLY_INFO = {'Homo sapiens': 'GRCh38', 'Mus musculus': 'GRCm39'}

# Analysis step versions
ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES = {'atac': 'igvf:single-cell-uniform-pipeline-chromap-atacseq-step-v1.1.0',
                                         'rna': 'igvf:single-cell-uniform-pipeline-kalliso-bustools-rnaseq-step-v1.1.0'
                                         }

# QC objects
# Metadata: Table column names
# Attachment: keys are portal properties and values are column IDs
# Metadata: keys are JSON keys, values are portal properties
TERRA_QC_OUTPUTS = {'rnaseq': {'gene_count': {'metadata': ['rna_kb_library_qc_metrics_json', 'rna_kb_run_info_json'],
                                              'description': 'RNAseq Kallisto Bustools QC metric',
                                              'attachment': {'rnaseq_kb_info': 'rna_kb_parameters_json'},
                                              'object_type': 'single_cell_rna_seq_quality_metric',
                                              'metadata_map': {'numRecords': 'n_records',   # starting here is inspect.json
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
                                                               }
                                              }
                               },
                    'atacseq': {'alignment': {'metadata': None,
                                              'description': 'ATACseq chromap alignment QC metric',
                                              'attachment': {'atac_bam_summary_stats': 'atac_bam_summary_stats'},
                                              'object_type': 'single_cell_atac_seq_quality_metric'
                                              },
                                'fragment': {'metadata': ['atac_fragments_metrics'],
                                             'description': 'ATACseq chromap fragments QC metric',
                                             'attachment': {'atac_fragment_summary_stats': 'atac_fragments_barcode_summary',
                                                            'atac_fragments_alignment_stats': 'atac_fragments_alignment_stats',
                                                            },
                                             'object_type': 'single_cell_atac_seq_quality_metric',
                                             'metadata_map': {'Number_of_reads': 'n_reads',
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
                                                              }
                                             }
                                }
                    }

# IGVF file url parsing regex for accession
IGVF_URL_PATH_REGEX = re.compile(
    r'\/.*-files\/(IGVF[A-Z0-9]+|TST[A-Z0-9]+)\/@@download')

# GS file path regex
GS_FILE_PATH_REGEX = re.compile(
    r'gs://([a-z0-9\-]+)/submissions/final-outputs/([a-z0-9\-]+)/single_cell_pipeline/([a-z0-9\-]+)/[a-z\-]+/[a-z\_]+/([a-z0-9\-]+)/.*')
