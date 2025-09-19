from typing import TypeAlias
import igvf_client

# Possible terms unsed in dataset descriptions and aliases
POSSIBLE_SINGLE_CELL_PIPELINE_PHRASES = ['scpipe',
                                         'single cell',
                                         'single cell uniform pipeline',
                                         'uniform pipeline',
                                         'sc-pipe',
                                         'scpipeline',
                                         'sc-uniform-pipeline',
                                         'sc-pipeline',
                                         'uniform-pipeline']

# Measurement sets will be excluded if they have any of these audits
EXCLUDED_NONCOMP_AUDITS = ['missing sequence specification',
                           'missing barcode replacement file',
                           'missing read names',
                           'missing barcode onlist']

EXCLUDED_ERROR_AUDITS = ['upload status not validated',
                         'unexpected barcode onlist',
                         'inconsistent sequence specifications',
                         'inconsistent preferred assay titles']


# Easier type annotation
FileObjTypes: TypeAlias = igvf_client.models.tabular_file.TabularFile \
    | igvf_client.models.matrix_file.MatrixFile \
    | igvf_client.models.alignment_file.AlignmentFile

QcObjTypes: TypeAlias = igvf_client.models.single_cell_rna_seq_quality_metric.SingleCellRnaSeqQualityMetric \
    | igvf_client.models.single_cell_atac_seq_quality_metric.SingleCellAtacSeqQualityMetric

IgvfApiType: TypeAlias = igvf_client.api.igvf_api.IgvfApi


# Some constances
RNASEQ_REF_FILES = ['/reference-files/IGVFFI5078MNED/',
                    '/reference-files/IGVFFI9561BASO/']

ATAC_REF_FILES = [
    sorted(['/reference-files/IGVFFI0653VCGH/',
           '/reference-files/IGVFFI7969JLFC/']),
    sorted(['/reference-files/IGVFFI9282QLXO/',
           '/reference-files/IGVFFI5593VLWB/'])
]

ATAC_ANALYSIS_STEP_VERSION = '/analysis-step-versions/894a5695-77b0-4770-b296-194ef8ab8115/'

RNA_ANALYSIS_STEP_VERSION = '/analysis-step-versions/5c1aba4e-2fbf-47a3-9c3c-592fb67495e3/'
