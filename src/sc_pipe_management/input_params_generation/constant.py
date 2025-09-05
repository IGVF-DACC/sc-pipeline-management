import dataclasses
import re


# Convert portal read_names to pipeline input param names
READ_NAME_TO_READ_TYPE_MAP = {'Read 1': 'read1',
                              'Read 2': 'read2',
                              'Barcode index': 'barcode'
                              }

# Matching assay term name to assay type
# NOTE: This list may need to be periodically updated.
ASSAY_NAMES_CONVERSION_REF = {
    '/assay-terms/OBI_0002764/': 'atac',    # single-cell ATAC-seq
    '/assay-terms/OBI_0002762/': 'atac',     # single-nucleus ATAC-seq
    '/assay-terms/OBI_0003109/': 'rna',     # single-nucleus RNA sequencing assay
    '/assay-terms/OBI_0002631/': 'rna',  # single-cell RNA sequencing assay
}

# Preferred assay titles that will require ATAC-RNA onlist mapping
ASSAYS_NEED_ONLIST_MAPPING = ['10x multiome', '10x multiome with MULTI-seq']

# For URL request
BASE_IGVF_PORTAL_URL = 'https://api.data.igvf.org'

# Read ID Regex
IGVF_ACCESSION_REGEX = re.compile(r'(IGVFFI[A-Z0-9]+|IGVFDS[A-Z0-9]+)')

# seqspec index tool
ASSAY_TYPE_TO_TOOL_FORMAT = {'rna': 'kb-single', 'atac': 'chromap'}

# File deprecated status
FILE_DEPRECATED_STATUSES = ['revoked', 'deleted', 'replaced', 'archived']


@dataclasses.dataclass(frozen=True)
class RunReferenceFiles:
    """Dataclass to hold pipeline run reference files."""
    chromap_index: str
    genome_fasta: str
    kb_index: str
    genome_ref: str


# Run reference files for different taxa
TAXA_TO_GENOME_REF_FILES = {
    'Homo sapiens': RunReferenceFiles(
        chromap_index='https://api.data.igvf.org/reference-files/IGVFFI7969JLFC/@@download/IGVFFI7969JLFC.tar.gz',
        genome_fasta='https://api.data.igvf.org/reference-files/IGVFFI0653VCGH/@@download/IGVFFI0653VCGH.fasta.gz',
        kb_index='https://api.data.igvf.org/reference-files/IGVFFI9561BASO/@@download/IGVFFI9561BASO.tar.gz',
        genome_ref='gs://broad-buenrostro-pipeline-genome-annotations/IGVF_human_v43/IGVF_human_v43_Homo_sapiens_genome_files_hg38_v43.tsv'
    ),
    'Mus musculus': RunReferenceFiles(
        chromap_index='https://api.data.igvf.org/reference-files/IGVFFI5593VLWB/@@download/IGVFFI5593VLWB.tar.gz',
        genome_fasta='https://api.data.igvf.org/reference-files/IGVFFI9282QLXO/@@download/IGVFFI9282QLXO.fasta.gz',
        kb_index='https://api.data.igvf.org/reference-files/IGVFFI5078MNED/@@download/IGVFFI5078MNED.tar.gz',
        genome_ref='gs://broad-buenrostro-pipeline-genome-annotations/IGVF_mouse_v32/IGVF_mouse_v32_Mus_musculus_genome_files_mm39_v32.tsv'
    )
}


class BadDataException(Exception):
    """Custom exception for bad data."""
    pass
