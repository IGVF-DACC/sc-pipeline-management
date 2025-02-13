import igvf_utils as iu
from igvf_utils.connection import Connection
from pathlib import Path
import pandas as pd
import igvf_and_terra_api_tools as api_tools
import logging
import multiprocessing
from functools import partial
import datetime
import igvf_utils.exceptions as iu_exceptions
import os
import dataclasses


# TODO:
# 1) Figure out how to post documents from non-local options (likely need to work from a Google container)
# 2) Build in some filtering or quality check if a pipeline run fails, so there is no file to post.


# NOTE: File aliases will just be analysis set accession + column header name

# A possbily hacky dict to match output file columns with data types and file format
TERRA_OUTPUT_TABLE_COLUMN_TYPES = {
    # ATAC
    'alignment_file': [('atac_bam', 'bam', 'Raw Aligned bam file from Chromap')],
    'document': [('atac_bam_log', 'txt', 'Raw Log file from aligner'),
                 ('atac_chromap_barcode_metadata', 'tsv',
                  'Raw Per barcode alignment statistics file from Chromap'),
                 ('atac_snapatac2_barcode_metadata', 'tsv',
                  'Filtered Per barcode statistics file from SnapATAC2'),
                 ('csv_summary', 'csv', 'Filtered CSV summary report'),
                 ('html_summary', 'html', 'Filtered HTML summary report'),
                 ('joint_barcode_metadata', 'csv',
                  'Filtered Joint per barcode statistics file'),
                 ('rna_barcode_metadata', 'tsv',
                  'Raw Per barcode alignment statistics file'),
                 ('rna_log', 'txt', 'Raw Log file from kb')
                 ],
    # ATAC
    'index_file': [('atac_filter_fragments_index', 'tbi', 'Raw Fragment file index from Chromap')],
    # ATAC
    'tabular_file': [('atac_filter_fragments', 'tsv', 'Raw Fragment file from Chromap')],
    'matrix_file': [('rna_aggregated_counts_h5ad', 'h5ad', 'Raw Aggregated(Ambiguous+Spliced+Unspliced) count matrix in h5ad format'),
                    ('rna_mtx_tar', 'tar',
                     'Raw Tarball containing four separated count matrices in mtx format: Spliced, Unspliced, Ambiguous, and Total'),
                    ('rna_mtxs_h5ad', 'h5ad',
                     'Raw h5ad containing four separated count matrices: Spliced, Unspliced, Ambiguous, and Total'),
                    ('rna_kb_output', 'tar',
                     'Raw Tarball containing all the logs and bus files generated from kb')
                    ]
}

# Hard code for columns with IGVF accessions
ACCESSION_HEADERS_BY_ASSAY_TYPES = {'atac': ['atac_read1_accessions', 'atac_read2_accessions', 'atac_barcode_accessions'],
                                    'rna': ['rna_read1_accessions', 'rna_read2_accessions']
                                    }

# NOTE: These may end up not being needed if genome tsv reference and assembly column is available in input
# Convert Terra table's info to portal enum
GENOME_ASSEMBLY_INFO = {'hg38': 'GRCh38', 'mm39': 'GRCm39'}

# Genome reference files map
GENOME_REFERENCE_FILES = {'GRCh38': {'rna': ['/reference-files/TSTFI89593580/'],
                                     'atac': ['/reference-files/TSTFI36924773/']
                                     },
                          'GRCm39': {'rna': ['/reference-files/TSTFI17059281/'],
                                     'atac': ['/reference-files/TSTFI77805766/']
                                     }
                          }


# Define how to return errors
@dataclasses.dataclass(frozen=True)
class PostResult:
    col_header: str
    accession: str | None
    error: Exception | None = None

    @classmethod
    def Failure(cls, col_header, error):
        return cls(col_header, None, error)

    @classmethod
    def Success(cls, col_header, accession):
        return cls(col_header, accession, None)

    def Description(self) -> str:
        if self.error:
            return 'POST Failed: ' + str(self.error)
        return self.accession


@dataclasses.dataclass
class RunResult:
    analysis_set_acc: str
    post_results: list[PostResult]
    finish_time: str


# Helper util functions
def check_single_or_multiome(terra_data_record: pd.Series) -> list:
    """Check based on input measurement set IDs to see if the pipeline is RNAseq-only, ATACseq-only, or multiome-seq

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)

    Returns:
        list: content of the pipeline run
    """
    pipeline_output_content = []
    if len(eval(terra_data_record['atac_MeaSetIDs'])) > 0:
        pipeline_output_content.append('ATACseq')
    if len(eval(terra_data_record['rna_MeaSetIDs'])) > 0:
        pipeline_output_content.append('RNAseq')
    return pipeline_output_content


def get_award_and_lab(terra_data_record: pd.Series, igvf_api) -> tuple:
    """Get the Lab and Award info from analysis accession IDs.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        igvf_api (_type_): _description_

    Returns:
        tuple: (award, lab)
    """
    analysis_set_objs = igvf_api.analysis_sets(
        accession=[terra_data_record['analysis_set_acc']]).graph[0]
    return (analysis_set_objs.award, analysis_set_objs.lab)


def get_seqfile_accs_from_table(terra_data_record: pd.Series, seqfile_acc_cols: list) -> list:
    """Get a list of sequence file accessions (for derived_from) from Terra table.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        seqfile_acc_cols (list): columns with sequence file accessions of interest

    Returns:
        list: A list of sequence file accessions
    """
    seqfile_accessions = []
    for seqfile_col in seqfile_acc_cols:
        seqfile_accessions.extend(eval(terra_data_record[seqfile_col]))
    # Only to get the IGVF seqfile accessions from prod or sandbox
    return [seqfile_accession for seqfile_accession in seqfile_accessions if seqfile_accession.startswith(('IGVF', 'TSTF'))]


def get_access_status(sequence_file_accessions: list, igvf_api) -> bool:
    """Get file controlled_access status. If any seq data is controlled access, the output data inherit that.

    Args:
        sequence_file_accessions (list): Input sequence file accessions
        igvf_api (_type_): _description_

    Returns:
        bool: True or False
    """
    access_status = 0
    for seq_file_acc in sequence_file_accessions:
        seqfile_obj = igvf_api.sequence_files(
            accession=[seq_file_acc]).graph[0]
        access_status += seqfile_obj.controlled_access
    if access_status > 0:
        return True
    else:
        return False


def get_seqfile_access_lvl_and_accessions(terra_data_record: pd.Series, assay_type: str, igvf_api) -> tuple:
    """Get sequence files' controlled access info and accession IDs

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        assay_type (str): rna or atac
        igvf_api (_type_): _description_

    Returns:
        tuple: (controlled_access_boolean, sequence file accessions)
    """
    seqfile_accessions = get_seqfile_accs_from_table(
        terra_data_record=terra_data_record, seqfile_acc_cols=ACCESSION_HEADERS_BY_ASSAY_TYPES[assay_type])
    output_data_controlled_access = get_access_status(
        sequence_file_accessions=seqfile_accessions, igvf_api=igvf_api)
    return (output_data_controlled_access, seqfile_accessions)


def get_gspath_and_alias(terra_data_record: pd.Series, col_name: str, lab: str) -> tuple:
    """Get an output file's Google cloud address and generate a portal alias

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_name (str): Output file column header
        lab (str): data generator/submitter lab

    Returns:
        tuple: (gs cloud link, portal alias)
    """
    gs_cloud_link = terra_data_record[col_name]
    file_alias = f'{lab.split("/")[-2]}:{terra_data_record["analysis_set_acc"]}_{col_name}_uniform-pipeline'
    return (gs_cloud_link, file_alias)


# Helper function to post data object via IGVF utils
def single_post_to_portal(igvf_data_payload: dict, igvf_utils_api, upload_file: bool = False) -> str:
    """POST a single data object to the portal

    Args:
        igvf_data_payload (dict): submission POST payload
        igvf_utils_api (_type_): IU connection api
        upload_file (bool, optional): _description_. Defaults to False.

    Returns:
        str: accession generated by POST
    """
    _schema_property = igvf_utils_api.get_profile_from_payload(
        igvf_data_payload).properties
    stdout = igvf_utils_api.post(
        igvf_data_payload, upload_file=upload_file, return_original_status_code=True)
    return stdout[0]['accession']


# # Post documents to the portal
# # TODO: needs to figure out how to post documents as it is only for local files
# def post_single_document(terra_data_record: pd.Series, col_header: str, curr_description: str, lab: str, award: str, igvf_utils_api, upload_file: bool) -> tuple:
#     """_summary_

#     Args:
#         terra_data_record (pd.Series): _description_
#         col_header (str): _description_
#         curr_description (str): _description_
#         lab (str): _description_
#         award (str): _description_
#         igvf_utils_api (_type_): _description_
#         upload_file (bool): _description_

#     Returns:
#         tuple: _description_
#     """
#     curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(
#         terra_data_record=terra_data_record, col_name=col_header, lab=lab)
#     if str(curr_gs_cloud_link).startswith('gs://'):
#         curr_doc_payload = dict(lab=lab, award=award,
#                                 aliases=[curr_file_alias],
#                                 attachment={'path': curr_gs_cloud_link},
#                                 document_type='quality control report',
#                                 description=curr_description
#                                 )
#         logging.info(f'Posting {curr_file_alias} to the portal.')
#         curr_post_res, curr_post_status = single_post_to_portal(
#             igvf_data_payload=curr_doc_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
#         return col_header, curr_post_res, curr_post_status
#     else:
#         return col_header, 'POST failed', 'File path is not a gs cloud link.'


# def post_all_documents_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_utils_api, upload_file: bool) -> dict:
#     """_summary_

#     Args:
#         terra_data_record (pd.Series): _description_
#         lab (str): _description_
#         award (str): _description_
#         igvf_utils_api (_type_): _description_
#         upload_file (bool): _description_

#     Returns:
#         dict: _description_
#     """
#     curr_post_summary = {}
#     for (col_header, _curr_file_format, curr_description) in TERRA_OUTPUT_TABLE_COLUMN_TYPES['document']:
#         if terra_data_record[col_header]:
#             col_header, curr_post_res, curr_post_status = post_single_document(
#                 terra_data_record, col_header, curr_description, lab, award, igvf_utils_api, upload_file)
#             curr_post_summary[col_header] = [curr_post_res, curr_post_status]
#     return curr_post_summary

# Posting RNA data to the portal


def single_tabular_file_post(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_assembly: str, igvf_utils_api, upload_file: bool) -> tuple:
    """Post a single tabular file (for RNA output) to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_assembly (str): The genome assembly
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        tuple: file name, posted output on the portal, POST status
    """
    try:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(
            terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        curr_seqfile_accs = get_seqfile_accs_from_table(
            terra_data_record=terra_data_record, seqfile_acc_cols=ACCESSION_HEADERS_BY_ASSAY_TYPES['rna'])
        # TODO: Missing analysis_step_version
        curr_mtx_file_payload = dict(award=award,
                                     lab=lab,
                                     aliases=[curr_file_alias],
                                     content_type='sparse gene count matrix',
                                     md5sum=curr_md5sum,
                                     file_format=curr_file_format,
                                     derived_from=curr_seqfile_accs,
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     principal_dimension='cell',
                                     secondary_dimensions=['gene'],
                                     # TODO: need a proper way to look this up
                                     reference_files=GENOME_REFERENCE_FILES[curr_assembly]['rna'],
                                     description=curr_description,
                                     _profile='matrix_file')
        # Post the metadata
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_mtx_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        # Upload files
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except Exception as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_all_rna_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_utils_api, upload_file: bool) -> dict:
    """POST all RNA data to the portal (all TabFiles)

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        igvf_utils_api (_type_): _description_
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        list[PostResult]
    """
    curr_assembly = GENOME_ASSEMBLY_INFO[terra_data_record['Genome']]
    with multiprocessing.Pool() as pool:
        results = pool.starmap(partial(single_tabular_file_post, terra_data_record, lab=lab, award=award, curr_assembly=curr_assembly, igvf_utils_api=igvf_utils_api, upload_file=upload_file),
                               [(col_header, curr_file_format, curr_description) for col_header, curr_file_format, curr_description in TERRA_OUTPUT_TABLE_COLUMN_TYPES['matrix_file']])
    return results


# Post ATAC data to the portal
def post_single_alignment_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_ctrl_access: bool, curr_seqfile_accs: list, igvf_utils_api, upload_file: bool) -> tuple:
    """Post one alignment file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_seqfile_accs (list): The sequence file accessions (for derived_from)
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        tuple: file name, posted output on the portal, POST status
    """
    try:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(
            terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        curr_assembly = GENOME_ASSEMBLY_INFO[terra_data_record['Genome']]
        # TODO: Missing analysis_step_version
        curr_alignment_payload = dict(award=award,
                                      lab=lab,
                                      aliases=[curr_file_alias],
                                      assembly=curr_assembly,
                                      controlled_access=curr_ctrl_access,
                                      file_format=curr_file_format,
                                      content_type='alignments',
                                      filtered=False,
                                      derived_from=curr_seqfile_accs,
                                      file_set=terra_data_record['analysis_set_acc'],
                                      md5sum=curr_md5sum,
                                      submitted_file_name=curr_gs_cloud_link,
                                      # TODO: need a proper way to look this up
                                      reference_files=GENOME_REFERENCE_FILES[curr_assembly]['atac'],
                                      redacted=False,
                                      description=curr_description,
                                      _profile='alignment_file'
                                      )
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_alignment_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except Exception as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_single_tabular_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_ctrl_access: bool, curr_alignment_file_accs: list, igvf_utils_api, upload_file: bool) -> tuple:
    """Post one tabular file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_alignment_file_accs (list): The alignment file accessions (for derived_from)
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        tuple: file name, posted output on the portal, POST status
    """
    try:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(
            terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        # TODO: Missing analysis_step_version
        curr_tab_file_payload = dict(award=award,
                                     lab=lab,
                                     aliases=[curr_file_alias],
                                     content_type='fragments',
                                     controlled_access=curr_ctrl_access,
                                     file_format=curr_file_format,
                                     md5sum=curr_md5sum,
                                     derived_from=list(
                                         set(curr_alignment_file_accs)),
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     description=curr_description,
                                     _profile='tabular_file'
                                     )
        logging.info(f'Posting {curr_file_alias} to the portal.')
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_tab_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except Exception as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_single_index_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_ctrl_access: bool, curr_fragfile_accs: list, igvf_utils_api, upload_file: bool) -> tuple:
    """Post one index file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_fragfile_accs (list): The fragment file accessions (for derived_from)
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        tuple: file name, posted output on the portal, POST status
    """
    try:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(
            terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        curr_md5sum = api_tools.calculate_gsutil_hash(
            file_path=curr_gs_cloud_link)
        # TODO: Missing analysis_step_version
        curr_index_file_payload = dict(award=award,
                                       lab=lab,
                                       aliases=[curr_file_alias],
                                       content_type='index',
                                       controlled_access=curr_ctrl_access,
                                       file_format=curr_file_format,
                                       md5sum=curr_md5sum,
                                       derived_from=list(
                                           set(curr_fragfile_accs)),
                                       submitted_file_name=curr_gs_cloud_link,
                                       file_set=terra_data_record['analysis_set_acc'],
                                       description=curr_description,
                                       _profile='index_file'
                                       )
        logging.info(f'Posting {curr_file_alias} to the portal.')
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_index_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, accession=curr_post_res)
    except Exception as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_all_atac_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_api, igvf_utils_api, upload_file: bool) -> dict:
    """Post all ATAC data to the portal (alignment, fragment, and index files).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        dict: {column_header: [posted output on the portal, POST status]}
    """
    curr_ctrl_access, curr_seqfile_accs = get_seqfile_access_lvl_and_accessions(
        terra_data_record=terra_data_record, assay_type='atac', igvf_api=igvf_api)
    curr_post_summary = []

    with multiprocessing.Pool() as pool:
        # Post alignment files
        alignment_results = pool.starmap(partial(post_single_alignment_file, terra_data_record, lab=lab, award=award, curr_ctrl_access=curr_ctrl_access, curr_seqfile_accs=curr_seqfile_accs, igvf_utils_api=igvf_utils_api, upload_file=upload_file), [
                                         (col_header, curr_file_format, curr_description) for col_header, curr_file_format, curr_description in TERRA_OUTPUT_TABLE_COLUMN_TYPES['alignment_file']])
        # Check for alignment post files
        curr_alignment_file_accs = [res.Description(
        ) for res in alignment_results if not res.Description().startswith('POST fail')]
        curr_post_summary += alignment_results

        # Post tabular files
        tabular_results = pool.starmap(partial(post_single_tabular_file, terra_data_record, lab=lab, award=award, curr_ctrl_access=curr_ctrl_access, curr_alignment_file_accs=curr_alignment_file_accs,
                                       igvf_utils_api=igvf_utils_api, upload_file=upload_file), [(col_header, curr_file_format, curr_description) for col_header, curr_file_format, curr_description in TERRA_OUTPUT_TABLE_COLUMN_TYPES['tabular_file']])
        curr_fragfile_accs = [res.Description(
        ) for res in tabular_results if not res.Description().startswith('POST fail')]
        curr_post_summary += tabular_results

        # Post index files
        index_results = pool.starmap(partial(post_single_index_file, terra_data_record, lab=lab, award=award, curr_ctrl_access=curr_ctrl_access, curr_fragfile_accs=curr_fragfile_accs, igvf_utils_api=igvf_utils_api, upload_file=upload_file), [
                                     (col_header, curr_file_format, curr_description) for col_header, curr_file_format, curr_description in TERRA_OUTPUT_TABLE_COLUMN_TYPES['index_file']])
        curr_post_summary += index_results
    return curr_post_summary


# Full data posting from a single pipeline run
def post_all_data_from_one_run(terra_data_record: pd.Series, igvf_api, igvf_utils_api, upload_file: bool) -> dict:
    """Post all single cell uniform pipeline output to the portal from Terra.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        dict: {column_header: [posted output on the portal, POST status]}
    """
    logging.info(
        f'Posting data from {terra_data_record["analysis_set_acc"]} to the portal.')
    # Check if the pipeline is RNAseq-only, ATACseq-only, or multiome-seq
    pipeline_run_output = check_single_or_multiome(
        terra_data_record=terra_data_record)
    # Posting
    post_results = []
    # Get lab and award
    pipeline_data_award, pipeline_data_lab = get_award_and_lab(
        terra_data_record=terra_data_record, igvf_api=igvf_api)
    # If RNAseq results are present, post them
    if 'RNAseq' in pipeline_run_output:
        post_results += post_all_rna_data_to_portal(terra_data_record=terra_data_record,
                                                    lab=pipeline_data_lab,
                                                    award=pipeline_data_award,
                                                    igvf_utils_api=igvf_utils_api,
                                                    upload_file=upload_file
                                                    )
    # If ATACseq results are present, post them
    if 'ATACseq' in pipeline_run_output:
        post_results += post_all_atac_data_to_portal(terra_data_record=terra_data_record,
                                                     lab=pipeline_data_lab,
                                                     award=pipeline_data_award,
                                                     igvf_api=igvf_api,
                                                     igvf_utils_api=igvf_utils_api,
                                                     upload_file=upload_file
                                                     )
    # TODO: Post all relevaent documents
    # document_output_post_log = post_all_documents_to_portal(terra_data_record=terra_data_record,
    #                                                         lab=pipeline_data_lab,
    #                                                         award=pipeline_data_award,
    #                                                         igvf_utils_api=igvf_utils_api,
    #                                                         upload_file=upload_file
    #                                                         )
    # output_post_log.update(document_output_post_log)
    return RunResult(
        analysis_set_acc=terra_data_record['analysis_set_acc'],
        post_results=post_results,
        finish_time=str(datetime.datetime.now()))


# Post the entire batch of pipeline outputs to the portal
def post_all_successful_runs(full_terra_data_table: pd.DataFrame, igvf_api, igvf_utils_api, upload_file: bool) -> list:
    """Post all successful runs from a Terra data table to the portal.

    Args:
        full_terra_data_table (pd.DataFrame): The Complete Terra data table
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal

    Returns:
        list: PostResult objects
    """
    results = []
    for _, curr_pipeline in full_terra_data_table.iterrows():
        results.append(post_all_data_from_one_run(
            terra_data_record=curr_pipeline, igvf_api=igvf_api, igvf_utils_api=igvf_utils_api, upload_file=upload_file))
    return results


def summarize_post_status(post_results: list) -> pd.DataFrame:
    """Summarize the POST status of all data from a single pipeline run.

    Args:
        post_status_df (pd.DataFrame): The POST status DataFrame

    Returns:
        pd.DataFrame: The POST status summary DataFrame
    """
    post_status_summary = {}
    for result in post_results:
        accession_results = post_status_summary.setdefault(
            result.analysis_set_acc, {'time_stamp': result.finish_time})
        for post_res in result.post_results:
            accession_results[post_res.col_header] = post_res.Description()
    post_status_summary_table = pd.DataFrame(
        post_status_summary).transpose().fillna('')
    post_status_summary_table['post_status_fail'] = post_status_summary_table.apply(
        lambda x: x.str.startswith(('POST Failed')), axis=1).sum(axis=1)
    return post_status_summary_table


def add_post_status_summary_to_output_data_table(full_terra_data_table: pd.DataFrame, post_status_df: pd.DataFrame) -> pd.DataFrame:
    """Add the POST status summary to the output data table.

    Args:
        full_terra_data_table (pd.DataFrame): The Complete Terra data table
        post_status_df (pd.DataFrame): The POST status DataFrame

    Returns:
        pd.DataFrame: The output data table with POST status summary
    """
    return pd.merge(left=full_terra_data_table, right=post_status_df[['time_stamp', 'post_status_fail']], left_on='analysis_set_acc', right_index=True)


def save_pipeline_postres_tables(pipeline_postres_table: pd.DataFrame, updated_full_data_table: pd.DataFrame, output_dir: str):
    """Save the pipeline input table to a TSV file.

    Args:
        pipeline_input_table (pd.DataFrame): The pipeline input table
        output_dir (str): The output directory
    """
    curr_datetime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    pipeline_postres_table.to_csv(os.path.join(
        output_dir, f'single-cell_uniform_pipeline_postres_detail_{curr_datetime}.tsv'), sep='\t')

    updated_full_data_table.to_csv(os.path.join(
        output_dir, f'full_terra_data_table_with_post_status_{curr_datetime}.tsv'), sep='\t')
