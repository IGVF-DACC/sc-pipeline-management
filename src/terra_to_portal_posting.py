import igvf_utils as iu
from igvf_utils.connection import Connection
from pathlib import Path
import pandas as pd
import igvf_and_terra_api_tools as api_tools
import re


# TODO: 
# 1) build a post message parsing to get alignment file to pass into index file and tabular files as derived_from
# 2) 


# NOTE: File aliases will just be analysis set accession + column header name

terra_output_table_column_types = {
    'alignment_file': [('atac_bam', 'bam')],    # ATAC
    'document': [('atac_bam_log', 'txt'),
                 ('atac_chromap_barcode_metadata', 'tsv'),
                 ('atac_snapatac2_barcode_metadata', 'tsv'),
                 ('csv_summary', 'csv'),
                 ('html_summary', 'html'),
                 ('joint_barcode_metadata', 'csv'),
                 ('rna_barcode_metadata', 'tsv'),
                 ('rna_log', 'txt')
                 ],
    'index_file': [('atac_filter_fragments_index', 'tbi')],     # ATAC
    'tabular_file': [('atac_filter_fragments', 'tsv')],     # ATAC
    'matrix_file': [('rna_aggregated_counts_h5ad', 'h5ad'),     # All RNA
                    ('rna_mtx_tar', 'tar'),
                    ('rna_mtxs_h5ad', 'h5ad'),
                    ('rna_kb_output', 'tar')
                    ]
}

# terra_data_filename_parsers = {
#     'rna_mtxs_h5ad': re.compile(r'\.rna\.align\.kb\.(.*)\.count\_matrix\.h5ad'),
#     'rna_kb_output': re.compile(r'\.rna\.align\.kb\.(.*)\.tar\.gz'),
#     'atac_bam': re.compile(r'\.atac\.align\.k4\.(genome_name)\.bam')
# }

accession_headers_by_assay_types = {'atac': ['atac_read1_accessions', 'atac_read2_accessions'],
                                    'rna': ['rna_read1_accessions', 'rna_read2_accessions']
                                    }

genome_assembly_info = {'hg38': ('GRCh38', 'IGVFFI0653VCGH'), 'mm39': ('GRCm39', 'IGVFFI9282QLXO')}


def get_lab_and_award(terra_data_record: pd.Series, igvf_api) -> tuple:
    analysis_set_objs = igvf_api.analysis_sets(accession=[terra_data_record['analysis_set_acc']]).graph[0]
    return (analysis_set_objs.award, analysis_set_objs.lab)


def get_seqfile_accs_from_table(terra_data_record: pd.Series, seqfile_acc_cols: list):
    seqfile_accessions = []
    for seqfile_col in seqfile_acc_cols:
        seqfile_accessions.extend(eval(terra_data_record[seqfile_col]))
    return seqfile_accessions


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
        seqfile_obj = igvf_api.sequence_files(accession=[seq_file_acc]).graph[0]
        access_status += seqfile_obj.controlled_access
    if access_status > 0:
        return True
    else:
        return False


def get_seqfile_access_lvl_and_access(terra_data_record: pd.Series, assay_type: str, igvf_api):
    seqfile_accessions = get_seqfile_accs_from_table(terra_data_record=terra_data_record, seqfile_acc_cols=accession_headers_by_assay_types[assay_type])
    output_data_controlled_access = get_access_status(sequence_file_accessions=seqfile_accessions, igvf_api=igvf_api)
    return (output_data_controlled_access, seqfile_accessions)


def get_gspath_and_alias(terra_data_record: pd.Series, col_name: str, lab: str) -> tuple:
    gs_cloud_link = terra_data_record[col_name]
    file_alias = f'{lab.split("/")[-2]}:{terra_data_record["analysis_set_acc"]}_{col_name}_uniform-pipeline'
    return (gs_cloud_link, file_alias)


def alignment_file_payload_template(igvf_award: str,
                                    igvf_lab: str,
                                    igvf_aliases: list,
                                    assembly: str,
                                    controlled_access: str,
                                    derived_from: list,
                                    file_set: str,
                                    md5sum: str,
                                    submitted_file_name: str,
                                    ref_files: list
                                    ):
    alignment_file_payload = {'award': igvf_award,
                              'lab': igvf_lab,
                              'aliases': igvf_aliases,
                              'assembly': assembly,
                              'content_type': 'aligments',
                              'controlled_access': controlled_access,
                              'derived_from': derived_from,     # NOTE: This is from the input sequence files
                              'file_format': 'bam',
                              'file_set': file_set,
                              'md5sum': md5sum,
                              'submitted_file_name': submitted_file_name,
                              'reference_files': ref_files
                              }
    return alignment_file_payload


def document_payload_template(igvf_award: str,
                              igvf_aliases: list,
                              submitted_file_name: str
                              ):
    document_payload = {'award': igvf_award,
                        'aliases': igvf_aliases,
                        'document_type': 'quality control report',
                        'attachment': submitted_file_name
                        }
    return document_payload


def index_file_payload_template(igvf_award: str,
                                igvf_aliases: list,
                                analysis_step_version: str,
                                content_type: str,
                                file_format: str,
                                file_set: str,
                                igvf_lab: str,
                                md5sum: str,
                                derived_from: list,
                                submitted_file_name: str,
                                ref_files: list
                                ):
    index_file_payload = {'award': igvf_award,
                          'aliases': igvf_aliases,
                          'analysis_step_version': analysis_step_version,
                          'content_type': content_type,    # NOTE: sparse gene count matrix
                          'derived_from': derived_from,    # TODO:this needs to be collected after alignment files are uploaded
                          'file_format': file_format,
                          'file_set': file_set,
                          'lab': igvf_lab,
                          'md5sum': md5sum,
                          'submitted_file_name': submitted_file_name,
                          'reference_files': ref_files
                          }
    return index_file_payload


def matrix_file_payload_template(igvf_award: str,
                                 igvf_aliases: list,
                                #  analysis_step_version: str,
                                 content_type: str,
                                 file_format: str,
                                 file_set: str,
                                 igvf_lab: str,
                                 md5sum: str,
                                 derived_from: list,
                                 submitted_file_name: str,
                                 principle_dim: str,
                                 secondary_dims: list
                                 ):
    matrix_file_payload = {'award': igvf_award,
                           'aliases': igvf_aliases,
                        #    'analysis_step_version': analysis_step_version,
                           'content_type': content_type,
                           'derived_from': derived_from,
                           'file_format': file_format,
                           'file_set': file_set,
                           'lab': igvf_lab,
                           'md5sum': md5sum,
                           'submitted_file_name': submitted_file_name,
                           'principle_dimension': principle_dim,   # TODO: Dependent on the output
                           'secondary_dimensions': secondary_dims,    # TODO: Dependent on the output
                           }
    return matrix_file_payload


def tabular_file_payload_template(igvf_award: str,
                                  igvf_aliases: list,
                                  content_type: str,
                                  controlled_access: bool,
                                  file_format: str,
                                  file_set: str,
                                  igvf_lab: str,
                                  md5sum: str,
                                  derived_from: list,
                                  submitted_file_name: str,
                                  transcriptome_annotation: str
                                  ):
    tab_file_payload = {'award': igvf_award,
                        'aliases': igvf_aliases,
                        'content_type': content_type,
                        'controlled_access': controlled_access,
                        'derived_from': derived_from,   # TODO:this needs to be collected after alignment files are uploaded
                        'file_format': file_format,
                        'file_set': file_set,
                        'lab': igvf_lab,
                        'md5sum': md5sum,
                        'submitted_file_name': submitted_file_name
                        }
    return tab_file_payload


def single_post_to_portal(file_type: str, igvf_data_payload: dict, igvf_utils_api, upload_file: bool = False):
    igvf_data_payload[Connection.PROFILE_KEY] = file_type
    _schema_property = igvf_utils_api.get_profile_from_payload(igvf_data_payload).properties
    stdout = igvf_utils_api.post(igvf_data_payload, upload_file=upload_file, return_original_status_code=True)
    return (stdout[0]['accession'], stdout[1])


def post_all_documents_to_portal(terra_data_record: pd.Series,
                                 lab: str,
                                 award: str,
                                 igvf_utils_api,
                                 upload_file: bool = False
                                 ):
    curr_post_summary = {}
    for (col_header, _curr_file_format) in terra_output_table_column_types['document']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        # curr_doc_payload = document_payload_template(igvf_award=award, igvf_aliases=curr_file_alias, submitted_file_name=curr_gs_cloud_link)
        curr_doc_payload = dict(igvf_award=award, igvf_aliases=curr_file_alias, submitted_file_name=curr_gs_cloud_link)
        curr_new_acc, curr_post_status = single_post_to_portal(file_type='document', igvf_data_payload=curr_doc_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]
    return curr_post_summary


def post_all_rna_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, seqfile_accs: list, igvf_utils_api, upload_file: bool = False):
    curr_post_summary = {}
    for (col_header, curr_file_format) in terra_output_table_column_types['matrix_file']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_md5sum = api_tools.calculate_gsutil_hash(file_path=curr_gs_cloud_link)
        curr_content_type = 'sparse gene count matrix'
        curr_seqfile_accs = get_seqfile_accs_from_table(terra_data_record=terra_data_record, seqfile_acc_cols=accession_headers_by_assay_types['rna'])
        # curr_mtx_file_payload = matrix_file_payload_template(igvf_award=award,
        #                                                      igvf_lab=lab,
        #                                                      igvf_aliases=curr_file_alias,
        #                                                      content_type=curr_content_type,
        #                                                      md5sum=curr_md5sum,
        #                                                      file_format=curr_file_format,
        #                                                      derived_from=curr_seqfile_accs,
        #                                                      submitted_file_name=curr_gs_cloud_link,
        #                                                      file_set=terra_data_record['analysis_set_acc'],
        #                                                      principle_dims='cell',
        #                                                      secondary_dims=['gene']
        #                                                      )
        curr_mtx_file_payload = dict(igvf_award=award,
                                     igvf_lab=lab,
                                     igvf_aliases=curr_file_alias,
                                     content_type=curr_content_type,
                                     md5sum=curr_md5sum,
                                     file_format=curr_file_format,
                                     derived_from=curr_seqfile_accs,
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     principle_dims='cell',
                                     secondary_dims=['gene'],
                                     
                                     )
        curr_new_acc, curr_post_status = single_post_to_portal(file_type='matrix_file', igvf_data_payload=curr_mtx_file_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]
    return curr_post_summary


def post_all_atac_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, seqfile_accs: list, igvf_api, igvf_utils_api, upload_file: bool = False):
    curr_post_summary = {}
    # Need to post alignment files first
    curr_alignment_file_accs = []
    for (col_header, curr_file_format) in terra_output_table_column_types['alignment_file']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_md5sum = api_tools.calculate_gsutil_hash(file_path=curr_gs_cloud_link)
        curr_alignfile_ctrl_access, curr_seqfile_accs = get_seqfile_access_lvl_and_access(terra_data_record=terra_data_record, assay_type='atac', igvf_api=igvf_api)
        # curr_alignment_payload = alignment_file_payload_template(igvf_award=lab,
        #                                                          igvf_lab=lab,
        #                                                          igvf_aliases=curr_file_alias,
        #                                                          assembly=genome_assembly_info[terra_data_record['Genome']][0],
        #                                                          controlled_access=curr_alignfile_ctrl_access,
        #                                                          derived_from=curr_seqfile_accs,
        #                                                          file_set=terra_data_record['analysis_set_acc'],
        #                                                          md5sum=curr_md5sum,
        #                                                          submitted_file_name=curr_gs_cloud_link,
        #                                                          ref_files=genome_assembly_info[terra_data_record['Genome']][1],
        #                                                          )
        curr_alignment_payload = dict(igvf_award=lab,
                                      igvf_lab=lab,
                                      igvf_aliases=curr_file_alias,
                                      assembly=genome_assembly_info[terra_data_record['Genome']][0],
                                      controlled_access=curr_alignfile_ctrl_access,
                                      derived_from=curr_seqfile_accs,
                                      file_set=terra_data_record['analysis_set_acc'],
                                      md5sum=curr_md5sum,
                                      submitted_file_name=curr_gs_cloud_link,
                                      ref_files=genome_assembly_info[terra_data_record['Genome']][1]
                                      )
        curr_new_acc, curr_post_status = single_post_to_portal(file_type='alignment_file', igvf_data_payload=curr_alignment_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_alignment_file_accs.append(curr_new_acc)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]

    # Post Tabular Files
    for (col_header, curr_file_format) in terra_output_table_column_types['tabular_file']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_md5sum = api_tools.calculate_gsutil_hash(file_path=curr_gs_cloud_link)
        # curr_tab_file_payload = tabular_file_payload_template(igvf_award=lab,
        #                                                       igvf_lab=lab,
        #                                                       igvf_aliases=curr_file_alias,
        #                                                       content_type='fragments',
        #                                                       file_format=curr_file_format,
        #                                                       md5sum=curr_md5sum,
        #                                                       derived_from=list(set(curr_alignment_file_accs)),
        #                                                       submitted_file_name=curr_gs_cloud_link,
        #                                                       transcriptome_annotation=
        #                                                       )
        curr_tab_file_payload = dict(igvf_award=lab,
                                     igvf_lab=lab,
                                     igvf_aliases=curr_file_alias,
                                     content_type='fragments',
                                     file_format=curr_file_format,
                                     md5sum=curr_md5sum,
                                     derived_from=list(set(curr_alignment_file_accs)),
                                     submitted_file_name=curr_gs_cloud_link
                                     )
        curr_new_acc, curr_post_status = single_post_to_portal(file_type='tabular_file', igvf_data_payload=curr_tab_file_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]
        return curr_post_summary
