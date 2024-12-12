import igvf_utils as iu
from igvf_utils.connection import Connection
from pathlib import Path
import pandas as pd
import igvf_and_terra_api_tools as api_tools


# TODO: 
# 1) Figure out how to post documents from non-local options (likely need to work from a Google container)
# 2) 


# NOTE: File aliases will just be analysis set accession + column header name

# A possbily hacky dict to match output file columns with data types and file format
terra_output_table_column_types = {
    'alignment_file': [('atac_bam', 'bam', '[Raw] Aligned bam file from Chromap')],    # ATAC
    'document': [('atac_bam_log', 'txt', '[Raw] Log file from aligner'),
                 ('atac_chromap_barcode_metadata', 'tsv', '[Raw] Per barcode alignment statistics file from Chromap'),
                 ('atac_snapatac2_barcode_metadata', 'tsv', '[Filtered] Per barcode statistics file from SnapATAC2'),
                 ('csv_summary', 'csv', '[Filtered] CSV summary report'),
                 ('html_summary', 'html', '[Filtered] HTML summary report'),
                 ('joint_barcode_metadata', 'csv', '[Filtered] Joint per barcode statistics file'),
                 ('rna_barcode_metadata', 'tsv', '[Raw] Per barcode alignment statistics file'),
                 ('rna_log', 'txt', '[Raw] Log file from kb')
                 ],
    'index_file': [('atac_filter_fragments_index', 'tbi', '	[Raw] Fragment file index from Chromap')],     # ATAC
    'tabular_file': [('atac_filter_fragments', 'tsv', '[Raw] Fragment file from Chromap')],     # ATAC
    'matrix_file': [('rna_aggregated_counts_h5ad', 'h5ad', '[Raw] Aggregated(Ambiguous+Spliced+Unspliced) count matrix in h5ad format'),
                    ('rna_mtx_tar', 'tar', '[Raw] Tarball containing four separated count matrices in mtx format: Spliced, Unspliced, Ambiguous, and Total'),
                    ('rna_mtxs_h5ad', 'h5ad', '[Raw] h5ad containing four separated count matrices: Spliced, Unspliced, Ambiguous, and Total'),
                    ('rna_kb_output', 'tar', '[Raw] Tarball containing all the logs and bus files generated from kb')
                    ]
}

# Hard code for columns with IGVF accessions
accession_headers_by_assay_types = {'atac': ['atac_read1_accessions', 'atac_read2_accessions'],
                                    'rna': ['rna_read1_accessions', 'rna_read2_accessions']
                                    }

# Convert Terra table's info to portal enum
genome_assembly_info = {'hg38': 'GRCh38', 'mm39': 'GRCm39'}


def get_award_and_lab(terra_data_record: pd.Series, igvf_api) -> tuple:
    """Get the Lab and Award info from analysis accession IDs.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        igvf_api (_type_): _description_

    Returns:
        tuple: (award, lab)
    """
    analysis_set_objs = igvf_api.analysis_sets(accession=[terra_data_record['analysis_set_acc']]).graph[0]
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


def get_seqfile_access_lvl_and_access(terra_data_record: pd.Series, assay_type: str, igvf_api) -> tuple:
    """Get sequence files' controlled access info and accession IDs

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        assay_type (str): rna or atac
        igvf_api (_type_): _description_

    Returns:
        tuple: (controlled_access_boolean, sequence file accessions)
    """
    seqfile_accessions = get_seqfile_accs_from_table(terra_data_record=terra_data_record, seqfile_acc_cols=accession_headers_by_assay_types[assay_type])
    output_data_controlled_access = get_access_status(sequence_file_accessions=seqfile_accessions, igvf_api=igvf_api)
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


def single_post_to_portal(igvf_data_payload: dict, igvf_utils_api, upload_file: bool = False) -> tuple:
    """POST a single data object to the portal

    Args:
        igvf_data_payload (dict): submission POST payload
        igvf_utils_api (_type_): IU connection api
        upload_file (bool, optional): _description_. Defaults to False.

    Returns:
        tuple: (accession generated by POST, POST status)
    """
    _schema_property = igvf_utils_api.get_profile_from_payload(igvf_data_payload).properties
    stdout = igvf_utils_api.post(igvf_data_payload, upload_file=upload_file, return_original_status_code=True)
    return (stdout[0]['accession'], stdout[1])


# TODO: Cannot post anything that is not local
def post_all_documents_to_portal(terra_data_record: pd.Series,
                                 lab: str,
                                 award: str,
                                 igvf_utils_api,
                                 upload_file: bool = False
                                 ) -> dict:
    """Post all documents generated to the portal

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): data generator/submitter lab
        award (str): IGVF award ID
        igvf_utils_api (_type_): IU connection api
        upload_file (bool, optional): _description_. Defaults to False.

    Returns:
        dict: {column_header: [new accession, post status]}
    """
    curr_post_summary = {}
    for (col_header, _curr_file_format, curr_description) in terra_output_table_column_types['document']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_doc_payload = dict(lab=lab, award=award,
                                aliases=[curr_file_alias],
                                attachment={'path': curr_gs_cloud_link},    # TODO: This will fail with a non-local file path
                                document_type='quality control report',
                                description=curr_description
                                )
        curr_new_acc, curr_post_status = single_post_to_portal(file_type='document', igvf_data_payload=curr_doc_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]
    return curr_post_summary


def post_all_rna_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_utils_api, upload_file: bool = False) -> dict:
    """Post all RNA output data to IGVF portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): data generator/submitter
        award (str): IGVF award ID
        igvf_utils_api (_type_): IU connect api
        upload_file (bool, optional): _description_. Defaults to False.

    Returns:
        dict: {column_header: [new accession, post status]}
    """
    curr_post_summary = {}
    for (col_header, curr_file_format, curr_description) in terra_output_table_column_types['matrix_file']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_md5sum = api_tools.calculate_gsutil_hash(file_path=curr_gs_cloud_link)
        curr_seqfile_accs = get_seqfile_accs_from_table(terra_data_record=terra_data_record, seqfile_acc_cols=accession_headers_by_assay_types['rna'])
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
                                     reference_files=['TSTFI36924773'],     # NOTE: This probably needs somewhere to either specify on Terra or by hard coder
                                     description=curr_description,
                                     _profile='matrix_file'
                                     )
        curr_new_acc, curr_post_status = single_post_to_portal(igvf_data_payload=curr_mtx_file_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]
    return curr_post_summary


def post_all_atac_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_api, igvf_utils_api, upload_file: bool = False) -> dict:
    """POST all ATAC data to IGVF portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): data generator/submitter
        award (str): IGVF award ID
        igvf_api (_type_): IGVF python client api
        igvf_utils_api (_type_): IU connect api
        upload_file (bool, optional): _description_. Defaults to False.

    Returns:
        dict: {column_header: [new accession, post status]}
    """
    # Get file access and sequence files
    curr_ctrl_access, curr_seqfile_accs = get_seqfile_access_lvl_and_access(terra_data_record=terra_data_record, assay_type='atac', igvf_api=igvf_api)
    curr_post_summary = {}
    # Need to post alignment files first
    curr_alignment_file_accs = []
    for (col_header, curr_file_format, curr_description) in terra_output_table_column_types['alignment_file']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_md5sum = api_tools.calculate_gsutil_hash(file_path=curr_gs_cloud_link)
        # curr_alignfile_ctrl_access, curr_seqfile_accs = get_seqfile_access_lvl_and_access(terra_data_record=terra_data_record, assay_type='atac', igvf_api=igvf_api)
        curr_alignment_payload = dict(award=award,
                                      lab=lab,
                                      aliases=[curr_file_alias],
                                      assembly=genome_assembly_info[terra_data_record['Genome']],
                                      controlled_access=curr_ctrl_access,
                                      file_format=curr_file_format,
                                      content_type='alignments',
                                      filtered=False,
                                      derived_from=curr_seqfile_accs,
                                      file_set=terra_data_record['analysis_set_acc'],
                                      md5sum=curr_md5sum,
                                      submitted_file_name=curr_gs_cloud_link,
                                      reference_files=['TSTFI36924773'],
                                      redacted=False,   # TODO: Need to figure this one out
                                      description=curr_description,
                                      _profile='alignment_file'
                                      )
        curr_new_acc, curr_post_status = single_post_to_portal(igvf_data_payload=curr_alignment_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_alignment_file_accs.append(curr_new_acc)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]

    # Post Tabular Files
    curr_fragfile_accs = []
    for (col_header, curr_file_format, curr_description) in terra_output_table_column_types['tabular_file']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_md5sum = api_tools.calculate_gsutil_hash(file_path=curr_gs_cloud_link)
        # curr_tabfile_ctrl_access, curr_seqfile_accs = get_seqfile_access_lvl_and_access(terra_data_record=terra_data_record, assay_type='atac', igvf_api=igvf_api)
        curr_tab_file_payload = dict(award=award,
                                     lab=lab,
                                     aliases=[curr_file_alias],
                                     content_type='fragments',
                                     controlled_access=curr_ctrl_access,
                                     file_format=curr_file_format,
                                     md5sum=curr_md5sum,
                                     derived_from=list(set(curr_alignment_file_accs)),
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     description=curr_description,
                                     _profile='tabular_file'
                                     )
        curr_new_acc, curr_post_status = single_post_to_portal(igvf_data_payload=curr_tab_file_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_fragfile_accs.append(curr_new_acc)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]

    # Post fragment index files
    for (col_header, curr_file_format, curr_description) in terra_output_table_column_types['index_file']:
        curr_gs_cloud_link, curr_file_alias = get_gspath_and_alias(terra_data_record=terra_data_record, col_name=col_header, lab=lab)
        curr_md5sum = api_tools.calculate_gsutil_hash(file_path=curr_gs_cloud_link)
        curr_index_file_payload = dict(award=award,
                                       lab=lab,
                                       aliases=[curr_file_alias],
                                       content_type='index',
                                       controlled_access=curr_ctrl_access,
                                       file_format=curr_file_format,
                                       md5sum=curr_md5sum,
                                       derived_from=list(set(curr_fragfile_accs)),
                                       submitted_file_name=curr_gs_cloud_link,
                                       file_set=terra_data_record['analysis_set_acc'],
                                       description=curr_description,
                                       _profile='index_file'
                                       )
        curr_new_acc, curr_post_status = single_post_to_portal(igvf_data_payload=curr_index_file_payload, igvf_utils_api=igvf_utils_api, upload_file=upload_file)
        curr_post_summary[col_header] = [curr_new_acc, curr_post_status]
    return curr_post_summary


def post_all_data_from_one_run(terra_data_table: pd.Series, igvf_api, igvf_utils_api, upload_file):
    output_post_log = {}
    pipeline_data_award, pipeline_data_lab = get_award_and_lab(terra_data_record=terra_data_table, igvf_api=igvf_api)
    rna_output_post_log = post_all_atac_data_to_portal(terra_data_record=terra_data_table,
                                                       lab=pipeline_data_lab,
                                                       award=pipeline_data_award,
                                                       igvf_api=igvf_api,
                                                       igvf_utils_api=igvf_utils_api,
                                                       upload_file=upload_file
                                                       )
    output_post_log.update(rna_output_post_log)
    atac_output_post_log = post_all_atac_data_to_portal(terra_data_record=terra_data_table,
                                                        lab=pipeline_data_lab,
                                                        award=pipeline_data_award,
                                                        igvf_api=igvf_api,
                                                        igvf_utils_api=igvf_utils_api,
                                                        upload_file=upload_file
                                                        )
    output_post_log.update(atac_output_post_log)
    return output_post_log
