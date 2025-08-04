import pandas as pd
import igvf_and_terra_api_tools as api_tools
import logging
import multiprocessing
from functools import partial
import datetime
import os
import dataclasses
import requests
import time
from itertools import chain

from constants import (
    TERRA_OUTPUT_TABLE_COLUMN_TYPES,
    ACCESSION_HEADERS_BY_ASSAY_TYPES,
    ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES,
    TERRA_QC_OUTPUTS,
    POSTING_LAB,
    POSTING_AWARD
)

from sc_pipe_management.accession.parse_terra_metadata import (
    parse_workflow_uuids_from_gs_path,
    parse_igvf_accessions_from_urls,
    get_genome_assembly,
    check_single_or_multiome,
    get_seqfile_accs_from_table,
    get_seqfile_access_lvl_and_accessions,
    get_file_aliases,
    download_qc_file_from_gcp,
    read_json_file,
    mk_qc_obj_aliases,
    get_existing_analysis_set_docs
)


# TODO:\
# fapi._session() can time out, how to handle that?
# Figure out how to resume/fix bad posts


# NOTE:
# 1) File aliases are lab:gcp-bucket_submissionID_workflowID_subworkflowID_OutputFileColumn
# 2) RNA doesn't have an assembly (matrix file)

# Define how to return errors
@dataclasses.dataclass(frozen=True)
class PostResult:
    col_header: str
    uuid: str | None
    error: Exception | None = None

    @classmethod
    def Failure(cls, col_header, error):
        return cls(col_header, None, error)

    @classmethod
    def Success(cls, col_header, uuid):
        return cls(col_header, uuid, None)

    def Description(self) -> str:
        if self.error:
            return 'POST Failed: ' + str(self.error)
        return self.uuid


@dataclasses.dataclass
class RunResult:
    analysis_set_acc: str
    post_results: list[PostResult]
    finish_time: str


class PortalConflictError(Exception):
    """Custom exception for when a new data object conflicts with an existing one in the portal, by md5 or aliases."""
    pass


def mk_doc_payload(lab: str, award: str, doc_aliases: list, local_file_path: str) -> dict:
    """Create a document payload for posting to the portal.

    Args:
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        doc_aliases (list): Document aliases
        local_file_path (str): Path to the document attachment

    Returns:
        dict: Document payload ready for posting
    """
    return dict(
        aliases=doc_aliases,
        lab=lab,
        award=award,
        document_type='pipeline parameters',
        description='Terra workflow configuration for the single-cell pipeline run',
        attachment={'path': local_file_path},
        _profile='document'
    )


# Functions to create payloads
def mk_anaset_docs_patching_payload(doc_aliases: list, analysis_set_acc: str, igvf_utils_api) -> dict:
    """Create a payload for patching the analysis set with the document accession.

    Args:
        doc_aliases (str): Document aliases
        analysis_set_acc (str): Analysis set accession
        igvf_utils_api (_type_): igvf utils API client

    Returns:
        dict: Payload for patching the analysis set
    """
    existing_docs = get_existing_analysis_set_docs(analysis_set_acc=analysis_set_acc,
                                                   igvf_utils_api=igvf_utils_api)
    if set(doc_aliases).issubset(set(existing_docs)):
        return {}
    new_docs_aliases = list(set(doc_aliases).union(set(existing_docs)))
    return {
        'documents': new_docs_aliases,
        igvf_utils_api.IGVFID_KEY: f"/analysis-sets/{analysis_set_acc}/",
        'uniform_pipeline_status': 'completed',
        '_profile': 'analysis_set'
    }


# Functions to do IGVF portal posting and error handling
def get_conflict_object_id(igvf_utils_api, igvf_data_payload: dict) -> str | None:
    """Get the ID of a conflicting object, if it exists by aliases or MD5."""
    if igvf_data_payload.get('aliases'):
        existing_file_by_alias = igvf_utils_api.get(
            f"aliases:{igvf_data_payload['aliases'][0]}")
        if existing_file_by_alias:
            return existing_file_by_alias['@id']
    if igvf_data_payload.get('md5sum'):
        existing_file_by_md5 = igvf_utils_api.get(
            f"md5:{igvf_data_payload['md5sum']}")
        if existing_file_by_md5:
            return existing_file_by_md5['@id']
    return None


def single_post_to_portal(igvf_data_payload: dict, igvf_utils_api, upload_file: bool = False, resumed_posting: bool = False) -> str:
    """POST a single data object to the portal, after checking for conflicts by MD5 or aliases.

    Args:
        igvf_data_payload (dict): submission POST payload
        igvf_utils_api (_type_): IGVF utils API client
        upload_file (bool, optional): whether to upload data file. Defaults to False.
        resumed_posting (bool, optional): whether to patch an existing post. Defaults to False. When False, if any new object submitted clashes with an existing object, the returned existing object record ID will have extra content, thus, subsequent linked objects will fail. If True, the existing object record id will be returned, so subsequent objects can be posted if there is any linkage such as derived_from or quality_metric_of. The True option is for cases where some part of the data finished posting, but some other files failed due to other reasons, so the unaccessioned objects can posted again without having to re-post the entire run.

    Returns:
        str: UUID generated by POST
    """
    _schema_property = igvf_utils_api.get_profile_from_payload(
        igvf_data_payload).properties
    # First check payload MD5 and aliases to see if the object already exists
    existing_object_id = get_conflict_object_id(
        igvf_utils_api, igvf_data_payload)
    # If conflict exists and not expecting it to, raise an error
    if (existing_object_id is not None) and (not resumed_posting):
        raise PortalConflictError(
            f'{igvf_data_payload.get("_profile").capitalize().replace("_", " ")} failed to accession due to conflict with an existing object: {existing_object_id}.')
    # If conflict exists and patching is expected, return the existing object ID
    if (existing_object_id is not None) and (resumed_posting):
        return existing_object_id
    # Post new object to portal if no conflict is expected
    stdout = igvf_utils_api.post(
        igvf_data_payload, upload_file=upload_file, return_original_status_code=True, truncate_long_strings_in_payload_log=True)
    return stdout[0]['uuid']


def midway_failing_chained_posting(will_fail_col_headers: list[str]) -> PostResult:
    """Return a PostResult indicating that the posting failed midway due to some objects already posted."""
    midway_fail_results = []
    for col_header in will_fail_col_headers:
        midway_fail_results.append(PostResult.Failure(
            col_header=col_header,
            error=f'Cannot post {col_header} because its linked object failed to post.'))
    return midway_fail_results


# Post workflow configuration document to the portal
def post_single_document_to_portal(terra_data_record: pd.Series, lab: str, award: str, config_file_path: str, doc_aliases: list, igvf_utils_api, resumed_posting: bool = False) -> list[PostResult]:
    """Post workflow configuration document to the portal and link it to the analysis set.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        igvf_utils_api (_type_): igvf utils API client
        config_file_path (str): Path to the workflow configuration file
        doc_aliases (list): Document aliases
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        list[PostResult]: list(PostResult(col_header='workflow_config_document', accession=IGVF UUID, error=error if any))
    """
    workflow_config_post_results = []
    try:
        # Make a document object payload
        doc_payload = mk_doc_payload(
            lab=lab,
            award=award,
            doc_aliases=doc_aliases,
            local_file_path=config_file_path
        )
        # Post the document to the portal
        curr_post_res = single_post_to_portal(igvf_data_payload=doc_payload,
                                              igvf_utils_api=igvf_utils_api,
                                              resumed_posting=resumed_posting,
                                              upload_file=False)
        workflow_config_post_results.append(PostResult.Success(
            col_header='workflow_config_document', uuid=curr_post_res))
        # Patch the analysis set with the document accession
        try:
            analysis_set_doc_patch_payload = mk_anaset_docs_patching_payload(
                doc_aliases=doc_aliases,
                analysis_set_acc=terra_data_record['analysis_set_acc'],
                igvf_utils_api=igvf_utils_api
            )
            if analysis_set_doc_patch_payload:
                igvf_utils_api.patch(mk_anaset_docs_patching_payload(
                    doc_aliases=doc_aliases,
                    analysis_set_acc=terra_data_record['analysis_set_acc'],
                    igvf_utils_api=igvf_utils_api
                ))
                workflow_config_post_results.append(PostResult.Success(
                    col_header='analysis_set_workflow_config_patch', uuid=terra_data_record["analysis_set_acc"]))
            else:
                print(
                    '>>>> No update needed for analysis set with workflow config document.')
                workflow_config_post_results.append(PostResult.Success(
                    col_header='analysis_set_workflow_config_patch', uuid='No update needed'))
        # Log error of patching analysis set with workflow config document
        except (requests.exceptions.HTTPError, PortalConflictError) as e:
            workflow_config_post_results.append(PostResult.Failure(
                col_header='analysis_set_workflow_config_patch', error=e))
    # Log error of posting workflow config document
    except (requests.exceptions.HTTPError, PortalConflictError) as e:
        workflow_config_post_results.append(PostResult.Failure(
            col_header='workflow_config_document', error=e))
    # Return the post results
    return workflow_config_post_results


# Post QC metrics to the portal
def post_single_qc_metric_to_portal(terra_data_record: pd.Series, qc_data_info: dict, qc_of_file_accs: list, lab: str, award: str, analysis_step_version: str, qc_prefix: str, qc_obj_aliases: list, igvf_utils_api, output_root_dir: str = '/igvf/data/', resumed_posting: bool = False) -> list[PostResult]:
    """Post a single QC metric to the portal. If the QC_OF objects fail to post due to conflicts or other errors, it will return a failure result.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        qc_data_info (dict): The QC data collection
        qc_of_file_accs (list): The QC file accessions
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        analysis_step_version (str): The analysis step version
        qc_prefix (str): The prefix for the QC object type (fragment, gene count, etc.)
        qc_obj_aliases (list): The list of QC object aliases
        igvf_utils_api (_type_): IGVF python client api
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        PostResult: list(PostResult(col_header=qc_obj_type, accession=..., error=...))
    """
    try:
        # NOTE: GCP VM persistent disk is /igvf/data
        curr_qc_file_dir = os.path.join(
            output_root_dir, 'qc_metrics', terra_data_record['analysis_set_acc'])
        if not os.path.exists(curr_qc_file_dir):
            os.makedirs(curr_qc_file_dir)
        # Base QC object payload
        qc_payload = dict(lab=lab,
                          award=award,
                          aliases=qc_obj_aliases,
                          analysis_step_version=analysis_step_version,
                          description=qc_data_info['description'],
                          quality_metric_of=qc_of_file_accs,
                          _profile=qc_data_info['object_type']
                          )
        # Add attachments (key is portal property name, value is table colname)
        # Do this first because some QC obj only has attachments
        if qc_data_info['attachment'] is not None:
            for key, value in qc_data_info['attachment'].items():
                # Download the file first
                curr_attachment_file = download_qc_file_from_gcp(
                    gs_file_path=terra_data_record[value], downloaded_dir=curr_qc_file_dir)
                qc_payload.update({key: {'path': curr_attachment_file}})
        # Add on other metadata
        if qc_data_info['metadata'] is not None:
            # Downloadec JSON file paths
            curr_metadata_files = [download_qc_file_from_gcp(
                gs_file_path=terra_data_record[metadata_col], downloaded_dir=curr_qc_file_dir) for metadata_col in qc_data_info['metadata']]
            # Parse Json files
            qc_json_parsed = read_json_file(
                json_file_paths=curr_metadata_files)
            # Update base qc payload
            for key, value in qc_data_info['metadata_map'].items():
                qc_payload.update({value: qc_json_parsed[key]})
        curr_post_res = single_post_to_portal(igvf_data_payload=qc_payload,
                                              igvf_utils_api=igvf_utils_api,
                                              upload_file=False,
                                              resumed_posting=resumed_posting
                                              )
        # Post result
        return PostResult.Success(col_header=f'{qc_prefix}:{qc_data_info["object_type"]}', uuid=curr_post_res)
    except (requests.exceptions.HTTPError, PortalConflictError) as e:
        return PostResult.Failure(col_header=f'{qc_prefix}:{qc_data_info["object_type"]}', error=e)


class PostService:
    """Service class to handle posting data to the portal."""

    def __init__(self, igvf_utils_api):
        self.igvf_utils_api = igvf_utils_api

    def post(self, payload: dict, upload_file: bool = False, resumed_posting: bool = False) -> str:
        """Post data to the portal."""
        return single_post_to_portal(
            igvf_data_payload=payload,
            igvf_utils_api=self.igvf_utils_api,
            upload_file=upload_file,
            resumed_posting=resumed_posting
        )


class TerraMetadata:
    """Class to hold Terra metadata for posting to the portal."""

    def __init__(self, workflow_id: str, workflow_name: str, workflow_version: str):
        self.workflow_id = workflow_id
        self.workflow_name = workflow_name
        self.workflow_version = workflow_version

    def get_workflow_info(self) -> dict:
        return string_value_1

    def get_file_path(self, col_header: str) -> str:
        """Get the file path for a given column header."""
        return f"gs://{self.workflow_id}/{col_header}"


class MatrixFile:

    def __init__(self, file_metadata: FileMetadata,):
        self.file_metadata = file_metadata

    def payload_mtx():
        return dict(prop1=self.file_metadata.file_format, pro2=, pro3=)


# Posting RNA data to the portal
def post_single_matrix_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_derived_from: list, curr_description: str, lab: str, award: str, curr_content_type: str, curr_file_specification: list, igvf_utils_api, upload_file: bool, resumed_posting: bool = False) -> PostResult:
    """Post a single tabular file (for RNA output) to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_assembly (str): The genome assembly
        curr_content_type (str): The content type
        curr_file_specification (list): The file format specification
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        PostResult: PostResult(col_header=col_header, accession=..., error=...)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_aliases = get_file_aliases(col_header=col_header,
                                             lab=lab,
                                             terra_data_record=terra_data_record,
                                             terra_uuids=parse_workflow_uuids_from_gs_path(
                                                 gs_path=curr_gs_cloud_link).aliases()
                                             )
        # Calculate md5sum
        curr_md5sum = api_tools.calculate_gsfile_hex_hash(
            file_path=curr_gs_cloud_link)
        # Only one reference file
        reference_files = parse_igvf_accessions_from_urls(
            igvf_file_urls=[terra_data_record['kb_index']])
        curr_mtx_file_payload = dict(award=award,
                                     lab=lab,
                                     analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES['rna'],
                                     aliases=curr_file_aliases,
                                     content_type=curr_content_type,
                                     md5sum=curr_md5sum,
                                     file_format=curr_file_format,
                                     derived_from=curr_derived_from,
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     principal_dimension='cell',
                                     secondary_dimensions=['gene'],
                                     filtered=False,
                                     reference_files=reference_files,
                                     description=curr_description,
                                     file_format_specifications=curr_file_specification,
                                     _profile='matrix_file')
        # Post the metadata
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_mtx_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False, resumed_posting=resumed_posting)
        # Upload files
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, uuid=curr_post_res)
    except (requests.exceptions.HTTPError, PortalConflictError) as e:
        return PostResult.Failure(col_header=col_header, error=e)


class Submission:
    terra_data_record: TerraDataRecord
    lab: str
    award: str


class TerraDataRecord:
    data: pd.Series


def post_all_rna_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/', resumed_posting: bool = False) -> list[PostResult]:
    """POST all RNA data as matrix files to the portal (H5AD, tarball, and QC metrics).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        igvf_utils_api (_type_): _description_
        upload_file (bool): Whether to upload the file to the portal
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        list[PostResult]: list[PostResult], e.g., [PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    # Get sequence file accessions
    curr_seqfile_accs = get_seqfile_accs_from_table(
        terra_data_record=terra_data_record, seqfile_acc_cols=ACCESSION_HEADERS_BY_ASSAY_TYPES['rna'])

    # Get seqspec and barcode replacement file accessions (barcode file can be NaN on Terra)
    curr_seqspec_and_barcode_replacement = parse_igvf_accessions_from_urls(
        igvf_file_urls=[terra_data_record['rna_seqspec_urls'], terra_data_record.get('barcode_replacement_file', None)])

    # Combine seqfile and seqspec into one list as derived_from
    all_derived_from = list(
        set(chain(curr_seqfile_accs, curr_seqspec_and_barcode_replacement)))

    # Prepare arguments for multiprocessing
    tasks = [
        (
            terra_data_record,
            col_header,
            file_info,
            all_derived_from,
            file_info['description'],
            lab,
            award,
            file_info['content_type'],
            file_info['file_format_specifications'],
            igvf_utils_api,
            upload_file,
            resumed_posting
        )
        for col_header, file_info in TERRA_OUTPUT_TABLE_COLUMN_TYPES['matrix_file'].items()
    ]

    # Use multiprocessing to post RNA matrix files
    with multiprocessing.Pool() as pool:
        results = pool.starmap(post_single_matrix_file, tasks)

    # Post RNA QC metrics
    qc_prefix = 'gene_count'
    # Get already posted RNA mtx file UUIDs
    posted_rna_accessions = [
        res.uuid for res in results if res.error is None]

    # If Post failed
    if not posted_rna_accessions:
        results.extend(midway_failing_chained_posting(
            will_fail_col_headers=[f'{qc_prefix}_qc']
        ))
        return results

    curr_workflow_config = parse_workflow_uuids_from_gs_path(
        gs_path=terra_data_record['rna_kb_h5ad'])
    # Generate QC object aliases
    qc_obj_aliases = mk_qc_obj_aliases(curr_workflow_config=curr_workflow_config,
                                       analysis_set_acc=terra_data_record['analysis_set_acc'],
                                       qc_prefix=qc_prefix,
                                       lab=lab)
    # Post QC metrics to the portal
    qc_metric_post_result = post_single_qc_metric_to_portal(terra_data_record=terra_data_record,
                                                            qc_data_info=TERRA_QC_OUTPUTS['rnaseq'][qc_prefix],
                                                            qc_of_file_accs=posted_rna_accessions,
                                                            lab=lab,
                                                            award=award,
                                                            analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                                                'rna'],
                                                            qc_prefix=qc_prefix,
                                                            qc_obj_aliases=qc_obj_aliases,
                                                            igvf_utils_api=igvf_utils_api,
                                                            output_root_dir=output_root_dir,
                                                            resumed_posting=resumed_posting
                                                            )
    # Add QC post report
    results.append(qc_metric_post_result)
    return results


# Post ATAC data to the portal
def post_single_alignment_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_ctrl_access: bool, curr_derived_from: list, igvf_utils_api, upload_file: bool, resumed_posting: bool = False) -> PostResult:
    """Post one alignment file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        PostResult: e.g., PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_aliases = get_file_aliases(col_header=col_header,
                                             lab=lab,
                                             terra_data_record=terra_data_record,
                                             terra_uuids=parse_workflow_uuids_from_gs_path(
                                                 gs_path=curr_gs_cloud_link).aliases()
                                             )
        curr_md5sum = api_tools.calculate_gsfile_hex_hash(
            file_path=curr_gs_cloud_link)
        # Get reference files
        reference_files = parse_igvf_accessions_from_urls(
            igvf_file_urls=[terra_data_record['chromap_index'], terra_data_record['genome_fasta']])
        # Make payload
        curr_alignment_payload = dict(award=award,
                                      lab=lab,
                                      analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES['atac'],
                                      aliases=curr_file_aliases,
                                      controlled_access=curr_ctrl_access,
                                      file_format=curr_file_format,
                                      content_type='alignments',
                                      filtered=False,
                                      derived_from=curr_derived_from,
                                      file_set=terra_data_record['analysis_set_acc'],
                                      md5sum=curr_md5sum,
                                      submitted_file_name=curr_gs_cloud_link,
                                      reference_files=reference_files,
                                      redacted=False,
                                      description=curr_description,
                                      _profile='alignment_file'
                                      )
        # Post
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_alignment_payload, igvf_utils_api=igvf_utils_api, upload_file=False, resumed_posting=resumed_posting)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, uuid=curr_post_res)
    except (requests.exceptions.HTTPError, PortalConflictError) as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_single_tabular_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_assembly: str, curr_ctrl_access: bool, curr_derived_from: list, curr_file_specification: list, igvf_utils_api, upload_file: bool, resumed_posting: bool = False) -> PostResult:
    """Post one fragment tabular file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_assembly (str): The genome assembly
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The sequence file accessions (for derived_from)
        curr_file_specification (list): The file format specification
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        PostResult: e.g., PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_aliases = get_file_aliases(col_header=col_header,
                                             lab=lab,
                                             terra_data_record=terra_data_record,
                                             terra_uuids=parse_workflow_uuids_from_gs_path(
                                                 gs_path=curr_gs_cloud_link).aliases()
                                             )
        curr_md5sum = api_tools.calculate_gsfile_hex_hash(
            file_path=curr_gs_cloud_link)
        curr_tab_file_payload = dict(award=award,
                                     lab=lab,
                                     analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES['atac'],
                                     aliases=curr_file_aliases,
                                     content_type='fragments',
                                     controlled_access=curr_ctrl_access,
                                     file_format=curr_file_format,
                                     md5sum=curr_md5sum,
                                     derived_from=curr_derived_from,
                                     submitted_file_name=curr_gs_cloud_link,
                                     file_set=terra_data_record['analysis_set_acc'],
                                     description=curr_description,
                                     filtered=False,
                                     assembly=curr_assembly,
                                     file_format_specifications=curr_file_specification,
                                     _profile='tabular_file'
                                     )
        # Work on schema dependency
        if curr_file_format == 'bed':
            curr_tab_file_payload['file_format_type'] = 'bed3+'
        # Post
        logging.info(f'Posting {curr_file_aliases} to the portal.')
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_tab_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False, resumed_posting=resumed_posting)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, uuid=curr_post_res)
    except (requests.exceptions.HTTPError, PortalConflictError) as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_single_index_file(terra_data_record: pd.Series, col_header: str, curr_file_format: str, curr_description: str, lab: str, award: str, curr_ctrl_access: bool, curr_derived_from: list, igvf_utils_api, upload_file: bool, resumed_posting: bool = False) -> PostResult:
    """Post one index file to the portal.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        col_header (str): A data table column header, which is the output file name
        curr_file_format (str): The output file format
        curr_description (str): The output file description
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The file that this index file is index for (for derived_from)
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        PostResult: e.g., PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)
    """
    try:
        # Get the gs cloud link
        curr_gs_cloud_link = terra_data_record[col_header]
        if not str(curr_gs_cloud_link).startswith('gs://'):
            raise Exception('File path is not a gs cloud link.')
        # Parse gs cloud link to get file alias
        curr_file_aliases = get_file_aliases(col_header=col_header,
                                             lab=lab,
                                             terra_data_record=terra_data_record,
                                             terra_uuids=parse_workflow_uuids_from_gs_path(
                                                 gs_path=curr_gs_cloud_link).aliases()
                                             )
        curr_md5sum = api_tools.calculate_gsfile_hex_hash(
            file_path=curr_gs_cloud_link)
        curr_index_file_payload = dict(award=award,
                                       lab=lab,
                                       analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                           'atac'],
                                       aliases=curr_file_aliases,
                                       content_type='index',
                                       controlled_access=curr_ctrl_access,
                                       file_format=curr_file_format,
                                       md5sum=curr_md5sum,
                                       derived_from=list(
                                           set(curr_derived_from)),
                                       submitted_file_name=curr_gs_cloud_link,
                                       file_set=terra_data_record['analysis_set_acc'],
                                       description=curr_description,
                                       _profile='index_file'
                                       )
        logging.info(f'Posting {curr_file_aliases} to the portal.')
        curr_post_res = single_post_to_portal(
            igvf_data_payload=curr_index_file_payload, igvf_utils_api=igvf_utils_api, upload_file=False, resumed_posting=resumed_posting)
        if upload_file:
            igvf_utils_api.upload_file(
                file_id=curr_post_res, file_path=curr_gs_cloud_link)
        return PostResult.Success(col_header=col_header, uuid=curr_post_res)
    except (requests.exceptions.HTTPError, PortalConflictError) as e:
        return PostResult.Failure(col_header=col_header, error=e)


def post_all_atac_alignment_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, curr_ctrl_access: bool, curr_derived_from: list, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/', resumed_posting: bool = False) -> list[PostResult]:
    """Post all ATAC alignment data to the portal (alignment, bam index files, and QC).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        curr_ctrl_access (bool): The controlled access status (computed based on SeqFile status)
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        award (str): The data submitter lab's award
        igvf_utils_api (_type_): IGVF utils api
        upload_file (bool): Whether to upload the file to the portal
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        list[PostResult]: [PostResult(
            col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    alignment_file_col = 'atac_bam'
    bam_index_col = 'atac_bam_index'
    qc_prefix = 'alignment'

    curr_post_summary = []
    # Post alignment bam
    alignment_result = post_single_alignment_file(terra_data_record=terra_data_record,
                                                  col_header=alignment_file_col,
                                                  curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                      'alignment_file'][alignment_file_col]['file_format'],
                                                  curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                      'alignment_file'][alignment_file_col]['description'],
                                                  lab=lab,
                                                  award=award,
                                                  curr_ctrl_access=curr_ctrl_access,
                                                  curr_derived_from=curr_derived_from,
                                                  igvf_utils_api=igvf_utils_api,
                                                  upload_file=upload_file,
                                                  resumed_posting=resumed_posting)
    curr_post_summary.append(alignment_result)

    # Get the alignment file post result. Return errors early if the post failed
    if alignment_result.error is not None:
        curr_post_summary.extend(midway_failing_chained_posting(
            will_fail_col_headers=[bam_index_col, f'{qc_prefix}_qc']))
        return curr_post_summary

    # Get the posted alignment UUID
    posted_alignment_uuid = alignment_result.accession
    # Post bam index .bam.bai files
    bam_index_result = post_single_index_file(terra_data_record=terra_data_record,
                                              col_header=bam_index_col,
                                              curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                  'index_file'][bam_index_col]['file_format'],
                                              curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                  'index_file'][bam_index_col]['description'],
                                              lab=lab,
                                              award=award,
                                              curr_ctrl_access=curr_ctrl_access,
                                              curr_derived_from=[
                                                  posted_alignment_uuid],
                                              igvf_utils_api=igvf_utils_api,
                                              upload_file=upload_file,
                                              resumed_posting=resumed_posting)
    curr_post_summary.append(bam_index_result)

    # Submit Alignment QC
    curr_workflow_config = parse_workflow_uuids_from_gs_path(
        gs_path=terra_data_record[alignment_file_col])
    # Generate QC object aliases
    qc_obj_aliases = mk_qc_obj_aliases(curr_workflow_config=curr_workflow_config,
                                       analysis_set_acc=terra_data_record['analysis_set_acc'],
                                       qc_prefix=qc_prefix,
                                       lab=lab)
    # Post QC metrics to the portal
    alignment_metric_post_result = post_single_qc_metric_to_portal(terra_data_record=terra_data_record,
                                                                   qc_data_info=TERRA_QC_OUTPUTS['atacseq'][qc_prefix],
                                                                   qc_of_file_accs=[
                                                                       posted_alignment_uuid],
                                                                   lab=lab,
                                                                   award=award,
                                                                   analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                                                       'atac'],
                                                                   qc_prefix=qc_prefix,
                                                                   qc_obj_aliases=qc_obj_aliases,
                                                                   igvf_utils_api=igvf_utils_api,
                                                                   output_root_dir=output_root_dir,
                                                                   resumed_posting=resumed_posting
                                                                   )
    # Add QC post report
    curr_post_summary.append(alignment_metric_post_result)
    # Return all post res
    return curr_post_summary


def post_all_atac_fragment_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, curr_assembly: str, curr_derived_from: list, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/', resumed_posting: bool = False) -> list[PostResult]:
    """Post all ATAC fragment data to the portal (fragment files, index files, and QC).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The controlled access status (computed based on SeqFile status)
        curr_assembly (str): The genome assembly
        curr_derived_from (list): The sequence file and seqspec accessions (for derived_from)
        curr_file_specification (list): The file format specifications
        igvf_api (_type_): IGVF python client api
        igvf_utils_api (_type_): IGVF utils api
        upload_file (bool): Whether to upload a file or not
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        list[PostResult]: [PostResult(
            col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    fragment_file_col = 'atac_fragments'
    fragment_index_file_col = 'atac_fragments_index'
    qc_prefix = 'fragment'

    curr_post_summary = []
    # Post fragment file
    fragment_file_result = post_single_tabular_file(terra_data_record=terra_data_record,
                                                    col_header=fragment_file_col,
                                                    curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                        'tabular_file'][fragment_file_col]['file_format'],
                                                    curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                        'tabular_file'][fragment_file_col]['description'],
                                                    lab=lab,
                                                    award=award,
                                                    curr_assembly=curr_assembly,
                                                    curr_ctrl_access=False,
                                                    curr_derived_from=curr_derived_from,
                                                    curr_file_specification=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                        'tabular_file'][fragment_file_col]['file_format_specifications'],
                                                    igvf_utils_api=igvf_utils_api,
                                                    upload_file=upload_file,
                                                    resumed_posting=resumed_posting)
    curr_post_summary.append(fragment_file_result)

    # Check if the fragment file post result has an error
    if fragment_file_result.error is not None:
        curr_post_summary.extend(midway_failing_chained_posting(
            will_fail_col_headers=[fragment_index_file_col, f'{qc_prefix}_qc']))
        return curr_post_summary

    # Get the posted alignment UUID
    posted_fragment_uuid = fragment_file_result.accession

    # Post fragment index file

    fragment_index_file_result = post_single_index_file(terra_data_record=terra_data_record,
                                                        col_header=fragment_index_file_col,
                                                        curr_file_format=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                            'index_file'][fragment_index_file_col]['file_format'],
                                                        curr_description=TERRA_OUTPUT_TABLE_COLUMN_TYPES[
                                                            'index_file'][fragment_index_file_col]['description'],
                                                        lab=lab,
                                                        award=award,
                                                        curr_ctrl_access=False,
                                                        curr_derived_from=[
                                                            posted_fragment_uuid],
                                                        igvf_utils_api=igvf_utils_api,
                                                        upload_file=upload_file,
                                                        resumed_posting=resumed_posting
                                                        )
    curr_post_summary.append(fragment_index_file_result)

    # Post Fragment QC
    curr_workflow_config = parse_workflow_uuids_from_gs_path(
        gs_path=terra_data_record['atac_fragments'])
    # Generate QC object aliases
    qc_obj_aliases = mk_qc_obj_aliases(curr_workflow_config=curr_workflow_config,
                                       analysis_set_acc=terra_data_record['analysis_set_acc'],
                                       qc_prefix=qc_prefix,
                                       lab=lab)
    fragment_metric_post_result = post_single_qc_metric_to_portal(terra_data_record=terra_data_record,
                                                                  qc_data_info=TERRA_QC_OUTPUTS['atacseq'][qc_prefix],
                                                                  qc_of_file_accs=[
                                                                      posted_fragment_uuid],
                                                                  lab=lab,
                                                                  award=award,
                                                                  analysis_step_version=ANALYSIS_STEP_VERSIONS_BY_ASSAY_TYPES[
                                                                      'atac'],
                                                                  qc_prefix=qc_prefix,
                                                                  qc_obj_aliases=qc_obj_aliases,
                                                                  igvf_utils_api=igvf_utils_api,
                                                                  output_root_dir=output_root_dir,
                                                                  resumed_posting=resumed_posting
                                                                  )
    # Add QC post report
    curr_post_summary.append(fragment_metric_post_result)
    # Return results
    return curr_post_summary


def execute_partial_function(func):
    """Helper function to execute a partial function."""
    return func()


def post_all_atac_data_to_portal(terra_data_record: pd.Series, lab: str, award: str, igvf_api, igvf_utils_api, upload_file: bool, output_root_dir: str = '/igvf/data/', resumed_posting: bool = False) -> list[PostResult]:
    """Post all ATAC data to the portal (alignment, fragment, and index files).

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        lab (str): The data submitter lab
        award (str): The data submitter lab's award
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        list[PostResult]: list[PostResult], e.g., [PostResult(col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)]
    """
    # Current assembly
    curr_assembly = get_genome_assembly(taxa=terra_data_record['taxa'])

    # Get file controlled access level and sequence files
    curr_ctrl_access, curr_seqfile_accs = get_seqfile_access_lvl_and_accessions(
        terra_data_record=terra_data_record, assay_type='atac', igvf_api=igvf_api)

    # Get seqspec file accessions
    curr_seqspec_accs = parse_igvf_accessions_from_urls(
        igvf_file_urls=[terra_data_record['atac_seqspec_urls']])

    # Derived from will be a combination of seqfile, seqspec accessions, and barcode replacement file
    # Beware that some assays will have ATAC barcodes and Read2 concatenated into one fastq
    all_derived_from = list(set(chain(curr_seqfile_accs, curr_seqspec_accs)))

    # Multiprocessing tasks for posting ATAC data
    tasks = [
        partial(post_all_atac_alignment_data_to_portal,
                terra_data_record=terra_data_record,
                lab=lab,
                award=award,
                curr_ctrl_access=curr_ctrl_access,
                curr_derived_from=all_derived_from,
                igvf_utils_api=igvf_utils_api,
                upload_file=upload_file,
                output_root_dir=output_root_dir,
                resumed_posting=resumed_posting),
        partial(post_all_atac_fragment_data_to_portal,
                terra_data_record=terra_data_record,
                lab=lab,
                award=award,
                curr_assembly=curr_assembly,
                curr_derived_from=all_derived_from,
                igvf_utils_api=igvf_utils_api,
                upload_file=upload_file,
                output_root_dir=output_root_dir,
                resumed_posting=resumed_posting)
    ]

    # Use multiprocessing to execute tasks
    atac_post_results = []
    with multiprocessing.Pool() as pool:
        results = pool.map(execute_partial_function, tasks)
        for result in results:
            atac_post_results.extend(result)
    return atac_post_results


# Full data posting from a single pipeline run
def post_all_data_from_one_run(terra_data_record: pd.Series, igvf_api, igvf_utils_api, upload_file: bool, config_file_path: str, doc_aliases: list, pipeline_data_lab: str = POSTING_LAB, pipeline_data_award: str = POSTING_AWARD, output_root_dir: str = '/igvf/data/', resumed_posting: bool = False) -> RunResult:
    """Post all single cell uniform pipeline output to the portal from Terra.

    Args:
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        config_file_path (str): The path to the pipeline run config file
        doc_aliases (list): The aliases for the config file
        pipeline_data_lab (str): The data submitter lab, default IGVF pipeline processing lab
        pipeline_data_award (str): The data submitter lab's award, default IGVF DACC award
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        RunResult: RunResult(analysis_set_acc=IGVFxxx, post_results=[PostResult(
            col_header='rna_kb_output_folder_tar_gz', accession='TSTFI00286528', error=None)], finish_time=timestamp)
    """
    logging.info(
        f'Posting data from {terra_data_record["analysis_set_acc"]} to the portal.')
    # Check if the pipeline is RNAseq-only, ATACseq-only, or multiome-seq
    pipeline_run_output = check_single_or_multiome(
        terra_data_record=terra_data_record)
    # Posting
    post_results = []
    # If RNAseq results are present, post them
    if 'RNAseq' in pipeline_run_output:
        post_results += post_all_rna_data_to_portal(terra_data_record=terra_data_record,
                                                    lab=pipeline_data_lab,
                                                    award=pipeline_data_award,
                                                    igvf_utils_api=igvf_utils_api,
                                                    upload_file=upload_file,
                                                    output_root_dir=output_root_dir,
                                                    resumed_posting=resumed_posting
                                                    )
    # If ATACseq results are present, post them
    if 'ATACseq' in pipeline_run_output:
        post_results += post_all_atac_data_to_portal(terra_data_record=terra_data_record,
                                                     lab=pipeline_data_lab,
                                                     award=pipeline_data_award,
                                                     igvf_api=igvf_api,
                                                     igvf_utils_api=igvf_utils_api,
                                                     upload_file=upload_file,
                                                     output_root_dir=output_root_dir,
                                                     resumed_posting=resumed_posting
                                                     )

    # Add posting run config document and patch analysis set with doc and status completed
    post_results += post_single_document_to_portal(terra_data_record=terra_data_record,
                                                   lab=pipeline_data_lab,
                                                   award=pipeline_data_award,
                                                   config_file_path=config_file_path,
                                                   doc_aliases=doc_aliases,
                                                   igvf_utils_api=igvf_utils_api,
                                                   resumed_posting=resumed_posting)

    # Return full post res
    return RunResult(
        analysis_set_acc=terra_data_record['analysis_set_acc'],
        post_results=post_results,
        finish_time=str(datetime.datetime.now()))


# Post the entire batch of pipeline outputs to the portal
def post_all_successful_runs(full_terra_data_table: pd.DataFrame, igvf_api, igvf_utils_api, upload_file: bool, config_file_collection: dict, output_root_dir: str = '/igvf/data/', resumed_posting: bool = False) -> list[PostResult]:
    """Post all successful runs from a Terra data table to the portal.

    Args:
        full_terra_data_table (pd.DataFrame): The Complete Terra data table
        igvf_api (_type_): IGVF utils api
        igvf_utils_api (_type_): IGVF python client api
        upload_file (bool): Whether to upload the file to the portal
        config_file_collection (dict): A dictionary of config file paths and doc aliases for each pipeline run
        output_root_dir (str, optional): The output directory to download the QC files to. Defaults to '/igvf/data/'.
        resumed_posting (bool, optional): Whether to patch an existing post. Defaults to False. See `single_post_to_portal` for more details.

    Returns:
        list[PostResult]: PostResult objects
    """
    midpt_post_anasets_file = os.path.join(
        output_root_dir, 'midpoint_posted_analysis_sets.txt')
    with open(midpt_post_anasets_file, 'a+') as f:
        results = []
        for _, curr_pipeline in full_terra_data_table.iterrows():
            curr_config_file_info = config_file_collection[curr_pipeline['analysis_set_acc']]
            pipeline_post_res = post_all_data_from_one_run(
                terra_data_record=curr_pipeline,
                igvf_api=igvf_api,
                igvf_utils_api=igvf_utils_api,
                upload_file=upload_file,
                config_file_path=curr_config_file_info.download_path,
                doc_aliases=curr_config_file_info.doc_aliases,
                output_root_dir=output_root_dir,
                resumed_posting=resumed_posting
            )
            # Once a full run is posted, write the AnaSet accession to a file
            # In case of interruption, we can skip those already posted
            f.write(f"{pipeline_post_res.analysis_set_acc}\n")
            f.flush()
            results.append(pipeline_post_res)
            time.sleep(0.1)  # Sleep to avoid hitting API rate limits
    return results


# Functions to summarize POST status and add it to the output data table
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


def save_pipeline_postres_tables(pipeline_postres_table: pd.DataFrame, updated_full_data_table: pd.DataFrame, output_root_dir: str):
    """Save the pipeline input table to a TSV file.

    Args:
        pipeline_input_table (pd.DataFrame): The pipeline input table
        output_dir (str): The output directory
    """
    curr_datetime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    # make output dir if it does not exist
    if os.path.exists(output_root_dir) is False:
        os.makedirs(output_root_dir)

    # The original table has an extra index column from merging
    updated_full_data_table.to_csv(os.path.join(
        output_root_dir, f'full_terra_data_table_with_post_status_{curr_datetime}.tsv'), sep='\t', index=False)

    # Change post result detail table index to match that of updated full data table for easy Terra upload
    postres_table_index_name = [
        col for col in updated_full_data_table.columns if col.startswith('entity:')][0]
    pipeline_postres_table.index.name = postres_table_index_name.replace(
        '_id', '_postres_id')
    # Output post result detail table in Terra format
    pipeline_postres_table.to_csv(os.path.join(
        output_root_dir, f'single-cell_uniform_pipeline_postres_detail_{curr_datetime}.tsv'), sep='\t')
