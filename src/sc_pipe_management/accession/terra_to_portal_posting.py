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

import sc_pipe_management.accession.igvf_payloads as igvf_payloads
import sc_pipe_management.accession.parse_terra_metadata as terra_parse
import sc_pipe_management.accession.constants as const


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


class PortalPostError(Exception):
    """Custom exception for when a POST to the portal fails."""
    pass


class IGVFPostService:
    """Service class to handle posting data to the portal."""

    def __init__(self, igvf_utils_api, data_obj_payload: igvf_payloads.Payload, upload_file: bool = False, resumed_posting: bool = False):
        """
        Initialize the IGVFPostService.
        """
        self.igvf_utils_api = igvf_utils_api
        self.data_obj_payload = data_obj_payload
        self.upload_file = upload_file
        self.resumed_posting = resumed_posting

    def _get_conflict_object_id(self) -> str | None:
        """Get the ID of a conflicting object, if it exists by aliases or MD5."""
        # Check by alias first (preferred method)
        if self.data_obj_payload.aliases:
            existing_file = self.igvf_utils_api.get(
                f"aliases:{self.data_obj_payload.aliases[0]}")
            return existing_file.get('@id') if existing_file else None

        # Fall back to MD5 check
        if self.data_obj_payload.md5sum:
            existing_file = self.igvf_utils_api.get(
                f"md5:{self.data_obj_payload.md5sum}")
            return existing_file.get('@id') if existing_file else None

        return None

    def _single_post_to_portal(self) -> PostResult:
        """POST a single data object to the portal, after checking for conflicts by MD5 or aliases."""
        try:
            payload_dict = self.data_obj_payload.get_payload()
            _schema_property = self.igvf_utils_api.get_profile_from_payload(
                payload_dict).properties
            # First check payload MD5 and aliases to see if the object already exists
            existing_object_id = self._get_conflict_object_id()
            # If conflict exists and not expecting it to, raise an error
            if (existing_object_id is not None) and (not self.resumed_posting):
                raise PortalConflictError(
                    f'{payload_dict.get("_profile").capitalize().replace("_", " ")} failed to accession due to conflict with an existing object: {existing_object_id}.')
            # If conflict exists and patching is expected, return the existing object ID
            if (existing_object_id is not None) and (self.resumed_posting):
                return existing_object_id
            # Post new object to portal if no conflict is expected, metadata only
            stdout = self.igvf_utils_api.post(
                payload_dict, upload_file=False, return_original_status_code=True, truncate_long_strings_in_payload_log=True)
            # Get the newly generated UUID
            new_uuid_generated = stdout[0]['uuid']
            if self.upload_file:
                self.igvf_utils_api.upload_file(
                    file_id=new_uuid_generated, file_path=self.data_obj_payload.submitted_file_name)
            return PostResult.Success(
                col_header=self.data_obj_payload.terra_output_name, uuid=new_uuid_generated)
        except (requests.exceptions.HTTPError, PortalConflictError) as e:
            # Handle errors
            return PostResult.Failure(
                col_header=self.data_obj_payload.terra_output_name, error=e)

    def _single_patch_to_portal(self, analysis_set_acc: str) -> PostResult:
        """PATCH a single data object to the portal, after checking for conflicts by MD5 or aliases."""
        try:
            self.igvf_utils_api.patch(self.data_obj_payload)
            PostResult.Success(
                col_header=self.data_obj_payload['_profile'], accession=analysis_set_acc)
        except (requests.exceptions.HTTPError, PortalConflictError) as e:
            return PostResult.Failure(
                col_header=self.data_obj_payload.terra_output_name, error=e)

    def get_chained_new_uuid_generated(self, lead_chain_post_results: list[PostResult]) -> list[PostResult]:
        """Post data to the portal, handling chained objects."""
        # Get a list of new UUIDs generated from the first order of post results
        new_uuid_generated = [
            res.uuid for res in lead_chain_post_results if res.error is None]

        # If Post failed (returned an empty list), return failure results
        if not new_uuid_generated:
            lead_chain_post_results.append(PostResult.Failure(
                col_header=self.data_obj_payload.terra_output_name,
                error='Cannot post QC metrics because the lead object failed to post.'))
            return lead_chain_post_results


class IGVFAccessioning:
    """Class to handle IGVF accessioning operations."""

    def __init__(self, igvf_utils_api, igvf_client_api, terra_metadata: terra_parse.TerraOutputMetadata, pipeline_params_info: igvf_payloads.PipelineParamsInfo, upload_file: bool = False, resumed_posting: bool = False, root_output_dir: str = '/igvf/data/'):
        self.igvf_utils_api = igvf_utils_api
        self.igvf_client_api = igvf_client_api
        self.upload_file = upload_file
        self.resumed_posting = resumed_posting
        self.terra_metadata = terra_metadata
        self.terra_data_record = terra_metadata.terra_data_record
        self.pipeline_params_info = pipeline_params_info
        self.root_output_dir = root_output_dir

    def _post_qc_metrics(self, file_post_res: list[PostResult], qc_info_map: const.QCInfoMap, qc_post_res_name: str, qc_prefix: str) -> PostResult:
        """Post QC metrics to the portal."""
        # Get a list of files that the QC metrics are for
        posted_files_uuids = [
            res.uuid for res in file_post_res if res.error is None]

        # If no files were posted successfully, return a failure result instead of posting QC metrics
        if not posted_files_uuids:
            return PostResult.Failure(
                col_header=qc_post_res_name,
                error='Cannot post QC metrics because the qc_of object(s) failed to post.'
            )

        # Get QC metrics payloads
        qc_payload = igvf_payloads.QCMetricsPayload(
            terra_metadata=self.terra_metadata,
            qc_info_map=qc_info_map,
            qc_prefix=qc_prefix,
            qc_of=posted_files_uuids,
            igvf_api=self.igvf_utils_api,
            root_output_dir=self.root_output_dir
        )

        # Post QC metrics (no file to upload)
        curr_qc_post_mthd = IGVFPostService(
            igvf_utils_api=self.igvf_utils_api,
            data_obj_payload=qc_payload,
            upload_file=False,
            resumed_posting=self.resumed_posting
        )

        return curr_qc_post_mthd._single_post_to_portal()

    def _post_matrix_files(self) -> list[PostResult]:
        """Post matrix files to the portal."""
        post_results = []

        # Get matrix file payloads
        matrix_payloads = []
        for terra_output_name in list(const.MATRIX_FILETYPES.keys()):
            curr_matrix_payload = igvf_payloads.MatrixFilePayload(
                terra_metadata=self.terra_metadata,
                terra_output_name=terra_output_name
            )
            matrix_payloads.append(curr_matrix_payload)

        # Post matrix files
        for matrix_payload in matrix_payloads:
            curr_mtx_post_mthd = IGVFPostService(
                igvf_utils_api=self.igvf_utils_api,
                data_obj_payload=matrix_payload,
                upload_file=self.upload_file,
                resumed_posting=self.resumed_posting
            )
            post_results.append(curr_mtx_post_mthd._single_post_to_portal())
        return post_results

    def _post_index_file(self, derived_from: list[PostResult], terra_output_name: str) -> PostResult:
        """Post index file for RNAseq to the portal."""
        derived_from_uuids = [
            res.uuid for res in derived_from if res.error is None]

        if not derived_from_uuids:
            return PostResult.Failure(
                col_header=terra_output_name,
                error='Cannot post index file because the derived_from object(s) failed to post.'
            )

        # Get index file payload
        index_file_payload = igvf_payloads.IndexFilePayload(
            terra_metadata=self.terra_metadata,
            derived_from=derived_from_uuids,
            terra_output_name=terra_output_name,
            igvf_api=self.igvf_client_api
        )

        # Post index file
        curr_index_post_mthd = IGVFPostService(
            igvf_utils_api=self.igvf_utils_api,
            data_obj_payload=index_file_payload,
            upload_file=self.upload_file,
            resumed_posting=self.resumed_posting
        )
        return curr_index_post_mthd._single_post_to_portal()

    def _post_alignment_file(self) -> PostResult:
        # Get alignment file payload
        alignment_payload = igvf_payloads.AlignmentFilePayload(
            terra_metadata=self.terra_metadata,
            igvf_api=self.igvf_client_api
        )

        # Post alignment file
        curr_alignment_post_mthd = IGVFPostService(
            igvf_utils_api=self.igvf_utils_api,
            data_obj_payload=alignment_payload,
            upload_file=self.upload_file,
            resumed_posting=self.resumed_posting
        )

        return curr_alignment_post_mthd._single_post_to_portal()

    def _post_fragment_file(self) -> PostResult:
        """Post ATACseq fragment file to the portal."""
        # Get fragment file payload
        fragment_payload = igvf_payloads.FragmentFilePayload(
            terra_metadata=self.terra_metadata,
            igvf_api=self.igvf_client_api
        )

        # Post fragment file
        curr_fragment_post_mthd = IGVFPostService(
            igvf_utils_api=self.igvf_utils_api,
            data_obj_payload=fragment_payload,
            upload_file=self.upload_file,
            resumed_posting=self.resumed_posting
        )

        return curr_fragment_post_mthd._single_post_to_portal()

    def post_all_rnaseq_output(self) -> list[PostResult]:
        """Post RNAseq output data to the portal."""
        post_results = []

        # Post matrix files
        post_results.extend(self._post_matrix_files())

        # Post QC metrics
        rna_qc_prefix = 'gene_count'
        qc_post_result = self._post_qc_metrics(
            file_post_res=post_results,
            qc_info_map=const.TERRA_QC_OUTPUTS['rnaseq'][rna_qc_prefix],
            qc_post_res_name='RNAseq QC Metrics',
            qc_prefix=rna_qc_prefix)
        post_results.append(qc_post_result)

        return post_results

    def post_all_atac_alignment_output(self) -> list[PostResult]:
        """Post ATACseq alignment output data to the portal."""
        post_results = []

        alignment_file_post_result = self._post_alignment_file()
        post_results.append(alignment_file_post_result)

        # Post Index files
        index_file_post_result = self._post_index_file(
            derived_from=[alignment_file_post_result], terra_output_name='atac_bam_index')
        post_results.append(index_file_post_result)

        # Post QC metrics
        alignment_qc_prefix = 'alignment'
        qc_post_result = self._post_qc_metrics(
            file_post_res=[alignment_file_post_result],
            qc_info_map=const.TERRA_QC_OUTPUTS['atacseq'][alignment_qc_prefix],
            qc_post_res_name='ATACseq Alignment QC Metrics',
            qc_prefix=alignment_qc_prefix)
        post_results.append(qc_post_result)

        return post_results

    def post_all_atac_fragment_output(self) -> list[PostResult]:
        """Post ATACseq fragment output data to the portal."""
        post_results = []

        # Post fragment files

        fragment_file_post_result = self._post_fragment_file()
        post_results.append(fragment_file_post_result)

        # Post Index files
        index_file_post_result = self._post_index_file(
            derived_from=[fragment_file_post_result], terra_output_name='atac_fragments_index')
        post_results.append(index_file_post_result)

        # Post QC metrics
        qc_post_result = self._post_qc_metrics(
            file_post_res=[fragment_file_post_result], qc_post_res_name='ATACseq Fragment QC Metrics')
        post_results.append(qc_post_result)

        return post_results

    def post_document(self) -> list[PostResult]:
        """Post workflow configuration document to the portal."""
        # Post workflow configuration document
        doc_payload = igvf_payloads.DocumentPayload(
            terra_metadata=self.terra_metadata,
            pipeline_params_info=self.pipeline_params_info,
            igvf_api=self.igvf_client_api,
        )

        # Post the document to the portal
        doc_post_mthd = IGVFPostService(
            igvf_utils_api=self.igvf_utils_api,
            data_obj_payload=doc_payload,
            upload_file=False,
            resumed_posting=self.resumed_posting
        )

        return doc_post_mthd.post_to_portal()

    def patch_analysis_set(self, document_uuid: str) -> PostResult:
        """Patch the analysis set with the new data."""
        # Patch the analysis set with the new data
        patch_payload = igvf_payloads.AnalysisSetPatchingPayload(
            terra_metadata=self.terra_metadata,
            input_params_doc_uuid=document_uuid,
            igvf_api=self.igvf_client_api
        )._get_patch_payload()

        # If no patch payload is generated, return a failure result
        if patch_payload is None:
            return PostResult.Failure(
                col_header='Analysis Set Patch',
                error='No patch payload generated, possibly due to no changes in the analysis set.'
            )

        patch_mthd = IGVFPostService(
            igvf_utils_api=self.igvf_utils_api,
            data_obj_payload=patch_payload,
            upload_file=False,
            resumed_posting=self.resumed_posting
        )

        return patch_mthd._single_patch_to_portal(analysis_set_acc=self.terra_metadata.anaset_accession)


def post_single_pipeline_run(terra_data_record: pd.Series,
                             pipeline_params_info: igvf_payloads.PipelineParamsInfo,
                             igvf_client_api,
                             igvf_utils_api,
                             output_root_dir: str = '/igvf/data/',
                             upload_file: bool = False,
                             resumed_posting: bool = False
                             ) -> list[PostResult]:
    """Post a single pipeline run to the portal.

    Args:
        terra_data_record (pd.Series): Single row of the Terra data table containing metadata for a pipeline run.
        pipeline_params_info (igvf_payloads.PipelineParamsInfo): Information about the pipeline parameters.
        igvf_client_api (_type_): IGVF client API instance.
        igvf_utils_api (_type_): IGVF utils API instance.
        output_root_dir (str, optional): Root directory for output files. Defaults to '/igvf/data/'.
        upload_file (bool, optional): Whether to upload files. Defaults to False.
        resumed_posting (bool, optional): Whether to resume posting. Defaults to False.

    Returns:
        list[PostResult]: List of PostResults.
    """
    # Terra metadata parsing
    terra_metadata = terra_parse.TerraOutputMetadata(
        terra_data_record=terra_data_record,
        igvf_client_api=igvf_client_api
    )
    # Initialize IGVFAccessioning
    igvf_accessioning = IGVFAccessioning(
        igvf_utils_api=igvf_utils_api,
        igvf_client_api=igvf_client_api,
        terra_metadata=terra_metadata,
        pipeline_params_info=pipeline_params_info,
        upload_file=upload_file,
        resumed_posting=resumed_posting,
        root_output_dir=output_root_dir
    )
    # Post all relevant results
    assays_to_post = terra_metadata.multiome_types

    post_results = []
    # Post RNAseq output if applicable
    if 'RNAseq' in assays_to_post:
        post_results += igvf_accessioning.post_all_rnaseq_output()

    # Post ATACseq alignment output if applicable
    if 'ATACseq' in assays_to_post:
        post_results += igvf_accessioning.post_all_atac_alignment_output()
        post_results += igvf_accessioning.post_all_atac_fragment_output()

    # Post input pipe params document
    document_post_result = igvf_accessioning.post_document()
    post_results.append(document_post_result)

    # Patch the analysis set with the new data
    anaset_patch_result = igvf_accessioning.patch_analysis_set(
        document_uuid=document_post_result.accession)
    post_results.append(anaset_patch_result)

    return post_results


def post_all_pipeline_runs_from_one_submission(terra_data_table: pd.DataFrame,
                                               pipeline_params_info: igvf_payloads.PipelineParamsInfo,
                                               igvf_client_api,
                                               igvf_utils_api,
                                               output_root_dir: str = '/igvf/data/',
                                               upload_file: bool = False,
                                               resumed_posting: bool = False
                                               ) -> list[RunResult]:
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
        list[RunResult]: RunResult objects
    """
    midpt_post_anasets_file = os.path.join(
        output_root_dir, 'midpoint_posted_analysis_sets.txt')
    with open(midpt_post_anasets_file, 'a+') as f:
        run_results = []
        for _, curr_pipeline in terra_data_table.iterrows():
            pipeline_post_res = post_single_pipeline_run(
                terra_data_record=curr_pipeline,
                pipeline_params_info=pipeline_params_info,
                igvf_client_api=igvf_client_api,
                igvf_utils_api=igvf_utils_api,
                output_root_dir=output_root_dir,
                upload_file=upload_file,
                resumed_posting=resumed_posting
            )
            # Once a full run is posted, write the AnaSet accession to a file
            # In case of interruption, we can skip those already posted
            f.write(f"{pipeline_post_res.analysis_set_acc}\n")
            f.flush()
            run_results.append(RunResult(
                analysis_set_acc=curr_pipeline['analysis_set_acc'],
                post_results=pipeline_post_res,
                finish_time=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            ))
            time.sleep(0.1)  # Sleep to avoid hitting API rate limits
    return run_results


# Functions to summarize POST status and add it to the output data table
def summarize_post_status(run_results: list[RunResult]) -> pd.DataFrame:
    """Summarize the POST status of all data from a single pipeline run.

    Args:
        post_status_df (pd.DataFrame): The POST status DataFrame

    Returns:
        pd.DataFrame: The POST status summary DataFrame
    """
    post_status_summary = {}
    for result in run_results:
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
