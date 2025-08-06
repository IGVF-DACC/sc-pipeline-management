import os
import datetime
import subprocess
import pandas as pd
import firecloud.api as fapi
from firecloud.errors import FireCloudServerError
import json
import dataclasses
from typing import Protocol
import logging

import sc_pipe_management.accession.parse_terra_metadata as terra_parse
import sc_pipe_management.accession.constants as const
import sc_pipe_management.igvf_and_terra_api_tools as api_tools


fapi._set_session()


@dataclasses.dataclass(frozen=True)
class QCFileDownloadInfo:
    paths_of_metadata_files: list[str]
    paths_of_attachment_files: list[str]


def _dump_json(input_json: dict, analysis_set_acc: str, output_root_dir: str = './run_config') -> str:
    """Save a JSON object to a file in the specified output directory.

    Args:
        input_json (dict): The JSON object to save (workflow config json)
        analysis_set_acc (str): Analysis set accession, used to name the output file
        output_root_dir (str, optional): Where the file will be saved to. Defaults to './run_config'.

    Returns:
        str: path to the saved JSON file
    """
    file_dir = os.path.join(
        output_root_dir, datetime.datetime.now().strftime("%m%d%Y"))
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
    output_file_path = os.path.join(
        file_dir, f'{analysis_set_acc}_run_config.json')
    with open(output_file_path, 'w') as file:
        json.dump(input_json, file)
    return output_file_path


def _get_file_aliases(col_header: str, lab: str, terra_data_record: pd.Series, terra_uuids: terra_parse.TerraJobUUIDs) -> list:
    """ Generate a file alias based on the output data name, lab, analysis set accession, and Terra UUIDs.

    Args:
        col_header (str): A data table column header, which is the output file name
        lab (str): The data submitter lab
        terra_data_record (pd.Series): One Terra pipeline (i.e., one row in the data table)
        terra_uuids (str): A string of the form GCPbucket_submissionID_workflowID_subworkflowID

    Returns:
        list: A file alias in the form of lab:analysis-set-acc_terra-uuids_col-header_uniform-pipeline
    """
    return [f'{lab.split("/")[-2]}:{terra_data_record["analysis_set_acc"]}_{terra_uuids.aliases()}_{col_header}_uniform-pipeline']


def _download_qc_file_from_gcp(gs_file_path: str, downloaded_dir: str) -> str:
    """Download Google cloud file to a local directory.

    Args:
        gs_file_path (str): GCP file path, e.g., gs://bucket/submissions/final-outputs/...
        downloaded_dir (str): The local directory to download the file to

    Raises:
        Exception: If the download fails, raise an exception with the error message

    Returns:
        str: Downloaded file path
    """
    if not gs_file_path.startswith('gs://'):
        raise Exception(
            f'File path {gs_file_path} is not a valid GCP file path.')
    if not os.path.exists(downloaded_dir):
        os.makedirs(downloaded_dir)
    file_name = os.path.basename(gs_file_path)
    downloaded_file_path = os.path.join(downloaded_dir, file_name)
    # Download
    gcp_download = subprocess.run(
        ['gcloud', 'storage', 'cp', gs_file_path, downloaded_file_path], capture_output=True)
    if gcp_download.returncode != 0:
        err_msg = ','.join([entry for entry in gcp_download.stderr.decode().split(
            '\n') if entry.startswith('ERROR:')])
        raise Exception(f'GCP download error: {err_msg}.')
    return downloaded_file_path


class Payload(Protocol):
    """Protocol for payload classes to implement."""

    def get_payload(self) -> dict[str, any]:
        """Method to get the payload as a dictionary."""
        pass

    @property
    def aliases(self) -> list[str]:
        """Property to get the aliases of the payload."""
        pass

    @property
    def md5sum(self) -> str | None:
        """Property to get the MD5 sum of the payload."""
        pass

    @property
    def terra_output_name(self) -> str:
        """Property to get the Terra output name of the payload."""
        pass

    @property
    def submitted_file_name(self) -> str | None:
        """Property to get the submitted file name of the payload."""
        pass


class MatrixFilePayload:
    """Class to create a matrix file payload for a given Terra output name."""

    def __init__(self, terra_metadata: terra_parse.TerraOutputMetadata, terra_output_name: str):
        # Data object lab and award
        self.lab = const.OUTPUT_SUBMITTER_INFO['lab']
        self.award = const.OUTPUT_SUBMITTER_INFO['award']
        # Terra output metadata for this pipeline run
        self.terra_metadata = terra_metadata
        # The actual pipeline output data Series (one row in the Terra data table)
        self.terra_data_record = self.terra_metadata.terra_data_record
        # The Terra UUIDs for this pipeline run
        self.terra_uuids = self.terra_metadata._parse_workflow_uuids_from_gs_path()
        # The matrix file data class object based on the Terra output name
        self._terra_output_name = terra_output_name
        self.file_obj_metadata = const.MATRIX_FILETYPES[self._terra_output_name]
        # The genome assembly info
        self.assembly = const.GENOME_ASSEMBLY_INFO.get(
            self.terra_metadata.taxa)

    @property
    def aliases(self) -> list[str]:
        """Property to get the aliases of the payload."""
        return _get_file_aliases(col_header=self.terra_output_name,
                                 lab=self.lab,
                                 terra_data_record=self.terra_metadata.terra_data_record,
                                 terra_uuids=self.terra_uuids)

    @property
    def submitted_file_name(self) -> str:
        """Property to get the submitted file name of the payload."""
        return self.terra_data_record[self.terra_output_name]

    @property
    def md5sum(self) -> str:
        """Property to get the MD5 sum of the payload."""
        return api_tools.calculate_gsfile_hex_hash(
            file_path=self.submitted_file_name
        )

    @property
    def terra_output_name(self) -> str:
        """Property to get the Terra output name of the payload."""
        return self._terra_output_name

    def get_payload(self):
        """Get the matrix file payload for the given Terra output name."""
        # Compute derived_from and reference files
        mtx_file_input_files = self.terra_metadata._get_input_file_accs_from_table(
            assay_type=self.file_obj_metadata.assay_type)
        # Make Mtx file payload
        mtx_file_payload = dict(award=self.award,
                                lab=self.lab,
                                analysis_step_version=self.file_obj_metadata.analysis_step_version,
                                aliases=self.aliases,
                                content_type=self.file_obj_metadata.content_type,
                                md5sum=self.md5sum,
                                file_format=self.file_obj_metadata.file_format,
                                derived_from=mtx_file_input_files.get_derived_from(),
                                submitted_file_name=self.submitted_file_name,
                                file_set=self.terra_metadata.anaset_accession,
                                principal_dimension='cell',
                                secondary_dimensions=['gene'],
                                filtered=False,
                                reference_files=mtx_file_input_files.reference_files,
                                description=self.file_obj_metadata.description,
                                file_format_specifications=self.file_obj_metadata.file_format_specifications,
                                _profile='matrix_file')
        return mtx_file_payload


class AlignmentFilePayload:
    """Class to create a tabular file payload for a given Terra output name."""

    def __init__(self, terra_metadata: terra_parse.TerraOutputMetadata, igvf_api):
        # Data object lab and award
        self.lab = const.OUTPUT_SUBMITTER_INFO['lab']
        self.award = const.OUTPUT_SUBMITTER_INFO['award']
        # IGVF client API for data access
        self.igvf_api = igvf_api
        # Data column name
        self._terra_output_name = 'atac_bam'
        # Terra output metadata for this pipeline run
        self.terra_metadata = terra_metadata
        # The actual pipeline output data Series (one row in the Terra data table)
        self.terra_data_record = self.terra_metadata.terra_data_record
        # The Terra UUIDs for this pipeline run
        self.terra_uuids = terra_metadata._parse_workflow_uuids_from_gs_path()
        # The tabular file data class object based on the Terra output name
        self.file_obj_metadata = const.ALIGNMENT_FILETYPES[self._terra_output_name]
        # Input file accessions
        self.input_file_accessions = self.terra_metadata._get_input_file_accs_from_table(
            assay_type=self.file_obj_metadata.assay_type)

    def _get_access_status(self) -> bool:
        """Get file controlled_access status. If any seq data is controlled access, the output data inherit that.

        Args:
            sequence_file_accessions (list): Input sequence file accessions
            igvf_api (_type_): _description_

        Returns:
            bool: True or False
        """
        access_status = 0
        for seq_file_acc in self.input_file_accessions.sequence_files:
            seqfile_obj = self.igvf_api.get_by_id(
                f'/sequence-files/{seq_file_acc}').actual_instance
            access_status += seqfile_obj.controlled_access
        if access_status > 0:
            return True
        else:
            return False

    @property
    def aliases(self) -> list[str]:
        """Property to get the aliases of the payload."""
        return _get_file_aliases(col_header=self.terra_output_name,
                                 lab=self.lab,
                                 terra_data_record=self.terra_metadata.terra_data_record,
                                 terra_uuids=self.terra_uuids)

    @property
    def submitted_file_name(self) -> str:
        """Property to get the submitted file name of the payload."""
        return self.terra_data_record[self.terra_output_name]

    @property
    def md5sum(self) -> str | None:
        """Property to get the MD5 sum of the payload."""
        return api_tools.calculate_gsfile_hex_hash(
            file_path=self.submitted_file_name
        )

    @property
    def terra_output_name(self) -> str:
        """Property to get the Terra output name of the payload."""
        return self._terra_output_name

    def get_payload(self):
        """Get the tabular file payload for the given Terra output name."""
        # Make Alignment file payload
        alignment_file_payload = dict(award=self.award,
                                      lab=self.lab,
                                      analysis_step_version=self.file_obj_metadata.analysis_step_version,
                                      aliases=self.aliases,
                                      controlled_access=self._get_access_status(),
                                      file_format=self.file_obj_metadata.file_format,
                                      content_type=self.file_obj_metadata.content_type,
                                      filtered=False,
                                      derived_from=self.input_file_accessions.get_derived_from(),
                                      file_set=self.terra_data_record['analysis_set_acc'],
                                      md5sum=self.md5sum,
                                      submitted_file_name=self.terra_data_record[self.terra_output_name],
                                      reference_files=self.input_file_accessions.reference_files,
                                      redacted=False,
                                      description=self.file_obj_metadata.description,
                                      _profile='alignment_file'
                                      )
        return alignment_file_payload


class FragmentFilePayload:
    """Class to create a fragment file payload for a given Terra output name."""

    def __init__(self, terra_metadata: terra_parse.TerraOutputMetadata, igvf_api):
        # Data object lab and award
        self.lab = const.OUTPUT_SUBMITTER_INFO['lab']
        self.award = const.OUTPUT_SUBMITTER_INFO['award']
        # IGVF client API for data access
        self.igvf_api = igvf_api
        # Data column name
        self._terra_output_name = 'atac_fragments'
        # Terra output metadata for this pipeline run
        self.terra_metadata = terra_metadata
        # The actual pipeline output data Series (one row in the Terra data table)
        self.terra_data_record = self.terra_metadata.terra_data_record
        # The Terra UUIDs for this pipeline run
        self.terra_uuids = terra_metadata._parse_workflow_uuids_from_gs_path()
        # The tabular file data class object based on the Terra output name
        self.file_obj_metadata = const.TABULAR_FILETYPES[self._terra_output_name]
        # Input file accessions
        self.input_file_accessions = self.terra_metadata._get_input_file_accs_from_table(
            assay_type=self.file_obj_metadata.assay_type)

    @property
    def aliases(self) -> list[str]:
        """Property to get the aliases of the payload."""
        return _get_file_aliases(col_header=self.terra_output_name,
                                 lab=self.lab,
                                 terra_data_record=self.terra_metadata.terra_data_record,
                                 terra_uuids=self.terra_uuids)

    @property
    def submitted_file_name(self) -> str:
        """Property to get the submitted file name of the payload."""
        return self.terra_data_record[self.terra_output_name]

    @property
    def md5sum(self) -> str | None:
        """Property to get the MD5 sum of the payload."""
        return api_tools.calculate_gsfile_hex_hash(
            file_path=self.submitted_file_name
        )

    @property
    def terra_output_name(self) -> str:
        """Property to get the Terra output name of the payload."""
        return self._terra_output_name

    def get_payload(self):
        """Get the tabular file payload for the given Terra output name."""
        # Make Fragment file payload
        fragment_file_payload = dict(award=self.award,
                                     lab=self.lab,
                                     analysis_step_version=self.file_obj_metadata.analysis_step_version,
                                     aliases=self.aliases,
                                     content_type=self.file_obj_metadata.content_type,
                                     controlled_access=False,
                                     file_format=self.file_obj_metadata.file_format,
                                     file_format_type='bed3+',
                                     md5sum=self.md5sum,
                                     derived_from=self.input_file_accessions.get_derived_from(),
                                     submitted_file_name=self.terra_data_record[self.terra_output_name],
                                     file_set=self.terra_data_record['analysis_set_acc'],
                                     description=self.file_obj_metadata.description,
                                     filtered=False,
                                     assembly=const.GENOME_ASSEMBLY_INFO.get(
                                         self.terra_metadata.taxa),
                                     file_format_specifications=self.file_obj_metadata.file_format_specifications,
                                     _profile='tabular_file'
                                     )
        return fragment_file_payload


class IndexFilePayload:
    """Class to create an index file payload for a given Terra output name."""

    def __init__(self, terra_metadata: terra_parse.TerraOutputMetadata, derived_from: list[str], terra_output_name: str, igvf_api):
        """Initialize the IndexFilePayload class."""
        # Data object lab and award
        self.lab = const.OUTPUT_SUBMITTER_INFO['lab']
        self.award = const.OUTPUT_SUBMITTER_INFO['award']
        # IGVF client API for data access
        self.igvf_api = igvf_api
        # Terra output metadata for this pipeline run
        self.terra_metadata = terra_metadata
        # The actual pipeline output data Series (one row in the Terra data table)
        self.terra_data_record = self.terra_metadata.terra_data_record
        # The Terra UUIDs for this pipeline run
        self.terra_uuids = terra_metadata._parse_workflow_uuids_from_gs_path()
        # The tabular file data class object based on the Terra output name
        self._terra_output_name = terra_output_name
        self.file_obj_metadata = const.INDEX_FILETYPES[self._terra_output_name]
        # The derived_from file accession
        self.derived_from = list(set(derived_from))

    @property
    def aliases(self) -> list[str]:
        """Property to get the aliases of the payload."""
        return _get_file_aliases(col_header=self.terra_output_name,
                                 lab=self.lab,
                                 terra_data_record=self.terra_metadata.terra_data_record,
                                 terra_uuids=self.terra_uuids)

    @property
    def submitted_file_name(self) -> str:
        """Property to get the submitted file name of the payload."""
        return self.terra_data_record[self.terra_output_name]

    @property
    def md5sum(self) -> str | None:
        """Property to get the MD5 sum of the payload."""
        return api_tools.calculate_gsfile_hex_hash(
            file_path=self.submitted_file_name
        )

    @property
    def terra_output_name(self) -> str:
        """Property to get the Terra output name of the payload."""
        return self._terra_output_name

    def _get_controlled_access_from_derived_from(self) -> bool:
        """Check if any of the derived_from files are controlled access."""
        for file_uuid in self.derived_from:
            file_obj = self.igvf_api.get_by_id(file_uuid).actual_instance
            if file_obj.controlled_access:
                return True
        return False

    def get_payload(self):
        """Get the tabular file payload for the given Terra output name."""
        # Make Index file payload
        index_file_payload = dict(award=self.award,
                                  lab=self.lab,
                                  analysis_step_version=self.file_obj_metadata.analysis_step_version,
                                  aliases=self.aliases,
                                  content_type='index',
                                  controlled_access=self._get_controlled_access_from_derived_from(),
                                  file_format=self.file_obj_metadata.file_format,
                                  md5sum=self.md5sum,
                                  derived_from=self.derived_from,
                                  submitted_file_name=self.terra_data_record[self.terra_output_name],
                                  file_set=self.terra_data_record['analysis_set_acc'],
                                  description=self.file_obj_metadata.description,
                                  _profile='index_file'
                                  )
        return index_file_payload


class QCMetricsPayload:
    """Class to create a QC metrics payload for a given Terra output name."""

    def __init__(self, terra_metadata: terra_parse.TerraOutputMetadata, qc_info_map: const.QCInfoMap, qc_prefix: str, qc_of: list[str], igvf_api, root_output_dir: str = '/igvf/data/'):
        # Data object lab and award
        self.lab = const.OUTPUT_SUBMITTER_INFO['lab']
        self.award = const.OUTPUT_SUBMITTER_INFO['award']
        # IGVF client API for data access
        self.igvf_api = igvf_api
        # QC info map for this pipeline run
        self.qc_info_map = qc_info_map
        # QC metrics of files
        self.qc_of = qc_of
        # Data column name
        self._terra_output_name = f'{qc_prefix}_metrics'
        # Terra output metadata for this pipeline run
        self.terra_metadata = terra_metadata
        # The actual pipeline output data Series (one row in the Terra data table)
        self.terra_data_record = self.terra_metadata.terra_data_record
        # The Terra UUIDs for this pipeline run
        self.terra_uuids = self.terra_metadata._parse_workflow_uuids_from_gs_path()
        # Output root directory
        self.output_root_dir = root_output_dir

    @staticmethod
    def _read_json_file(json_file_paths: list[str]) -> dict:
        """Read a JSON file and return its content as a dictionary.
        Args:
            json_file_path (str): The path to the JSON file.
        Returns:
            dict: The content of the JSON file as a dictionary.
        """
        json_dict = {}
        for json_file_path in json_file_paths:
            with open(json_file_path, 'r') as file:
                json_dict.update(json.load(file))
        return json_dict

    def _get_qc_files(self) -> QCFileDownloadInfo:
        # Set up download directory for QC files
        curr_qc_file_dir = os.path.join(
            self.output_root_dir, 'qc_metrics', self.terra_data_record['analysis_set_acc'])
        if not os.path.exists(curr_qc_file_dir):
            os.makedirs(curr_qc_file_dir)
        # Download attachment files
        download_attachment_files = []
        if self.qc_info_map.attachment is not None:
            for key, value in self.qc_info_map.attachment.items():
                # Download the file first
                curr_attachment_file = _download_qc_file_from_gcp(
                    gs_file_path=self.terra_data_record[value], downloaded_dir=curr_qc_file_dir)
                download_attachment_files.append(
                    {key: {'path': curr_attachment_file}})
        # Download QC files that will be converted to payloads
        downloaded_metadata_files = []
        if self.qc_info_map.metadata is not None:
            # Downloaded JSON file paths
            for metadata_qc_file_name in self.qc_info_map.metadata:
                downloaded_metadata_files.append(_download_qc_file_from_gcp(
                    gs_file_path=self.terra_data_record[metadata_qc_file_name], downloaded_dir=curr_qc_file_dir))
        return QCFileDownloadInfo(
            paths_of_metadata_files=downloaded_metadata_files,
            paths_of_attachment_files=download_attachment_files
        )

    @property
    def aliases(self) -> list[str]:
        return [f'{self.lab.split("/")[-2]}:{self.terra_data_record["analysis_set_acc"]}_{self.terra_uuids.aliases()}_{self._terra_output_name}_uniform-pipeline']

    @property
    def submitted_file_name(self) -> str:
        """Property to get the submitted file name of the payload."""
        return None

    @property
    def md5sum(self):
        """Property to get the MD5 sum of the payload."""
        return None

    @property
    def terra_output_name(self) -> str:
        """Property to get the Terra output name of the payload."""
        return self._terra_output_name

    def get_payload(self):
        curr_qc_file_dir = os.path.join(
            self.output_root_dir, 'qc_metrics', self.terra_data_record['analysis_set_acc'])
        if not os.path.exists(curr_qc_file_dir):
            os.makedirs(curr_qc_file_dir)
        # QC object payload
        qc_payload = dict(lab=self.lab,
                          award=self.award,
                          aliases=self.aliases,
                          analysis_step_version=self.qc_info_map.analysis_step_version,
                          description=self.qc_info_map.description,
                          quality_metric_of=self.qc_of,
                          _profile=self.qc_info_map.object_type
                          )
        # Downloaded QC files
        qc_file_download_info = self._get_qc_files()
        # Add attachments (key is portal property name, value is table colname)
        # Do this first because some QC obj only has attachments
        if qc_file_download_info.paths_of_attachment_files:
            for item in qc_file_download_info.paths_of_attachment_files:
                # Download the file first
                qc_payload.update(item)
        # Add on other metadata
        if qc_file_download_info.paths_of_metadata_files:
            qc_json_parsed = self._read_json_file(
                json_file_paths=qc_file_download_info.paths_of_metadata_files)
            # Update base qc payload
            for key, value in self.qc_info_map.metadata_map.items():
                qc_payload.update({value: qc_json_parsed[key]})
        return qc_payload


class PipelineParamsInfo:
    def __init__(self, terra_datatable: pd.DataFrame, igvf_client_api, terra_namespace: str = 'DACC_ANVIL', terra_workspace: str = 'IGVF Single-Cell Data Processing', output_root_dir: str = '/igvf/data/'):
        self.terra_namespace = terra_namespace
        self.terra_workspace = terra_workspace
        self.output_root_dir = output_root_dir
        # Initialize with Terra metadata class object
        self.terra_datatable = terra_datatable
        # IGVF client API for data access
        self.igvf_client_api = igvf_client_api

    def _get_single_input_params(self, terra_metadata: terra_parse.TerraOutputMetadata) -> str:
        """Get the workflow input configuration JSON file path for a given submission and workflow ID."""
        terra_uuids = terra_metadata._parse_workflow_uuids_from_gs_path()
        # Firecloud API call to get the workflow metadata
        input_params_request = fapi.get_workflow_metadata(namespace=self.terra_namespace,
                                                          workspace=self.terra_workspace,
                                                          submission_id=terra_uuids.submission_id,
                                                          workflow_id=terra_uuids.workflow_id
                                                          )
        if input_params_request.status_code != 200:
            raise FireCloudServerError(code=input_params_request.status_code,
                                       message=f'Error fetching input params data for {terra_metadata.anaset_accession}.')
        # Parse the response JSON to get the workflow input configuration
        workflow_data_res = input_params_request.json()
        try:
            input_params = workflow_data_res.get('inputs')
            # Remove one of the keys that is not needed in the config
            input_params['igvf_credentials'] = 'Google cloud path to a txt file with credentials to access the IGVF data portal.'
            # Set up output directory for the JSON file
            json_output_dir = os.path.join(
                self.output_root_dir, terra_metadata.anaset_accession)
            # Dumpt the config to a local JSON file
            local_file_path = _dump_json(input_json=input_params,
                                         analysis_set_acc=terra_metadata.anaset_accession,
                                         output_root_dir=json_output_dir)
            return local_file_path
        except KeyError:
            raise FireCloudServerError(code=input_params_request.status_code,
                                       message=f'Error parsing workflow input params for{terra_metadata.anaset_accession}. No inputs found.')

    def get_all_input_params(self) -> dict:
        """Get the workflow input configuration JSON for all submissions and workflow IDs."""
        all_input_params = {}
        try:
            for _, terra_data_record in self.terra_datatable.iterrows():
                terra_metadata = terra_parse.TerraOutputMetadata(
                    terra_data_record, self.igvf_client_api)
                input_params = self._get_single_input_params(terra_metadata)
                all_input_params[terra_metadata.anaset_accession] = input_params
        except FireCloudServerError as e:
            logging.debug(
                f'Error fetching input params for {terra_metadata.anaset_accession}: {e.message}')
        return all_input_params


class DocumentPayload:
    """Class to create a document payload for a given Terra output name."""

    def __init__(self, terra_metadata: terra_parse.TerraOutputMetadata, anaset_input_params_file_paths: dict, igvf_api):
        # Data object lab and award
        self.lab = const.OUTPUT_SUBMITTER_INFO['lab']
        self.award = const.OUTPUT_SUBMITTER_INFO['award']
        # IGVF client API for data access
        self.igvf_api = igvf_api
        # The Terra UUIDs for this pipeline run
        self.terra_uuids = terra_metadata._parse_workflow_uuids_from_gs_path()
        self._terra_output_name = 'pipeline_parameters'
        # Pipeline run input parameters for this pipeline run
        self.input_params_file_path = anaset_input_params_file_paths[
            terra_metadata.anaset_accession]

    def _mk_doc_aliases(self) -> list:
        """Create a list of document aliases for the workflow configuration."""
        return [f'igvf-dacc-processing-pipeline:{self.terra_uuids.input_param_aliases()}_pipeline_config']

    @property
    def aliases(self) -> list[str]:
        """Property to get the aliases of the payload."""
        return self._mk_doc_aliases()

    @property
    def submitted_file_name(self) -> str:
        """Property to get the submitted file name of the payload."""
        return None

    @property
    def md5sum(self):
        """Property to get the MD5 sum of the payload."""
        return None

    @property
    def terra_output_name(self) -> str:
        """Property to get the Terra output name of the payload."""
        return self._terra_output_name

    def get_payload(self) -> dict:
        """Get the document payload for the given Terra output name."""
        # local file path for the workflow configuration JSON
        doc_payload = dict(award=self.award,
                           lab=self.lab,
                           aliases=self.aliases,
                           document_type='pipeline parameters',
                           attachment={'path': self.input_params_file_path},
                           description='Terra workflow configuration for the single-cell pipeline run',
                           _profile='document')
        return doc_payload


class AnalysisSetPatchingPayload:
    """Class to create a patching payload for an analysis set."""

    def __init__(self, terra_metadata: terra_parse.TerraOutputMetadata, input_params_doc_uuid: str, igvf_utils_api):
        self.terra_metadata = terra_metadata
        self.anaset_accession = terra_metadata.anaset_accession
        self.input_params_doc_uuid = input_params_doc_uuid
        self.igvf_utils_api = igvf_utils_api

    def _get_existing_analysis_set_docs(self) -> list:
        """Get existing document UUID linked to the analysis set."""
        analysis_set_obj = self.igvf_utils_api.get(
            f'/analysis-sets/{self.anaset_accession}')
        if not analysis_set_obj.get('documents'):
            return []
        return sorted([doc_uuid.split('/')[-2] for doc_uuid in analysis_set_obj['documents']])

    def get_patch_payload(self) -> dict:
        """Get the patch payload for the analysis set."""
        anaset_document_uuids = self._get_existing_analysis_set_docs()
        # If Analysis set has no documents, create a documents array
        if self.input_params_doc_uuid in anaset_document_uuids:
            return None
        anaset_document_uuids.append(self.input_params_doc_uuid)
        return {
            'documents': anaset_document_uuids,
            self.igvf_utils_api.IGVFID_KEY: f"/analysis-sets/{self.anaset_accession}/",
            'uniform_pipeline_status': 'completed',
            '_profile': 'analysis_set'
        }
