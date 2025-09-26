import json
import argparse
import sys
import os
import logging

# Add the project root directory to Python path
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, src_path)

import sc_pipe_management.wrangler_utils.constant as const
import sc_pipe_management.igvf_and_terra_api_tools as api_tools


# Configure logging to show INFO level messages
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def _get_active_file_objs(file_ids: list[str], igvf_client_api: const.IgvfApiType) -> list[const.FileObjTypes]:
    """Get active file objects from a list of file IDs."""
    if not file_ids:
        return []
    all_file_objs = [
        igvf_client_api.get_by_id(file_id).actual_instance
        for file_id in file_ids]
    if not all_file_objs:
        return []
    return [
        file_obj for file_obj in all_file_objs
        if file_obj.status in ['in progress', 'released']]


def _check_if_controlled_access(analysis_file_obj: const.FileObjTypes, igvf_client_api: const.IgvfApiType) -> bool:
    """Get the raw data access level from the analysis file object."""
    all_input_seqfile_objs = []
    for file_set in analysis_file_obj.input_file_sets:
        if file_set.startswith('/measurement-sets/'):
            file_set_obj = igvf_client_api.get_by_id(file_set).actual_instance
            all_input_seqfile_objs.extend([
                igvf_client_api.get_by_id(file_id).actual_instance
                for file_id in file_set_obj.files
                if file_id.startswith('/sequence-files/')
            ])
    return any([seqfile_obj.controlled_access is True for seqfile_obj in all_input_seqfile_objs])


def _is_sublist_in_nested(posted_reference_files: list[str], reference_list: list[str] | list[list[str]]) -> bool:
    """Check if sublist is contained in main_list (which can be 1-layer or nested)."""
    if all(isinstance(item, str) for item in reference_list):
        return set(posted_reference_files).issubset(set(reference_list))
    else:
        for nested_item in reference_list:
            if isinstance(nested_item, list):
                if set(posted_reference_files).issubset(set(nested_item)):
                    return True
        return False


class BaseFileChecker:
    """Base class for file-specific checks."""

    def __init__(self, analysis_set_acc: str, igvf_client_api: const.IgvfApiType):
        self.analysis_set_acc = analysis_set_acc
        self.igvf_client_api = igvf_client_api

    def _check_basic_file_properties(self, file_obj: const.FileObjTypes, expected_asv: str, controlled_access: bool = None) -> list[str]:
        """Check basic file properties and return all errors found."""
        errors = []

        # Check analysis step version first
        if file_obj.analysis_step_version is None:
            errors.append(
                f'{file_obj.type[0]} {file_obj.accession} has no analysis step version.')
        elif file_obj.analysis_step_version != expected_asv:
            errors.append(
                f'{file_obj.type[0]} {file_obj.accession} has unexpected analysis step version.')

        # Check upload status
        if file_obj.upload_status != 'validated':
            errors.append(
                f'{file_obj.type[0]} {file_obj.accession} has invalid upload status.')

        # Check controlled access if provided
        if controlled_access is not None:
            if controlled_access and file_obj.controlled_access is False:
                errors.append(
                    f'URGENT: {file_obj.type[0]} file {file_obj.accession} should be controlled access True, not False.')
            elif not controlled_access and file_obj.controlled_access is True:
                errors.append(
                    f'URGENT: {file_obj.type[0]} file {file_obj.accession} should be controlled access False, not True.')

        return errors

    def _check_reference_files(self, file_obj: const.FileObjTypes, reference_list: list[str] | list[list[str]]):
        """Check if the reference files are in the expected list."""
        if file_obj.reference_files is None:
            raise ValueError(
                f'{file_obj.type[0]} {file_obj.accession} has no reference files.')
        if not _is_sublist_in_nested(file_obj.reference_files, reference_list):
            raise ValueError(
                f'{file_obj.type[0]} {file_obj.accession} has unexpected reference files: {file_obj.reference_files}. Expected one of: {reference_list}.')

    def _check_file_name(self, file_obj: const.FileObjTypes, expected_file_substring: str):
        """Check if the file name contains the expected substring."""
        if expected_file_substring not in file_obj.submitted_file_name:
            raise ValueError(
                f'{file_obj.type[0]} {file_obj.accession} has unexpected file name: {file_obj.submitted_file_name}.')


class DocumentChecker(BaseFileChecker):
    """Checker for pipeline parameters documents."""

    def check_pipeline_parameters_document(self, pipeline_parameters: list[str]):
        """Check pipeline parameters document."""
        if not pipeline_parameters:
            raise ValueError(
                f'No pipeline parameters document found for analysis set {self.analysis_set_acc}.')

        for doc in pipeline_parameters:
            doc_obj = self.igvf_client_api.get_by_id(doc).actual_instance
            if doc_obj.document_type != 'pipeline parameters':
                raise ValueError(
                    f'Unexpected pipeline parameters document type: {doc_obj.document_type}.')
            if doc_obj.attachment is None:
                raise ValueError(
                    f'Pipeline parameters document {doc} has no attachment.')


class QCMetricsChecker(BaseFileChecker):
    """Checker for quality control metrics."""

    def check_qc_metrics_basic(self, file_obj: const.FileObjTypes):
        """Check basic QC metrics properties."""
        if file_obj.quality_metrics is None:
            raise ValueError(
                f'{file_obj.type[0]} {file_obj.accession} has no quality metrics.')
        if len(file_obj.quality_metrics) != 1:
            raise ValueError(
                f'{file_obj.type[0]} {file_obj.accession} has more than 1 quality metrics.')

        qc_metrics_obj = self.igvf_client_api.get_by_id(
            file_obj.quality_metrics[0]).actual_instance

        if qc_metrics_obj.analysis_step_version is None:
            raise ValueError(
                f'QC metrics object linked to {file_obj.type[0]} {file_obj.accession} has no analysis step version.')

        return qc_metrics_obj

    def check_rna_qc_metrics(self, file_obj: const.FileObjTypes):
        """Check RNA-seq specific QC metrics."""
        qc_metrics_obj = self.check_qc_metrics_basic(file_obj)

        if qc_metrics_obj.type[0] != 'SingleCellRnaSeqQualityMetric':
            raise ValueError(
                f'Unexpected QC metrics object type: {qc_metrics_obj.type} linked to {file_obj.type[0]} {file_obj.accession}.')
        if qc_metrics_obj.rnaseq_kb_info is None:
            raise ValueError(
                f'QC metrics object linked to {file_obj.type[0]} {file_obj.accession} has no `rnaseq_kb_info` attachment.')
        if qc_metrics_obj.n_reads is None:
            raise ValueError(
                f'QC metrics object linked to {file_obj.type[0]} {file_obj.accession} may be missing run metrics such as `n_reads`.')

    def check_atac_qc_metrics(self, file_obj: const.FileObjTypes):
        """Check ATAC-seq specific QC metrics."""
        qc_metrics_obj = self.check_qc_metrics_basic(file_obj)

        if qc_metrics_obj.type[0] != 'SingleCellAtacSeqQualityMetric':
            raise ValueError(
                f'Unexpected QC metrics object type: {qc_metrics_obj.type} linked to {file_obj.type[0]} {file_obj.accession}.')

        # File type specific checks
        if file_obj.type[0] == 'AlignmentFile':
            if qc_metrics_obj.atac_bam_summary_stats is None:
                raise ValueError(
                    f'QC metrics object linked to {file_obj.type[0]} {file_obj.accession} has no `atac_bam_summary_stats` attachment.')
        elif file_obj.type[0] == 'TabularFile':
            if qc_metrics_obj.atac_fragment_summary_stats is None:
                raise ValueError(
                    f'QC metrics object linked to {file_obj.type[0]} {file_obj.accession} has no `atac_fragment_summary_stats` attachment.')
            if qc_metrics_obj.atac_fragments_alignment_stats is None:
                raise ValueError(
                    f'QC metrics object linked to {file_obj.type[0]} {file_obj.accession} has no `atac_fragments_alignment_stats` attachment.')
            if qc_metrics_obj.n_mapped_reads is None:
                raise ValueError(
                    f'QC metrics object linked to {file_obj.type[0]} {file_obj.accession} may be missing run metrics such as `n_mapped_reads`.')


class IndexFileChecker(BaseFileChecker):
    """Checker for index files."""

    def _get_index_file_ids(self, file_obj: const.FileObjTypes) -> list[str]:
        """Extract index file IDs from the file object."""
        return [
            file_id for file_id in file_obj.input_file_for
            if file_id.startswith('/index-files/')
        ]

    def _validate_index_file_count(self, file_obj: const.FileObjTypes, active_index_file_objs: list[const.FileObjTypes]):
        """Validate that exactly one index file is present and active."""
        if not active_index_file_objs:
            raise ValueError(
                f'{file_obj.type[0]} {file_obj.accession} is missing an index file.')
        if len(active_index_file_objs) != 1:
            raise ValueError(
                f'{file_obj.type[0]} {file_obj.accession} has more than 1 active index file linked.')

    def _check_bai_index_file(self, file_obj: const.FileObjTypes, index_file_obj: const.FileObjTypes):
        """Check BAI index file specific properties for AlignmentFile."""
        if index_file_obj.file_format != 'bai':
            raise ValueError(
                f'Index file {index_file_obj.accession} linked to AlignmentFile {file_obj.accession} has unexpected file format {index_file_obj.file_format}.')
        if file_obj.controlled_access != index_file_obj.controlled_access:
            raise ValueError(
                f'URGENT: Index file {index_file_obj.accession} linked to AlignmentFile {file_obj.accession} do not have matching controlled_access status.')
        self._check_file_name(
            file_obj=index_file_obj,
            expected_file_substring=f'{self.analysis_set_acc}.bam.bai')

    def _check_tbi_index_file(self, file_obj: const.FileObjTypes, index_file_obj: const.FileObjTypes):
        """Check TBI index file specific properties for TabularFile."""
        if index_file_obj.file_format != 'tbi':
            raise ValueError(
                f'Index file {index_file_obj.accession} linked to TabularFile {file_obj.accession} has unexpected file format {index_file_obj.file_format}.')
        self._check_file_name(
            file_obj=index_file_obj,
            expected_file_substring=f'{self.analysis_set_acc}.fragments.tsv.gz.tbi')

    def _check_index_file_by_type(self, file_obj: const.FileObjTypes, index_file_obj: const.FileObjTypes):
        """Check index file properties based on the main file type."""
        if file_obj.type[0] == 'AlignmentFile':
            self._check_bai_index_file(file_obj, index_file_obj)
        elif file_obj.type[0] == 'TabularFile':
            self._check_tbi_index_file(file_obj, index_file_obj)

    def check_index_file_presence_and_format(self, file_obj: const.FileObjTypes):
        """Check if index file is present and has correct format."""
        # Only check for files that should have index files
        if file_obj.type[0] not in ['AlignmentFile', 'TabularFile']:
            return

        # Get index file IDs
        index_file_ids = self._get_index_file_ids(file_obj)

        # Get active index file objects
        active_index_file_objs = _get_active_file_objs(
            file_ids=index_file_ids, igvf_client_api=self.igvf_client_api)

        # Validate count
        self._validate_index_file_count(
            file_obj=file_obj, active_index_file_objs=active_index_file_objs)

        # Check specific index file properties
        index_file_obj = active_index_file_objs[0]
        self._check_index_file_by_type(file_obj, index_file_obj)


class RNASeqFileChecker(BaseFileChecker):
    """Checker for RNA-seq files."""

    def check_rna_file_count(self, active_rna_file_objs: list[const.FileObjTypes]):
        """Check if exactly 2 RNA-seq files are present."""
        if len(active_rna_file_objs) != 2:
            raise ValueError(
                f'Analysis set {self.analysis_set_acc} does not have exactly 2 non-deprecated RNA-seq output data.')

    def check_rna_file_format_and_content(self, rna_file_obj: const.FileObjTypes):
        """Check RNA file format and content type."""
        if rna_file_obj.file_format == 'h5ad':
            if rna_file_obj.content_type != 'sparse gene count matrix':
                raise ValueError(
                    f'RNA-seq file {rna_file_obj.accession} has unexpected content type for h5ad format.')
            self._check_file_name(
                rna_file_obj, expected_file_substring=f'{self.analysis_set_acc}.h5ad')
        elif rna_file_obj.file_format == 'tar':
            if rna_file_obj.content_type != 'kallisto single cell RNAseq output':
                raise ValueError(
                    f'RNA-seq file {rna_file_obj.accession} has unexpected content type for tar format.')
            self._check_file_name(
                rna_file_obj, expected_file_substring=f'{self.analysis_set_acc}.tar.gz')


class ATACFileChecker(BaseFileChecker):
    """Checker for ATAC-seq files."""

    def check_atac_file_count(self, active_file_objs: list[const.FileObjTypes], file_type: str):
        """Check if exactly 1 ATAC file is present."""
        if len(active_file_objs) != 1:
            raise ValueError(
                f'Analysis set {self.analysis_set_acc} has more than 1 {file_type} file.')

    def check_alignment_file_specifics(self, alignment_file_obj: const.FileObjTypes):
        """Check alignment file specific properties."""
        if alignment_file_obj.file_format != 'bam':
            raise ValueError(
                f'ATAC alignment file {alignment_file_obj.accession} has unexpected file format.')
        if alignment_file_obj.content_type != 'alignments':
            raise ValueError(
                f'ATAC alignment file {alignment_file_obj.accession} has unexpected content type.')
        self._check_file_name(
            alignment_file_obj, expected_file_substring=f'{self.analysis_set_acc}.bam')

    def check_fragment_file_specifics(self, fragment_file_obj: const.FileObjTypes):
        """Check fragment file specific properties."""
        if fragment_file_obj.file_format != 'bed':
            raise ValueError(
                f'ATAC fragment file {fragment_file_obj.accession} has unexpected file format.')
        if fragment_file_obj.file_format_type != 'bed3+':
            raise ValueError(
                f'ATAC fragment file {fragment_file_obj.accession} has unexpected file format type.')
        if fragment_file_obj.content_type != 'fragments':
            raise ValueError(
                f'ATAC fragment file {fragment_file_obj.accession} has unexpected content type.')
        self._check_file_name(
            fragment_file_obj, expected_file_substring=f'{self.analysis_set_acc}.fragments.tsv.gz')

    def check_fragment_controlled_access(self, fragment_file_obj: const.FileObjTypes):
        """Check fragment file specific controlled access requirements."""
        if fragment_file_obj.controlled_access is True:
            raise ValueError(
                f'URGENT: Fragment file {fragment_file_obj.accession} should not be controlled_access files.')


class QualityCheckAnalysisSet:
    """Main class to run quality checks on analysis sets."""

    def __init__(self, analysis_set_acc: str, igvf_client_api: const.IgvfApiType):
        self.analysis_set_acc = analysis_set_acc
        self.igvf_client_api = igvf_client_api
        self.analysis_set_obj = self.igvf_client_api.get_by_id(
            f'/analysis-sets/{self.analysis_set_acc}').actual_instance
        self.analysis_set_files = self.analysis_set_obj.files
        self.controlled_access = _check_if_controlled_access(
            self.analysis_set_obj, self.igvf_client_api)
        self.pipeline_parameters = self.analysis_set_obj.pipeline_parameters

        # Initialize checkers
        self.doc_checker = DocumentChecker(analysis_set_acc, igvf_client_api)
        self.qc_checker = QCMetricsChecker(analysis_set_acc, igvf_client_api)
        self.index_checker = IndexFileChecker(
            analysis_set_acc, igvf_client_api)
        self.rna_checker = RNASeqFileChecker(analysis_set_acc, igvf_client_api)
        self.atac_checker = ATACFileChecker(analysis_set_acc, igvf_client_api)

    def _check_has_files(self):
        """Check if the analysis set has files."""
        if not self.analysis_set_obj.files:
            raise ValueError(
                f'URGENT: Analysis set {self.analysis_set_acc} has no files.')

    def _check_rnaseq_data(self):
        """Check RNA-seq data for the analysis set."""
        all_errors = []

        rna_file_ids = [output_file for output_file in self.analysis_set_files
                        if output_file.startswith('/matrix-files/')]
        active_rna_file_objs = _get_active_file_objs(
            file_ids=rna_file_ids, igvf_client_api=self.igvf_client_api)

        try:
            self.rna_checker.check_rna_file_count(active_rna_file_objs)
        except ValueError as e:
            all_errors.append(str(e))

        for rna_file_obj in active_rna_file_objs:
            # Collect basic file property errors
            basic_errors = self.rna_checker._check_basic_file_properties(
                rna_file_obj, const.RNA_ANALYSIS_STEP_VERSION)
            all_errors.extend(basic_errors)

            try:
                self.rna_checker._check_reference_files(
                    rna_file_obj, const.RNASEQ_REF_FILES)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.qc_checker.check_rna_qc_metrics(rna_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.rna_checker.check_rna_file_format_and_content(
                    rna_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

        # Raise all errors if any found
        if all_errors:
            raise ValueError(";".join(all_errors))

    def _check_atac_alignment_data(self):
        """Check ATAC alignment data for the analysis set."""
        all_errors = []

        alignment_file_ids = [file_id for file_id in self.analysis_set_files
                              if file_id.startswith('/alignment-files/')]
        active_alignment_file_objs = _get_active_file_objs(
            file_ids=alignment_file_ids, igvf_client_api=self.igvf_client_api)

        try:
            self.atac_checker.check_atac_file_count(
                active_alignment_file_objs, "alignment")
        except ValueError as e:
            all_errors.append(str(e))

        if active_alignment_file_objs:
            alignment_file_obj = active_alignment_file_objs[0]

            # Collect basic file property errors
            basic_errors = self.atac_checker._check_basic_file_properties(
                alignment_file_obj, const.ATAC_ANALYSIS_STEP_VERSION, self.controlled_access)
            all_errors.extend(basic_errors)

            try:
                self.atac_checker.check_alignment_file_specifics(
                    alignment_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.atac_checker._check_reference_files(
                    alignment_file_obj, const.ATAC_REF_FILES)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.qc_checker.check_atac_qc_metrics(alignment_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.index_checker.check_index_file_presence_and_format(
                    alignment_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

        # Raise all errors if any found
        if all_errors:
            raise ValueError(";".join(all_errors))

    def _check_atac_fragment_data(self):
        """Check ATAC fragment data for the analysis set."""
        all_errors = []

        fragment_file_ids = [file_id for file_id in self.analysis_set_files
                             if file_id.startswith('/tabular-files/')]
        active_fragment_file_objs = _get_active_file_objs(
            file_ids=fragment_file_ids, igvf_client_api=self.igvf_client_api)

        try:
            self.atac_checker.check_atac_file_count(
                active_fragment_file_objs, "fragment")
        except ValueError as e:
            all_errors.append(str(e))

        if active_fragment_file_objs:
            fragment_file_obj = active_fragment_file_objs[0]

            # Collect basic file property errors
            basic_errors = self.atac_checker._check_basic_file_properties(
                fragment_file_obj, const.ATAC_ANALYSIS_STEP_VERSION)
            all_errors.extend(basic_errors)

            try:
                self.atac_checker.check_fragment_file_specifics(
                    fragment_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.atac_checker.check_fragment_controlled_access(
                    fragment_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.atac_checker._check_reference_files(
                    fragment_file_obj, const.ATAC_REF_FILES)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.qc_checker.check_atac_qc_metrics(fragment_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

            try:
                self.index_checker.check_index_file_presence_and_format(
                    fragment_file_obj)
            except ValueError as e:
                all_errors.append(str(e))

        # Raise all errors if any found
        if all_errors:
            raise ValueError(";".join(all_errors))

    def run_all_checks(self):
        """Run all checks for the analysis set and collect all errors."""
        errors = []

        # Check if it has files
        try:
            self._check_has_files()
        except ValueError as e:
            errors.append(str(e))

        # Check pipeline parameters document
        try:
            self.doc_checker.check_pipeline_parameters_document(
                self.pipeline_parameters)
        except ValueError as e:
            errors.append(str(e))

        # Check if RNA-seq data is present
        rna_file_ids = [output_file for output_file in self.analysis_set_files
                        if output_file.startswith('/matrix-files/')]
        if rna_file_ids:
            try:
                self._check_rnaseq_data()
            except ValueError as e:
                errors.append(str(e))

        # Check if ATAC data is present
        atac_file_ids = [output_file for output_file in self.analysis_set_files
                         if output_file.startswith(('/alignment-files/', '/tabular-files/'))]
        if atac_file_ids:
            try:
                self._check_atac_alignment_data()
            except ValueError as e:
                errors.append(str(e))

            try:
                self._check_atac_fragment_data()
            except ValueError as e:
                errors.append(str(e))

        # Return all errors as a single string
        if errors:
            # Parse the errors first by splitting on spaces
            parsed_errors = []
            for error in errors:
                parsed_errors.extend(error.split(';'))

            numbered_errors = {f"Error_{i+1}": error for i,
                               error in enumerate(parsed_errors)}
            return numbered_errors
        else:
            return "All checks passed."


def main(list_of_analysis_set_acc: list[str], igvf_client_api: const.IgvfApiType, output_file_path: str) -> dict[str, str]:
    """Main function to run quality checks on a list of analysis set accessions."""
    results = {}
    for analysis_set_acc in list_of_analysis_set_acc:
        qc_checker = QualityCheckAnalysisSet(
            analysis_set_acc=analysis_set_acc, igvf_client_api=igvf_client_api)
        result = qc_checker.run_all_checks()
        results[analysis_set_acc] = result

    # Write results to JSON file
    with open(output_file_path, 'w') as f:
        json.dump(results, f, indent=4)
    logging.info(f'QC results written to {output_file_path}.')

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument('--igvf_endpoint', '-m', type=str, choices=['sandbox', 'prod', 'staging'], default='prod',
                        help="""The IGVF endpoint, sandbox, prod, or staging.""")

    group.add_argument('--qc_analysis_set_file', '-i', type=str,
                       help="""A file containing a list of analysis set accessions to be checked.""")
    group.add_argument('--qc_analysis_sets', '-r', type=str,
                       help="""A comma-separated list of analysis set accessions to be checked.""")

    parser.add_argument('--output_file_path', '-o', type=str, default='./qc_results.json',
                        help="""Output file with quality check results included.""")

    args = parser.parse_args()

    # Get the IGVF client API, production or sandbox
    igvf_api_keys = api_tools.set_up_api_keys(
        igvf_endpoint=args.igvf_endpoint)

    # Set up IGVF client API
    igvf_portal_api = api_tools.get_igvf_client_auth(igvf_api_keys=igvf_api_keys,
                                                     igvf_endpoint=args.igvf_endpoint)

    # Get a list of analysis set accessions from the arg or file, but not both
    if args.qc_analysis_set_file:
        with open(args.qc_analysis_set_file, 'r') as f:
            analysis_set_accs_to_check = ','.join(f.read().splitlines())
    elif args.qc_analysis_sets:
        analysis_set_accs_to_check = args.qc_analysis_sets

    # Run main function
    main(list_of_analysis_set_acc=analysis_set_accs_to_check.split(','),
         igvf_client_api=igvf_portal_api,
         output_file_path=args.output_file_path)
