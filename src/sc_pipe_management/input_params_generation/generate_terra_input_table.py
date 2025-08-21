def get_atac_seqfile_read_titles(read_names: list) -> str:
    """Get ATACseq read titles. Will always have read1, read2, and barcode index even if read2 and barcode may be concatenated.

    Args:
        read_names (list): seqfile.read_names values

    Returns:
        str: assay-type_read-names, e.g., rna_read1
    """
    assay_type = 'atac'
    read_file_titles = []
    for read_name in read_names:
        read_file_titles.append(
            '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    return read_file_titles


def get_rna_seqfile_read_titles(read_names: list) -> str:
    """Get RNA seq file read titles. If reads are not concatenated, will always have read1, read2, and barcode index.
    If reads are concatenated, will only have read1 and read2. The barcode fastaq will be ignored.

    Args:
        read_names (list): seqfile.read_names values

    Returns:
        str: assay-type_read-names, e.g., rna_read1
    """
    assay_type = 'rna'
    read_file_titles = []
    if len(read_names) == 1:
        for read_name in read_names:
            read_file_titles.append(
                '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    else:
        # NOTE: Hard coded for cases if read2 and barcode are concatenated. Ignore barcode index.
        if sorted(read_names) == sorted(['Read 2', 'Barcode index']):
            read_file_titles.append(
                '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP['Read 2']]))
        else:
            read_file_titles.append(
                '_'.join([assay_type, const.READ_NAME_TO_READ_TYPE_MAP[read_name]]))
    return read_file_titles


def get_seqfile_readnames(seqfile_item, assay_type: str) -> str:
    """HACK: Use sequence file metadata to infer FASTQ file read types.

    Args:
        seqfile_item (_type_): Return item from the API
        assay_type (str): translated from assay term into atac or rna

    Returns:
        str: e.g., atac_read1, rna_read2
    """
    read_file_titles = []
    curr_read_names = seqfile_item.read_names
    if not curr_read_names:
        raise const.BadDataException(
            'Error: No read names found in the sequence file item.')
    if assay_type == 'atac':
        return get_atac_seqfile_read_titles(read_names=curr_read_names)
    elif assay_type == 'rna':
        return get_rna_seqfile_read_titles(read_names=curr_read_names)


def download_file_via_https(seqspec_file_url: str, partial_root_dir: str) -> str:
    """Download file via href url on the portal and save it locally."""
    if not os.path.exists(partial_root_dir):
        os.makedirs(partial_root_dir)
    # TODO: Needs to update this to use API tools
    username = os.getenv('IGVF_API_KEY')
    password = os.getenv('IGVF_SECRET_KEY')
    session = requests.Session()
    session.auth = (username, password)
    response = session.get(seqspec_file_url)
    if response.status_code == 200:
        curr_output_file = os.path.join(
            partial_root_dir, 'temp', seqspec_file_url.split('/')[-1])
        with open(curr_output_file, 'wb') as file:
            file.write(response.content)
        return curr_output_file
    else:
        raise const.BadDataException(
            f'Error: Download failed with status code: {response.status_code}.')

    def _check_emptyinclusion_list(self) -> bool:
        """Check if the final inclusion list created is empty."""
        with open(self.output_barcode_list_file, 'r') as file_obj:
            first_char = file_obj.read(1)
            if not first_char.isalpha():
                return False


@dataclasses.dataclass(frozen=True)
class SeqSpecQualityCheckResult:
    """Dataclass to hold seqspec quality check result."""
    seqspec_onlist_errors: str | None
    seqspec_index_errors: str | None
    seqspec_discrepancy_errors: str | None
    seqspec_measet_discrepancy_errors: str | None


class SeqSpecQualityCheck:
    """Class to perform quality check on seqspec metadata."""

    def __init__(self, anaset_metadata: portal_parsing.InputAnalysisSetMetadata, igvf_api: igvf_client.api.igvf_api.IgvfApi):
        self.anaset_metadata = anaset_metadata
        self.igvf_api = igvf_api
        self.rna_measet_metadata = self.anaset_metadata.rna_input_info
        self.atac_measet_metadata = self.anaset_metadata.atac_input_info

    def _check_onlist_files(self) -> str | None:
        """Check if the onlist files are valid."""
        if not self.seqspec_metadata.onlist_files:
            return 'No onlist files found in the seqspec file.'
        return None

    def _check_index_read_ids(self) -> str | None:
        """Check if the index read IDs are valid."""
        if not self.seqspec_metadata.get_index_read_ids():
            return 'No index read IDs found in the seqspec file.'
        return None

    def _check_discrepancies(self) -> str | None:
        """Check for discrepancies in the seqspec metadata."""
        # Placeholder for actual discrepancy checks
        return None

    def _check_measet_discrepancies(self) -> str | None:
        """Check for discrepancies in the measet."""
        # Placeholder for actual measet discrepancy checks
        return None

    def run_quality_check(self) -> SeqSpecQualityCheckResult:
        """Run quality check on the seqspec metadata."""
        onlist_errors = self._check_onlist_files()
        index_errors = self._check_index_read_ids()
        discrepancy_errors = self._check_discrepancies()
        measet_discrepancy_errors = self._check_measet_discrepancies()

        return SeqSpecQualityCheckResult(
            seqspec_onlist_errors=onlist_errors,
            seqspec_index_errors=index_errors,
            seqspec_discrepancy_errors=discrepancy_errors,
            seqspec_measet_discrepancy_errors=measet_discrepancy_errors
        )
