import dataclasses
import subprocess

# Import seqspec modules
from seqspec.utils import load_spec
from seqspec.seqspec_index import run_index
import seqspec

# Import IGVF API client
import igvf_client

# Import custom modules
import sc_pipe_management.input_params_generation.constant as const
import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing


@dataclasses.dataclass(frozen=True)
class SeqSpecMetadata:
    """Dataclass to hold seqspec metadata."""
    seqspec_file_path: str
    modality: str
    ordered_read_ids: list
    onlist_files: list[str]

    def get_index_read_ids(self) -> str:
        """Get the index read IDs based on the usage purpose."""
        return ','.join(self.ordered_read_ids)

    def get_onlist_read_ids(self) -> str:
        """Get the read IDs for onlist usage purpose."""
        return self.ordered_read_ids[-1]


class GetSeqSpecMetadata:
    """Class to get seqspec metadata from a seqspec file."""

    def __init__(self, seqspec_file_path: str, igvf_api: igvf_client.api.igvf_api.IgvfApi):
        self.seqspec_file_path = seqspec_file_path
        self.igvf_api = igvf_api
        self.imported_seqspec = load_spec(self.seqspec_file_path)

    def _get_seqspec_modality(self) -> str:
        modalities = self.imported_seqspec.modalities
        if len(modalities) != 1:
            raise const.BadDataException(
                f'Error: Seqspec file {self.seqspec_file_path.split("/")[-2]} has multiple modalities. Only one modality is supported.')
        return modalities[0]

    def _find_igvf_acc_in_single_read_spec(self, single_read_spec: seqspec.Read.Read) -> str | None:
        """Find IGVF accession in a single !Read spec. Trying read_id first, then file_ids, and finally the URL."""
        # If read_id is an IGVF accession, return it
        read_id = single_read_spec.read_id
        if const.READ_ID_REGEX.match(read_id):
            return const.READ_ID_REGEX.search(read_id).group(1)
        # If read_id is not an IGVF accession, check the file_ids
        file_specs = single_read_spec.files
        for file_spec in file_specs:
            file_id = file_spec.file_id
            # If the read_id is not an IGVF accession, check the file_id
            if const.READ_ID_REGEX.match(file_id):
                return const.READ_ID_REGEX.search(file_id).group(1)
            # If file_id is not an IGVF accession, check the URL
            for item in file_spec.url.split('/'):
                if const.READ_ID_REGEX.match(item):
                    return const.READ_ID_REGEX.search(item).group(1)
        return None

    def _parse_read_ids_from_seqspec_file(self) -> list[tuple[str, str | None]]:
        """Retrieve read_id values of a specific modality from a seqspec."""
        sequence_specs = self.imported_seqspec.sequence_spec
        read_ids = []
        for spec in sequence_specs:
            if spec.modality == self._get_seqspec_modality():
                igvf_accession = self._find_igvf_acc_in_single_read_spec(
                    single_read_spec=spec)
                read_ids.append((spec.read_id, igvf_accession))
        return read_ids

    def _generate_ordered_read_ids(self, seqfiles_metadata: list[portal_parsing.SeqFileMetadata], assay_type: str) -> list[str]:
        """Get the read1,read2,barcode files order from the seqspec read_id and portal metadata."""
        # Get read IDs from the seqspec file
        read_ids = self._parse_read_ids_from_seqspec_file()
        id_by_name = {}
        for (read_id, igvf_accession) in read_ids:
            # If no IGVF accession is found, add to unknown
            if igvf_accession is None:
                raise const.BadDataException(
                    f'Error: No IGVF accession found for read ID {read_id}.')
            for seqfile_metadata in seqfiles_metadata:
                if seqfile_metadata.file_accession == igvf_accession:
                    read_names = seqfile_metadata.get_read_names_from_seqfile(
                        assay_type=assay_type)
                    for read_name in read_names:
                        if read_name in id_by_name:
                            raise const.BadDataException(
                                f'Error: Multiple sequence files found for the same {read_name}.')
                        id_by_name[read_name] = read_id
        # Ensure the order is read1, read2, barcode, and remove duplicates if any
        result = []
        seen = set()
        for read_name in ['read1', 'read2', 'barcode']:
            val = id_by_name.get(read_name)
            if val is not None and val not in seen:
                result.append(val)
                seen.add(val)
        return result

    def _parse_onlist_files(self, onlist_files: list[str]) -> list[str]:
        # If it's a list with one string containing commas, split it
        if len(onlist_files) == 1 and ',' in onlist_files[0]:
            return sorted(onlist_files[0].split(','))
        return sorted(onlist_files)

    def _get_onlist_files(self) -> list[str]:
        """Get the onlist files from the seqspec file."""
        curr_run_log = subprocess.run(['seqspec', 'file', '-m', self._get_seqspec_modality(), '-s', 'region-type',
                                      '-i', 'barcode', '-f', 'index', '-k', 'url', self.seqspec_file_path], capture_output=True)
        if curr_run_log.returncode != 0:
            raise const.BadDataException(
                f'Error: seqspec file {self.seqspec_file_path} does not have onlist files.')
        onlist_files = curr_run_log.stdout.decode('utf-8').strip().split('\n')
        return self._parse_onlist_files(onlist_files)

    def generate_seqspec_metadata(self, seqfiles_metadata: list[portal_parsing.SeqFileMetadata], assay_type: str) -> SeqSpecMetadata:
        """Generate seqspec metadata from the seqspec file and seqfiles metadata."""
        # Generate ordered read IDs
        ordered_read_ids = self._generate_ordered_read_ids(
            seqfiles_metadata=seqfiles_metadata, assay_type=assay_type)
        # Return the seqspec metadata
        return SeqSpecMetadata(
            seqspec_file_path=self.seqspec_file_path,
            modality=self._get_seqspec_modality(),
            ordered_read_ids=ordered_read_ids,
            onlist_files=self._get_onlist_files()
        )


@dataclasses.dataclass(frozen=True)
class SeqSpecToolOutput:
    """Dataclass to hold seqspec tool output."""
    read_index_string: str
    final_barcode_file: str
    errors: str | None


class GetSeqSpecToolOutput:
    """Class to generate seqspec metadata for a given seqspec file."""

    def __init__(self, seqspec_metadata: SeqSpecMetadata, onlist_method: str, output_barcode_list_file: str):
        self.seqspec_metadata = seqspec_metadata
        self.onlist_method = onlist_method
        self.output_barcode_list_file = output_barcode_list_file
        self.modality = self.seqspec_metadata.modality
        self.seqspec_file_path = self.seqspec_metadata.seqspec_file_path
        self.seqspec_index_read_ids = self.seqspec_metadata.get_index_read_ids()
        self.seqspec_onlist_read_id = self.seqspec_metadata.get_onlist_read_ids()

    def _get_seqspec_index(self) -> str:
        """Get the index read IDs based on the usage purpose."""
        curr_run_log = subprocess.run(['seqspec', 'index', '-m', self.modality, '-t',
                                       const.ASSAY_TYPE_TO_TOOL_FORMAT[self.modality], '-s', 'read', '-i', self.seqspec_index_read_ids, self.seqspec_file_path], capture_output=True)
        if curr_run_log.returncode != 0:
            raise const.BadDataException(
                f'Error: seqspec file {self.seqspec_file_path} does not have index read IDs.')
        if self.modality == 'rna':
            return curr_run_log.stdout.decode('utf-8').strip()
        elif self.modality == 'atac':
            return curr_run_log.stdout.decode('utf-8').strip().split(' ')[-1]

    def _get_seqspec_onlist(self) -> str:
        """Get the onlist files from the seqspec file."""
        if self.onlist_method == 'no combination':
            curr_run_log = subprocess.run(['seqspec', 'onlist', '-m', self.modality, '-s',
                                           'region-type', '-i', 'barcode', '-o', self.output_barcode_list_file, self.seqspec_file_path], capture_output=True)
        # If combinatorial, seqspec onlist will generated a new file
        else:
            # NOTE: Per sc FG mtg, on Mar 10, 2025, all combinatorial onlist files will be generated as product.
            curr_run_log = subprocess.run(['seqspec', 'onlist', '-m', self.modality, '-s', 'read', '-i',
                                           self.seqspec_onlist_read_id, '-f', 'product', '-o', self.output_barcode_list_file, self.seqspec_file_path], capture_output=True)
        if curr_run_log.returncode != 0:
            raise const.BadDataException(
                f'Error: Seqspec onlist generation command error. Seqspec tool debug msg: {curr_run_log.stderr.decode("utf-8")}.')
        return self.output_barcode_list_file

    def generate_seqspec_tool_output(self) -> SeqSpecToolOutput:
        """Generate seqspec tool output for a given seqspec file."""
        read_index_string = None
        onlist_file_path = None
        try:
            # Get the index read IDs
            read_index_string = self._get_seqspec_index()
            # Get the onlist files
            onlist_file_path = self._get_seqspec_onlist()
            return SeqSpecToolOutput(
                read_index_string=read_index_string,
                final_barcode_file=onlist_file_path,
                errors=None
            )
        except const.BadDataException as e:
            return SeqSpecToolOutput(
                read_index_string=read_index_string,
                final_barcode_file=onlist_file_path,
                errors=str(e)
            )
