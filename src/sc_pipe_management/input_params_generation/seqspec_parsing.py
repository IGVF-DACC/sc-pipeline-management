import dataclasses
import subprocess
from argparse import Namespace
from pathlib import Path
import io
from contextlib import redirect_stdout

# Import seqspec modules
import seqspec
from seqspec.utils import load_spec
from seqspec.seqspec_index import run_index
from seqspec.seqspec_onlist import run_onlist
from seqspec.seqspec_file import run_file

# Import IGVF API client
import igvf_client

# Import custom modules
import sc_pipe_management.input_params_generation.constant as const
import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing


def _capture_stdout_from_command(command_func, args) -> str:
    """Helper function to capture stdout from a seqspec command. The output will contain newlines."""
    buf = io.StringIO()
    with redirect_stdout(buf):
        command_func(None, args)
    return buf.getvalue().strip()


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

    def get_onlist_read_id(self) -> str:
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
        if const.IGVF_ACCESSION_REGEX.match(read_id):
            return const.IGVF_ACCESSION_REGEX.search(read_id).group(1)
        # If read_id is not an IGVF accession, check the file_ids
        file_specs = single_read_spec.files
        for file_spec in file_specs:
            file_id = file_spec.file_id
            # If the read_id is not an IGVF accession, check the file_id
            if const.IGVF_ACCESSION_REGEX.match(file_id):
                return const.IGVF_ACCESSION_REGEX.search(file_id).group(1)
            # If file_id is not an IGVF accession, check the URL
            for item in file_spec.url.split('/'):
                if const.IGVF_ACCESSION_REGEX.match(item):
                    return const.IGVF_ACCESSION_REGEX.search(item).group(1)
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

    def _get_onlist_files(self) -> list[str]:
        """Get the onlist files from the seqspec file."""
        # Set up seqspec file args
        file_args = Namespace(
            yaml=Path(self.seqspec_file_path),
            modality=self._get_seqspec_modality(),
            selector='region-type',
            ids='barcode',
            format='index',
            key='url',
            output=None,
            fullpath=None
        )
        # Run it while capturing stdout
        stdout = _capture_stdout_from_command(run_file, file_args)
        # Error handling if no onlist file found
        if not stdout:
            raise const.BadDataException(
                f"No onlist file found in the seqspec YAML file {Path(self.seqspec_file_path).name}.")
        else:
            return sorted(stdout.strip().split(','))

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
        self.seqspec_onlist_read_id = self.seqspec_metadata.get_onlist_read_id()

    def _get_seqspec_index(self) -> str:
        """Get the index read IDs based on the usage purpose."""
        index_args = Namespace(
            yaml=Path(self.seqspec_file_path),
            output=None,
            subregion_type=None,
            tool=const.ASSAY_TYPE_TO_TOOL_FORMAT[self.modality],
            selector='read',
            rev=False,
            region=False,
            modality=self.modality,
            ids=self.seqspec_index_read_ids,
            r=None
        )
        index_res = _capture_stdout_from_command(run_index, index_args)
        if const.ASSAY_TYPE_TO_TOOL_FORMAT[self.modality] == 'chromap':
            return index_res.split('--read-format')[1].strip()
        return index_res

    def _get_seqspec_onlist(self) -> str:
        """Get the onlist files from the seqspec file."""
        if self.onlist_method == 'no combination':
            onlist_args = Namespace(
                yaml=Path(self.seqspec_file_path),
                modality=self.modality,
                selector="region-type",
                id="barcode",
                output=self.output_barcode_list_file,
                format=None,
                r=None  # r is a deprecated arg, required but default to None
            )
        else:
            # If combinatorial, seqspec onlist will generated a new file
            # NOTE: Per FG mtg, on Mar 10, 2025, all combinatorial onlist files will be generated as product.
            onlist_args = Namespace(
                yaml=Path(self.seqspec_file_path),
                modality=self.modality,
                selector='read',
                id=self.seqspec_onlist_read_id,
                output=self.output_barcode_list_file,
                format='product',
                r=None
            )
        # Run the onlist command with the args
        return _capture_stdout_from_command(run_onlist, onlist_args)

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
        # Custom defined error and error types in seqspec.utils
        except (const.BadDataException, ValueError, IndexError) as e:
            return SeqSpecToolOutput(
                read_index_string=read_index_string,
                final_barcode_file=onlist_file_path,
                errors=str(e)
            )
