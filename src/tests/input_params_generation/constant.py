import sys
import os
import dataclasses

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.input_params_generation.portal_metadata_parsing as portal_parsing
import sc_pipe_management.input_params_generation.seqspec_parsing as seqspec_parsing


@dataclasses.dataclass(frozen=True)
class SeqSpecTestConstants:
    seqspec_file_path: str
    modality: str
    ordered_read_ids: list[str]
    onlist_files: list[str]
    read_index: str
    index_input: str
    onlist_input: str
    onlist_method: str
    final_barcode_file: str
    final_barcode_file_md5: str

    def convert_to_seqspec_metadata(self) -> seqspec_parsing.SeqSpecMetadata:
        return seqspec_parsing.SeqSpecMetadata(
            seqspec_file_path=self.seqspec_file_path,
            modality=self.modality,
            ordered_read_ids=self.ordered_read_ids,
            onlist_files=self.onlist_files
        )


SEQSPEC_PARSE_TEST_RESULTS = {
    # Analysis Set: IGVFDS9664YLAD
    '10x multiome': {
        'atac': SeqSpecTestConstants(
            seqspec_file_path='src/tests/test_files/IGVFFI1365BOUI.yaml.gz',
            modality='atac',
            ordered_read_ids=["IGVFFI4665EVGC",
                              "IGVFFI4986VMDU", "IGVFFI0391NHGA"],
            onlist_files=[
                'https://api.data.igvf.org/tabular-files/IGVFFI7587TJLC/@@download/IGVFFI7587TJLC.tsv.gz'
            ],
            onlist_method='no combination',
            read_index='bc:8:23:-,r1:0:49,r2:0:49',
            index_input='IGVFFI4665EVGC,IGVFFI4986VMDU,IGVFFI0391NHGA',
            onlist_input='IGVFFI0391NHGA',
            final_barcode_file='src/tests/test_files/10x multiome_atac_final_barcode_list.txt',
            final_barcode_file_md5='59ed4e3b7c64af921bde8b8497924c02'
        ),
        'rna': SeqSpecTestConstants(
            seqspec_file_path='src/tests/test_files/IGVFFI2433YUEJ.yaml.gz',
            modality='rna',
            ordered_read_ids=["IGVFFI6948UZJO", "IGVFFI7951DAQB"],
            onlist_files=[
                'https://api.data.igvf.org/tabular-files/IGVFFI8751YQRY/@@download/IGVFFI8751YQRY.tsv.gz'
            ],
            onlist_method='no combination',
            read_index='0,0,16:0,16,28:1,0,90',
            index_input='IGVFFI6948UZJO,IGVFFI7951DAQB',
            onlist_input='IGVFFI7951DAQB',
            final_barcode_file='src/tests/test_files/10x multiome_rna_final_barcode_list.txt',
            final_barcode_file_md5='5b223bbf2d45ba5e8bb55787872321c8'
        )
    },
    # Analysis set: IGVFDS0223KLTB
    'parse splitseq': {
        'rna': SeqSpecTestConstants(
            seqspec_file_path='src/tests/test_files/IGVFFI2264BQQD.yaml.gz',
            modality='rna',
            ordered_read_ids=["IGVFFI4633WXZI.fastq.gz",
                              "IGVFFI9320MADV.fastq.gz"],
            onlist_files=[
                "https://api.data.igvf.org/tabular-files/IGVFFI0924TKJO/@@download/IGVFFI0924TKJO.tsv.gz",
                "https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz",
                "https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz"
            ],
            onlist_method='multi',
            read_index='1,10,18,1,48,56,1,78,86:1,0,10:0,0,140',
            index_input='IGVFFI4633WXZI.fastq.gz,IGVFFI9320MADV.fastq.gz',
            onlist_input='IGVFFI9320MADV.fastq.gz',
            final_barcode_file='src/tests/test_files/parse splitseq_rna_final_barcode_list.txt',
            final_barcode_file_md5='7f1e51743561883977e5a905bb345ac8'
        )
    },
    # Analysis set: IGVFDS0657NHTA
    'shareseq': {
        'atac': SeqSpecTestConstants(
            seqspec_file_path='src/tests/test_files/IGVFFI8012OCZQ.yaml.gz',
            modality='atac',
            ordered_read_ids=["IGVFFI1918YZDJ.fastq.gz",
                              "IGVFFI4773YQEF.fastq.gz"],
            onlist_files=[
                "https://api.data.igvf.org/tabular-files/IGVFFI6680XCDT/@@download/IGVFFI6680XCDT.tsv.gz",
                "https://api.data.igvf.org/tabular-files/IGVFFI2668ZFEB/@@download/IGVFFI2668ZFEB.tsv.gz",
                "https://api.data.igvf.org/tabular-files/IGVFFI2668ZFEB/@@download/IGVFFI2668ZFEB.tsv.gz"
            ],
            onlist_method='product',
            read_index='bc:115:122,bc:153:160,bc:191:198,r1:0:99,r2:0:99',
            index_input='IGVFFI1918YZDJ.fastq.gz,IGVFFI4773YQEF.fastq.gz',
            onlist_input='IGVFFI4773YQEF.fastq.gz',
            final_barcode_file='src/tests/test_files/shareseq_atac_final_barcode_list.txt',
            final_barcode_file_md5='f7eabd710222e594459698a6c624de7b'
        ),
        'rna': SeqSpecTestConstants(
            seqspec_file_path='src/tests/test_files/IGVFFI5825ATCM.yaml.gz',
            modality='rna',
            ordered_read_ids=["IGVFFI0777CMJH.fastq.gz",
                              "IGVFFI7330UNCB.fastq.gz"],
            onlist_files=[
                "https://api.data.igvf.org/tabular-files/IGVFFI6680XCDT/@@download/IGVFFI6680XCDT.tsv.gz",
                "https://api.data.igvf.org/tabular-files/IGVFFI2668ZFEB/@@download/IGVFFI2668ZFEB.tsv.gz",
                "https://api.data.igvf.org/tabular-files/IGVFFI2668ZFEB/@@download/IGVFFI2668ZFEB.tsv.gz"
            ],
            onlist_method='product',
            read_index='1,115,123,1,153,161,1,191,199:1,0,10:0,0,100',
            index_input='IGVFFI0777CMJH.fastq.gz,IGVFFI7330UNCB.fastq.gz',
            onlist_input='IGVFFI7330UNCB.fastq.gz',
            final_barcode_file='src/tests/test_files/shareseq_rna_final_barcode_list.txt',
            final_barcode_file_md5='f7eabd710222e594459698a6c624de7b'
        )
    },
    # Analysis set: IGVFDS7784IUJF
    '10x multiome double': {
        'rna': SeqSpecTestConstants(
            seqspec_file_path='src/tests/test_files/IGVFFI8568RJIP.yaml.gz',
            modality='rna',
            ordered_read_ids=["IGVFFI4117HDFW", "IGVFFI9675ROSO"],
            onlist_files=[
                "https://api.data.igvf.org/tabular-files/IGVFFI8751YQRY/@@download/IGVFFI8751YQRY.tsv.gz"
            ],
            onlist_method='no combination',
            read_index='0,0,16:0,16,28:1,0,101',
            index_input='IGVFFI4117HDFW,IGVFFI9675ROSO',
            onlist_input='IGVFFI9675ROSO',
            final_barcode_file='src/tests/test_files/10x multiome double_rna_final_barcode_list.txt',
            final_barcode_file_md5='5b223bbf2d45ba5e8bb55787872321c8'
        )
    },
    # Analysis set: IGVFDS4790PPFI
    'rnaseq with bad read_id': {
        'rna': SeqSpecTestConstants(
            seqspec_file_path='src/tests/test_files/IGVFFI1382HZDV.yaml.gz',
            modality='rna',
            ordered_read_ids=["RNA Read 1", "RNA Read 2"],
            onlist_files=[
                'https://api.data.igvf.org/tabular-files/IGVFFI8751YQRY/@@download/IGVFFI8751YQRY.tsv.gz'
            ],
            onlist_method='no combination',
            read_index='0,0,16:0,16,28:1,0,90',
            index_input='RNA Read 1,RNA Read 2',
            onlist_input='RNA Read 2',
            final_barcode_file='src/tests/test_files/rnaseq with bad read_id_rna_final_barcode_list.txt',
            final_barcode_file_md5='5b223bbf2d45ba5e8bb55787872321c8'
        )
    }
}


@dataclasses.dataclass(frozen=True)
class TestSeqFileMetadata:
    file_accession: str
    read_names: list[str]


TEST_SEQFILES_METADATA = [portal_parsing.SeqFileMetadata(file_accession='IGVFFI6948UZJO', file_set='/measurement-sets/IGVFDS9139TBCB/', illumina_read_type='R1', sequencing_run=1, lane=4, flowcell_id='223JHYLT3', file_url='https://api.data.igvf.org/sequence-files/IGVFFI6948UZJO/@@download/IGVFFI6948UZJO.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI2433YUEJ/@@download/IGVFFI2433YUEJ.yaml.gz'], read_names=['Read 1']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI7951DAQB', file_set='/measurement-sets/IGVFDS9139TBCB/', illumina_read_type='R2', sequencing_run=1, lane=4, flowcell_id='223JHYLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI7951DAQB/@@download/IGVFFI7951DAQB.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI2433YUEJ/@@download/IGVFFI2433YUEJ.yaml.gz'], read_names=['Read 2']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI0391NHGA', file_set='/measurement-sets/IGVFDS1013HCXI/', illumina_read_type='R2', sequencing_run=1, lane=8, flowcell_id='223JKJLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI0391NHGA/@@download/IGVFFI0391NHGA.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI1365BOUI/@@download/IGVFFI1365BOUI.yaml.gz'], read_names=['Barcode index']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI4665EVGC', file_set='/measurement-sets/IGVFDS1013HCXI/', illumina_read_type='R1', sequencing_run=1, lane=8, flowcell_id='223JKJLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI4665EVGC/@@download/IGVFFI4665EVGC.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI1365BOUI/@@download/IGVFFI1365BOUI.yaml.gz'], read_names=['Read 1']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI4986VMDU', file_set='/measurement-sets/IGVFDS1013HCXI/', illumina_read_type='R3', sequencing_run=1, lane=8, flowcell_id='223JKJLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI4986VMDU/@@download/IGVFFI4986VMDU.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI1365BOUI/@@download/IGVFFI1365BOUI.yaml.gz'], read_names=['Read 2']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI4633WXZI', file_set='/measurement-sets/IGVFDS5272EQRG/', illumina_read_type='R1', sequencing_run=2, lane=1, flowcell_id='AAC2JNKHV',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI4633WXZI/@@download/IGVFFI4633WXZI.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI2264BQQD/@@download/IGVFFI2264BQQD.yaml.gz'], read_names=['Read 1']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI9320MADV', file_set='/measurement-sets/IGVFDS5272EQRG/', illumina_read_type='R2', sequencing_run=2, lane=1, flowcell_id='AAC2JNKHV',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI9320MADV/@@download/IGVFFI9320MADV.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI2264BQQD/@@download/IGVFFI2264BQQD.yaml.gz'], read_names=['Read 2']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI0777CMJH', file_set='/measurement-sets/IGVFDS6617HQDI/', illumina_read_type='R1', sequencing_run=7395, lane=None, flowcell_id='22KFLYLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI0777CMJH/@@download/IGVFFI0777CMJH.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI5825ATCM/@@download/IGVFFI5825ATCM.yaml.gz'], read_names=['Read 1']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI7330UNCB', file_set='/measurement-sets/IGVFDS6617HQDI/', illumina_read_type='R2', sequencing_run=7395, lane=None, flowcell_id='22KFLYLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI7330UNCB/@@download/IGVFFI7330UNCB.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI5825ATCM/@@download/IGVFFI5825ATCM.yaml.gz'], read_names=['Read 2', 'Barcode index']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI1918YZDJ', file_set='/measurement-sets/IGVFDS8487QXIK/', illumina_read_type='R1', sequencing_run=7383, lane=None, flowcell_id='22KFLYLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI1918YZDJ/@@download/IGVFFI1918YZDJ.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI8012OCZQ/@@download/IGVFFI8012OCZQ.yaml.gz'], read_names=['Read 1']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI4773YQEF', file_set='/measurement-sets/IGVFDS8487QXIK/', illumina_read_type='R2', sequencing_run=7383, lane=None, flowcell_id='22KFLYLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI4773YQEF/@@download/IGVFFI4773YQEF.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI8012OCZQ/@@download/IGVFFI8012OCZQ.yaml.gz'], read_names=['Read 2', 'Barcode index']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI4117HDFW', file_set='/measurement-sets/IGVFDS8023XWXR/', illumina_read_type='R1', sequencing_run=358, lane=None, flowcell_id=None,
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI4117HDFW/@@download/IGVFFI4117HDFW.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI8568RJIP/@@download/IGVFFI8568RJIP.yaml.gz'], read_names=['Read 1']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI9675ROSO', file_set='/measurement-sets/IGVFDS8023XWXR/', illumina_read_type='R2', sequencing_run=358, lane=None, flowcell_id=None,
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI9675ROSO/@@download/IGVFFI9675ROSO.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI8568RJIP/@@download/IGVFFI8568RJIP.yaml.gz'], read_names=['Read 2']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI8200NDMC', file_set='/measurement-sets/IGVFDS1474QEZX/', illumina_read_type='R2', sequencing_run=1, lane=7, flowcell_id='22Y7YHLT3',
                                                         file_url='https://api.data.igvf.org/sequence-files/IGVFFI8200NDMC/@@download/IGVFFI8200NDMC.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI1382HZDV/@@download/IGVFFI1382HZDV.yaml.gz'], read_names=['Read 2']),
                          portal_parsing.SeqFileMetadata(file_accession='IGVFFI9019AATO', file_set='/measurement-sets/IGVFDS1474QEZX/', illumina_read_type='R1', sequencing_run=1, lane=7, flowcell_id='22Y7YHLT3', file_url='https://api.data.igvf.org/sequence-files/IGVFFI9019AATO/@@download/IGVFFI9019AATO.fastq.gz', seqspec_urls=['https://api.data.igvf.org/configuration-files/IGVFFI1382HZDV/@@download/IGVFFI1382HZDV.yaml.gz'], read_names=['Read 1'])]
