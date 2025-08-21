import sys
import os

# Get the absolute path to the 'src' directory (one level up from the 'tests' directory)
src_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the 'src' directory to sys.path if it's not already there
if src_path not in sys.path:
    sys.path.insert(0, src_path)

import sc_pipe_management.input_params_generation.seqspec_parsing as seqspec_parsing


# Seqspec files for different assays (Splitseq is not released, needs to download when running the test)
SEQSPEC_FILES_BY_ASSAY_TITLES = {
    '10x multiome': {'atac': 'src/tests/test_files/IGVFFI1365BOUI.yaml.gz',
                     'rna': 'src/tests/test_files/IGVFFI2433YUEJ.yaml.gz'},  # 10x multiome RNA seqspec
    'parse splitseq': {'atac': None,
                       'rna': 'src/tests/test_files/IGVFFI2264BQQD.yaml.gz'},  # SPLiT-seq RNA seqspec (read_id has .fastq.gz suffix)
    'shareseq': {'atac': 'src/tests/test_files/IGVFFI8012OCZQ.yaml.gz',  # ShareSeq does not have atac seqspec
                 'rna': 'src/tests/test_files/IGVFFI5825ATCM.yaml.gz'},  # ShareSeq RNA seqspec
    '10x multiome double': {'atac': None,  # No need to test ATAC on double cDNA reads
                            'rna': 'src/tests/test_files/IGVFFI8568RJIP.yaml.gz'},  # This seqspec has cDNA on both R1 and R2
    'rnaseq with bad read_id': {'atac': None,  # ShareSeq does not have atac seqspec
                                'rna': 'src/tests/test_files/IGVFFI1382HZDV.yaml.gz'}  # This seqspec has a bad read_id
}

SEQSPEC_FILES_BY_ASSAY_TITLES = {
    # Standard 10x multiome seqspecs
    '10x multiome': {'atac': seqspec_parsing.SeqSpecMetadata(
        seqspec_file_path='src/tests/test_files/IGVFFI1365BOUI.yaml.gz',
        modality='atac',
        ordered_read_ids=["IGVFFI4665EVGC",
                          "IGVFFI4986VMDU", "IGVFFI0391NHGA"],
        onlist_files=[
            'https://api.data.igvf.org/tabular-files/IGVFFI7587TJLC/@@download/IGVFFI7587TJLC.tsv.gz']
    ),
        'rna': seqspec_parsing.SeqSpecMetadata(
        seqspec_file_path='src/tests/test_files/IGVFFI2433YUEJ.yaml.gz',
        modality='rna',
        ordered_read_ids=[
            "IGVFFI6948UZJO", "IGVFFI7951DAQB"],
        onlist_files=[
            'https://api.data.igvf.org/tabular-files/IGVFFI8751YQRY/@@download/IGVFFI8751YQRY.tsv.gz']
    )},
    # SPLiT-seq seqspec with read_id having .fastq.gz suffix
    'parse splitseq': {'atac': None,
                       'rna': seqspec_parsing.SeqSpecMetadata(
                           seqspec_file_path='src/tests/test_files/IGVFFI2264BQQD.yaml.gz',
                           modality='rna',
                           ordered_read_ids=["IGVFFI4633WXZI.fastq.gz",
                                             "IGVFFI9320MADV.fastq.gz"],
                           onlist_files=[
                               "https://api.data.igvf.org/tabular-files/IGVFFI0924TKJO/@@download/IGVFFI0924TKJO.tsv.gz",
                               "https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz",
                               "https://api.data.igvf.org/tabular-files/IGVFFI1138MCVX/@@download/IGVFFI1138MCVX.tsv.gz"
                           ]
                       )},
    # ShareSeq seqspecs
    'shareseq': {'atac': seqspec_parsing.SeqSpecMetadata(
        seqspec_file_path='src/tests/test_files/IGVFFI8012OCZQ.yaml.gz',
        modality='atac',
        ordered_read_ids=["IGVFFI1918YZDJ.fastq.gz",
                          "IGVFFI4773YQEF.fastq.gz"],
        onlist_files=[
            "https://api.data.igvf.org/tabular-files/IGVFFI1918YZDJ/@@download/IGVFFI1918YZDJ.tsv.gz",
            "https://api.data.igvf.org/tabular-files/IGVFFI4773YQEF/@@download/IGVFFI4773YQEF.tsv.gz"]
    ),
        'rna': 'src/tests/test_files/IGVFFI5825ATCM.yaml.gz'},  # ShareSeq RNA seqspec
    '10x multiome double': {'atac': None,  # No need to test ATAC on double cDNA reads
                            'rna': 'src/tests/test_files/IGVFFI8568RJIP.yaml.gz'},  # This seqspec has cDNA on both R1 and R2
    'rnaseq with bad read_id': {'atac': None,  # ShareSeq does not have atac seqspec
                                'rna': 'src/tests/test_files/IGVFFI1382HZDV.yaml.gz'}  # This seqspec has a bad read_id
}

# Expected outputs from seqspec
SEQSPEC_OUTPUT_BY_ASSAY_TITLES = {
    '10x multiome': {
        'atac': {'read_index': 'bc:8:23:-,r1:0:49,r2:0:49',
                 'index_input': 'IGVFFI4665EVGC,IGVFFI4986VMDU,IGVFFI0391NHGA',
                 'onlist_input': 'IGVFFI0391NHGA',
                 'onlist_final_list': '59ed4e3b7c64af921bde8b8497924c02',
                 'onlist_method': 'no combination'
                 },
        'rna': {'read_index': '0,0,16:0,16,28:1,0,90',
                'index_input': 'IGVFFI6948UZJO,IGVFFI7951DAQB',
                'onlist_input': 'IGVFFI7951DAQB',
                'onlist_final_list': '5b223bbf2d45ba5e8bb55787872321c8',
                'onlist_method': 'no combination'
                }
    },
    'parse splitseq': {
        'atac': None,
        'rna': {'read_index': '1,10,18,1,48,56,1,78,86:1,0,10:0,0,140',
                'index_input': 'IGVFFI4633WXZI.fastq.gz,IGVFFI9320MADV.fastq.gz',
                'onlist_input': 'IGVFFI9320MADV.fastq.gz',
                'onlist_final_list': '7f1e51743561883977e5a905bb345ac8',
                'onlist_method': 'product'
                }
    },
    '10x multiome double': {
        'atac': None,
        'rna': {'read_index': '0,0,16:0,16,28:1,0,101',
                'index_input': 'IGVFFI4117HDFW,IGVFFI9675ROSO',
                'onlist_input': 'IGVFFI9675ROSO',
                'onlist_final_list': '5b223bbf2d45ba5e8bb55787872321c8',
                'onlist_method': 'no combination'
                }
    },
    'shareseq': {
        'atac': {'read_index': 'bc:115:122,bc:153:160,bc:191:198,r1:0:99,r2:0:99',
                 'index_input': 'IGVFFI1918YZDJ.fastq.gz,IGVFFI4773YQEF.fastq.gz',
                 'onlist_input': 'IGVFFI4773YQEF.fastq.gz',
                 'onlist_final_list': 'f7eabd710222e594459698a6c624de7b',
                 'onlist_method': 'product'
                 },
        'rna': {'read_index': '1,115,123,1,153,161,1,191,199:1,0,10:0,0,100',
                'index_input': 'IGVFFI0777CMJH.fastq.gz,IGVFFI7330UNCB.fastq.gz',
                'onlist_input': 'IGVFFI7330UNCB.fastq.gz',
                'onlist_final_list': 'f7eabd710222e594459698a6c624de7b',
                'onlist_method': 'product'
                },
    },
    'rnaseq with bad read_id': {
        'atac': None,  # ShareSeq does not have atac seqspec
        'rna': {'read_index': '0,0,16:0,16,28:1,0,90',
                'index_input': 'RNA Read 1,RNA Read 2',
                'onlist_input': 'RNA Read 2',
                'onlist_final_list': '5b223bbf2d45ba5e8bb55787872321c8',
                'onlist_method': 'no combination'
                }
    }
}
