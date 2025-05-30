# Parsing Seqspec YAMLs for input parameters

## Seqspec index

Both RNAseq and ATACseq read format index will be generated by passing in a list of seqspec `read_id` of a specific order.

The general rule is that if Read1, Read2, and Barcode Index are 3 separate fastq file, all 3 `read_id` are used. However, if Read2 and Barcode Index are combined into the same fastq file, only Read2 `read_id` will be used and Barcode Index `read_id` will be disregarded.

### ATACseq read format

```bash
# Read1 fastq, Read 2 fastq, and Barcode index fastq are 3 separate files
>$ seqspec index -m atac -t chromap -s read -i Read1_read-id,Read2_read-id,BarcodeIndex_read-id seqspec.yaml.gz

# Read1 is one fastq file, and Read2 + Barcode Index are combined into the same fastq file.
>$ seqspec index -m atac -t chromap -s read -i Read1_read-id,Read2_read-id seqspec.yaml.gz
```

### RNAseq read format

```bash
# Read1 fastq, Read 2 fastq, and Barcode index fastq are 3 separate files
>$ seqspec index -m rna -t kb -s read -i Read1_read-id,Read2_read-id,BarcodeIndex_read-id seqspec.yaml.gz

# Read1 is one fastq file. Barcode Index is either combined into Read2 fastq or not specified.
>$ seqspec index -m rna -t kb -s read -i Read1_read-id,Read2_read-id seqspec.yaml.gz
```

## Seqspec onlist

If there will be a barcode onlist combinatorial step, Both RNAseq and ATACseq read format index will be generated by passing in either Barcode Index `read-id` or Read2 `read-id` depending on whether Barcode Index is a separate fastq file, combined with Read2 fastq, or not specified.

If there is no combinatorial step, a standard `seqspec onlist` command will be passed.

Currently, the only combinatorial style is `-f product` regardless of what the Measurement Set `onlist_method` value is.

### ATACseq onlist

```bash
# If no barcode combinatorial step
>$ seqspec onlist -m atac -s region-type -i barcode -o output_list.txt seqspec.yaml.gz

# If combinatorial, with Barcode Index as a separate fastq file
>$ seqspec onlist -m atac -f product -s read -i BarcodeIndex_read-id -o output_list.txt seqspec.yaml.gz

# If combinatorial, with Read2 and Barcode Index in the same fastq file
>$ seqspec onlist -m atac -f product -s read -i Read2_read-id -o output_list.txt seqspec.yaml.gz
```

### RNAseq onlist

```bash
# If no barcode combinatorial step
>$ seqspec onlist -m rna -s region-type -i barcode -o output_list.txt seqspec.yaml.gz

# If combinatorial, with Barcode Index as a separate fastq file
>$ seqspec onlist -m rna -f product -s read -i BarcodeIndex_read-id -o output_list.txt seqspec.yaml.gz

# If combinatorial, with a fastq file that is Read2+Barcode Index OR no Barcode Index specified
>$ seqspec onlist -m rna -f product -s read -i Read2_read-id -o output_list.txt seqspec.yaml.gz
```
