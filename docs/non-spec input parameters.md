# Generating non-seqspec related input parameters

## Read fastq file accessions and download urls

The fastq read files will be retrieved from the IGVF data portal based on the linked AnalysisSet and input MeasurementSet. The files will be sorted in the order of `illumina_read_type`, `sequencing run`, then `lane`, whichever is available. If an AnalysisSet has multiple input MeasurementSet of the same assay type (i.e., RNAseq or ATACseq), then the sorting order will be `file_set`, `illumina_read_type`, `sequencing run`, then `lane`.

For ATACseq, if Read2 and Barcode Index are in the same fastq file, then this file will be passed to both `atac_read2_accessions` and `atac_barcode_accessions`.

For RNAseq, `rna_barcode_accessions` will only be populated if there is a separate Barcode Index fastq file. No file accession will be passed if no Barcode Index is specified or Read2 and Barcode Index are combined into the same fastq file.

## Reference files

Genome fasta, and genome and transcriptome index files are used as input references. These files are available on IGVF portal.

### Human

* Genome fasta: <https://data.igvf.org/reference-files/IGVFFI0653VCGH/>
* Genome index: <https://data.igvf.org/reference-files/IGVFFI7969JLFC/>
* Transcriptome index: <https://data.igvf.org/reference-files/IGVFFI9561BASO/>

### Mouse

* Genome fasta: <https://data.igvf.org/reference-files/IGVFFI9282QLXO/>
* Genome index: <https://data.igvf.org/reference-files/IGVFFI5593VLWB/>
* Transcriptome index: <https://data.igvf.org/reference-files/IGVFFI5078MNED/>

## kb strand

The current `kb strand` value is hard coded in as `forward`.

## Whether to create onlist mapping

The boolean value for whether to create onlist mapping between RNAseq and ATACseq is determined by the input MeasurementSet `preferred_assay_title`. If the `preferred_assay_title` is in the following list, the boolean value will be `True`.

Currently accepted `preferred_assay_titles`:

* 10x multiome
* 10x multiome with MULTI-seq

## Subpool ID

The `subpool ID` uses the input samples' accessions. If there are multiple input samples, the sample IDs will be joined by hyphens `-`.

## Prefix

The prefix is added onto output file names. The current `prefix` value is the same as the input AnalysisSet accession.
