# sc-pipeline-management

## Current status

* Runs in command line via Python script.
  * **Pipeline input**: `run_portal_to_terra_input_generation.py` script will perform the following tasks:
    1. Fetch, generate, and output all the input parameters using a placeholder Analysis Set accession as a TSV file.
    2. Automatically upload this TSV file onto Terra workspace.
    3. Automatically upload the generated final barcode inclusion lists onto Terra workspace disk.
    4. **NOTE**: The starting of the workflow is still manually triggered at this moment on Terra.
  * **IGVF accessioning**: `run_terra_to_portal_data_post.py` script will perform the following taskss:
    1. Post and create the data objects for RNAseq data, ATACseq data, and QC metric data to the same placeholder Analysis Set that generated the input parameters.
    2. Upload the pipeline output files to the associated data objects.
* Runs in a Jupyter Notebook.

## Purpose

* To automate the process of generating single cell pipeline input parameters from and POSTing pipeline output data to IGVF portal
* More details on the codes are under the `docs/` folder.

## Requirement

### 1. Key Python packages

See `requirements.txt` for package details.

* Python 3.11
* igvf-client
* igvf-utils
* firecloud (0.16.37)
* seqspec

#### Installation

* All packages can be install directly from the requirements.txt using `pip install -r requirements.txt`
* If FireCloud installation cannot complete using `requirements.txt`, it can be install manually. Download the either the zip or tar.gz source file from [FireCloud release v0.16.37](https://github.com/broadinstitute/fiss/releases/tag/v0.16.37) and unzip. Run `python setup.py install` from the same directory.

### 2. Additional requirement

* **Set up Google cloud (Terra) access credentials.** FireCloud uses the [Google Cloud SDK](https://cloud.google.com/sdk/) to manage authorization. To use the firecloud CLI or API, you must install the SDK and login locally with `gcloud auth login`.
  * This is needed for 1) running input parameters generation, and 2) accessioning to IGVF portal.
  * Install Google CLI following this link: <https://cloud.google.com/sdk/docs/install-sdk>.

* **IGVF portal access credentials.** Generate API credentials under the IGVF profiles page after logging into the portal. Store the public API key and secret API key as environmental variables. See [IGVF API credential configration guide](https://github.com/IGVF-DACC/igvf_utils/wiki/Configuration) for more details.
  * This is needed for **all** the scripts in this repo.

## The input parameters that the codes will generate

1. Read 1 FASTQ file accessions and URLs for ATACseq and RNAseq
2. Read 2 FASTQ file accessions and URLs for ATACseq and RNAseq
3. Barcode FASTQ file accessions and URLs for ATACseq and RNAseq (optional for RNAseq)
4. Seqspec YAML file URLs
5. Read format index for ATACseq and RNAseq
6. Final barcode onlist for ATACseq and RNAseq
7. kb_strand (forward only)
8. onlist_mapping boolean value
9. Genome index file URL
10. Transcriptome index file URL
11. Genome FASTA file URL
12. Genome ref file URL (fallback option if no reference file indicated)
13. IGVF portal Analysis Set accession
14. Subpool ID (associated IGVF sample accessions)

## Pipeline output and submission to IGVF portal

* The overall data posted and relations to other objects on the IGVF portal is detailed in the [uniform pipeline workflow on IGVF portal](https://data.igvf.org/workflows/IGVFWF6403DVII/).
* The main outputs for different assay types are as follows:
  * ATACseq
    * Alignment (bam)
    * Bam index (bai)
    * Fragment files (tsv/bed3+)
    * Fragment index (tbi)
    * A set of QC metrics
  * RNAseq
    * Sparse gene count matrix (H5AD)
    * Kallisto output tarball (tar)
    * A set of QC metrics

<!-- markdownlint-disable MD033 -->
<img src="https://api.data.igvf.org/documents/6f250d2a-aff6-4ae5-ab71-9994a1bb584c/@@download/attachment/sc-pipeline_RNAseq_data.png" alt="RNAseq" width="500">

<img src="https://api.data.igvf.org/documents/9273ea79-5120-41c9-85b4-f2e76caf65e8/@@download/attachment/sc-pipeline_ATACseq_data.png" alt="ATACseq" width="500">
<!-- markdownlint-enable MD033 -->

## Generate input parameters for Terra pipeline run

```python
# Output the input data table for Terra to a local path using analysis set IDs from an input string
python3 src/sc_pipe_management/run_portal_to_terra_input_generation.py --igvf_endpoint=prod --input_analysis_set=IGVFDS5316CVFR --terra_etype=single_cell_run

# Output the input data table for Terra to a local path using a txt file with one analysis file per line
python3 src/sc_pipe_management/run_portal_to_terra_input_generation.py --igvf_endpoint=prod --input_analysis_set_file=/local_dirs/analysis_set_accessions.txt --terra_etype=single_cell_run
```

## Accession new data to IGVF portal

```python
# Post terra output data to IGVF portal and output a post result table to default path
python3 src/sc_pipe_management/run_terra_to_portal_data_post.py --post_endpoint=prod --terra_etype=single_cell_run --upload_file=True

# Post output data with some runs excluded
python3 src/sc_pipe_management/run_terra_to_portal_data_post.py --post_endpoint=prod --terra_etype=single_cell_run --excluded_accs=/local_dirs/excluded_analysis_set_accessions.txt --upload_file

# Resume posting output data after an unexpected pause or unposted data errors
python3 src/sc_pipe_management/run_terra_to_portal_data_post.py --post_endpoint=prod --terra_etype=single_cell_run --upload_file=True --resumed_posting
```

## Utils tools

### Create analysis sets to start pipeline runs

* The src/JupyterNotebook/analysis_set_setup.ipynb notebook contains plug-and-use templates and other utils to create new analysis sets based on a list of measurement sets.

### Quality check posted results

* This tool helps data wranglers, submitters, and other users who need to verify accessioning results in an automated fashion.

#### Run option 1: Python script

```python
# If Running against a list of accessions
python3 src/sc_pipe_management/wrangler_utils/check_accession_results.py --igvf_endpoint prod --qc_analysis_set_file /local_dir/analysis_set_accessions.txt --output_file_path /local_dir/qa_results.json

# If Running against one analysis set
python3 src/sc_pipe_management/wrangler_utils/check_accession_results.py --igvf_endpoint prod --qc_analysis_sets IGVFDS5316CVFR --output_file_path /local_dir/qa_results.json
```

#### Run option 2: Jupyter notebook

* The src/JupyterNotebook/qa_script.ipynb notebook contains plug-and-use templates to QA posted datasets.
