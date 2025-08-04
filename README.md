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
* More details on the codes are under the `doc` folder.

## Requirement

### 1. Key Python packages

See `requirements.txt` for package details.

* Python 3.11
* igvf-client
* igvf-utils
* firecloud
* seqspec

#### Installation

* The igvf-client, igvf-utils, and seqspec packages can be install directly from the requirements.txt using `pip install -r requirements.txt`
* The firecloud package needs to be installed manually due to setuptools incompatibilities. Download the either the zip or tar.gz source file from [FireCloud release v0.16.37](https://github.com/broadinstitute/fiss/releases/tag/v0.16.37) and unzip. Run `python setup.py install`.

### 2. Additional requirement

* **Set up Google cloud (Terra) access credentials.** FireCloud uses the [Google Cloud SDK](https://cloud.google.com/sdk/) to manage authorization. To use the firecloud CLI or API, you must install the SDK and login locally with `gcloud auth login`.
  * Install Google CLI following this link: <https://cloud.google.com/sdk/docs/install-sdk>.
* **IGVF portal access credentials.** Generate API credentials under the IGVF profiles page after logging into the portal. Store the public API key and secret API key as environmental variables. See [IGVF API credential configration guide](https://github.com/IGVF-DACC/igvf_utils/wiki/Configuration) for more details.

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

## Run example

### Note

* The `example_run_notebook.ipynb` contains plug-in-and-use codes for generating input data tables and posting output data to the IGVF data portal.
* The [AnalysisSet on IGVF Sandbox](https://sandbox.igvf.org/analysis-sets/TSTDS33660419/) is an example.

### Run option 1: via Python Script

```bash
# Output the input data table for Terra to a local path using analysis set IDs from an input string
python3 src/run_portal_to_terra_input_generation.py --igvf_endpoint=prod --input_analysis_set=IGVFDS5316CVFR --terra_etype=pipeline_test_run

# Output the input data table for Terra to a local path using a txt file with one analysis file per line
python3 src/run_portal_to_terra_input_generation.py --igvf_endpoint=prod --input_analysis_set_file=/local_dirs/analysis_set_accessions.txt --terra_etype=pipeline_test_run

# post terra output data to IGVF portal and output a post result table to a local path
python3 src/run_terra_to_portal_data_post.py --post_endpoint=prod --terra_namespace=DACC_ANVIL --terra_workspace='Playground for IGVF Single-Cell Data Processing' --terra_etype=DACC_single_cell_run_1 --upload_file=True --output_dir="$(pwd)/terra_input_datatables"
```

### Run option 2: via Jupyter Notebook

#### Step 1. Construct and upload pipeline input parameters from IGVF portal

* Requires an Analysis Set to be created in advance.

```python
# Set up IGVF Python client, IGVF utils, and FireCloud credentials
import igvf_and_terra_api_tools as api_tools
import portal_to_terra_input_from_anaset as portal2terra_transfer
import terra_to_portal_posting as terra2portal_transfer

from igvf_utils.connection import Connection
import firecloud.api as fapi

# Terra workspace info
terra_namespace = 'DACC_ANVIL'
terra_workspace = 'Playground for IGVF Single-Cell Data Processing'

# Analysis Set for single cell pipeline run
igvf_sandbox_analysis_sets = ['TSTDS33660419', 'TSTDS12024147']
local_final_barcode_dir = os.path.join(os.getcwd(), 'final_barcode_list/')
gs_barcode_list_bucket = 'gs://unittest_mock_bucket/submissions/final_barcode_onlist/'

# Set up IGVF query API (Sandbox or Production)
igvf_endpoint = 'sandbox'
igvf_api_keys = api_tools.set_up_api_keys(igvf_endpoint=igvf_endpoint)
igvf_api_sandbox = api_tools.get_igvf_client_auth(igvf_site=igvf_endpoint,
                                                  igvf_api_keys=igvf_keys)

# IGVF util (Sandbox or Production)
iu_conn_sandbox = api_tools.get_igvf_utils_connection(igvf_api_keys=igvf_api_keys,
                                                      igvf_utils_mode=igvf_endpoint,
                                                      submission_mode=True)

# To refresh FireCloud sessions
fapi._set_session()

# Construct the input table from portal based on analysis set accessions
portal_to_terra_input_table = portal2terra_transfer.generate_pipeline_input_table(query_analysis_set_accs=unittest_analysis_sets,
                                                                                  igvf_api=igvf_api_sandbox,
                                                                                  terra_etype='unittest_pipeline_tester',
                                                                                  local_barcode_file_dir=local_final_barcode_dir,
                                                                                  gs_barcode_list_bucket=gs_barcode_list_bucket
                                                                                  )

# Very optional if want to output the table locally for inspection
portal2terra_transfer.save_pipeline_input_table(pipeline_input_table=portal_to_terra_input_table, output_dir='./')

# Upload this to Terra under the name 'DACC_Tester'
portal_to_terra_entity_type = 'DACC_single_cell_run_1'
api_tools.upload_portal_input_tsv_to_terra(terra_namespace=terra_namespace,
                                            terra_workspace=terra_workspace,
                                            terra_etype=portal_to_terra_entity_type,
                                            porta_input_table=portal_to_terra_input_table,
                                            verbose=True)
```

#### Step 2. Post pipeline run results to IGVF portal

* This process will also upload a POST result report table to Terra.

```python
# Set up IGVF Python client, IGVF utils, and FireCloud credentials
import igvf_and_terra_api_tools as api_tools
import terra_to_portal_posting as terra2portal_transfer

from igvf_utils.connection import Connection
import firecloud.api as fapi

# Set up IGVF query API (Sandbox or Production)
igvf_endpoint = 'sandbox'
# Get local API keys
igvf_api_keys = api_tools.set_up_api_keys(igvf_endpoint=igvf_endpoint)
# Get IGVF client API
igvf_api_sandbox = api_tools.get_igvf_client_auth(igvf_site=igvf_endpoint,
                                                  igvf_api_keys=igvf_keys)
# IGVF util (Sandbox or Production)
iu_conn_sandbox = api_tools.get_igvf_utils_connection(igvf_api_keys=igvf_api_keys,
                                                      igvf_utils_mode=igvf_endpoint,
                                                      submission_mode=True)

# To refresh FireCloud sessions
fapi._set_session()

# Terra workspace info
terra_namespace = 'DACC_ANVIL'
terra_workspace = 'Playground for IGVF Single-Cell Data Processing'

# Set up params for retrieving post-pipeline data table from Terra
terra_to_portal_output_etype = 'DACC_single_cell_run_1'
terra_to_portal_postres_etype = 'DACC_single_cell_run_1_postres'

# Generate the input data table for Terra
terra_to_portal_post_datatable = api_tools.get_terra_tsv_data_table(terra_namespace=terra_namespace,
                                                                    terra_workspace=terra_workspace,
                                                                    terra_etype=terra_to_portal_output_etype)

# Post all rows in the data table to IGVF portal
# Run all posting jobs
terra_to_portal_post_runs = terra2portal_transfer.post_all_successful_runs(igvf_api=igvf_api_sandbox, igvf_utils_api=iu_conn_sandbox, upload_file=False, full_terra_data_table=terra_to_portal_post_datatable)

# Get the full post success summary
terra_to_portal_post_summary = terra2portal_transfer.summarize_post_status(post_results=terra_to_portal_post_runs)

# Update the original output table with brief post summary
terra_to_portal_post_summary = terra2portal_transfer.add_post_status_summary_to_output_data_table(full_terra_data_table=terra_to_portal_post_datatable, post_status_df=terra_to_portal_post_summary)

# Optional, if want to output the post results table locally
terra2portal_transfer.save_pipeline_postres_table(pipeline_postres_table=posting_all_data_report, output_dir='./')

# Optional, upload the run results to Terra
api_tools.upload_output_post_res_to_terra(terra_namespace=terra_namespace,
                                          terra_workspace=terra_workspace,
                                          terra_etype=terra_to_portal_postres_etype,
                                          verbose=True
                                         )
```
