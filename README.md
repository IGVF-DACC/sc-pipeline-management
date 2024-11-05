# sc-pipeline-management

## Current status

* Runs in a Jupyter Notebook on local machines
* Runs in command line via Python script

## Purpose

* To automate the process of generating single cell pipeline input parameters from and POSTing pipeline output data to IGVF portal

## Requirement

### 1. Python packages

* Python 3.11
* igvf-client==63.0.0
* igvf-utils==1.0.0
* firecloud==0.16.37
* seqspec==devel

### 2. Additional requirement

* Google cloud (Terra) access credentials
* IGVF portal access credentials

## Run example

### Note

* The `example_run_notebook.ipynb` contains plug-and-use codes for uploading.
* The [AnalysisSet on IGVF Sandbox](https://sandbox.igvf.org/analysis-sets/TSTDS33660419/) is an example.

### Run option 1: via Python Script

```bash
# output the input data table for Terra to a local path
python3 src/run_portal_to_terra_input_generation.py --igvf_endpoint=production --input_analysis_set=IGVFDS5316CVFR --output_dir=.

# post terra output data to IGVF portal and output a post result table to a local path
python3 src/run_terra_to_portal_data_post.py --post_endpoint=production --terra_namespace=DACC_ANVIL --terra_workspace='Playground for IGVF Single-Cell Data Processing' --terra_etype=DACC_single_cell_run_1 --upload_file=True --output_dir=.
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

# Set up IGVF query API (Sandbox or Production)
igvf_api_sandbox = api_tools.get_igvf_auth_and_api(igvf_site='sandbox')

# IGVF util (Sandbox or Production)
iu_conn_sandbox = Connection(igvf_mode='sandbox', submission=True)

# To refresh FireCloud sessions
fapi._set_session()

# Construct the input table from portal based on analysis set accessions
portal_to_terra_input_table = portal2terra_transfer.generate_pipeline_input_table(query_analysis_set_accs=igvf_sandbox_analysis_sets, igvf_api=igvf_api_sandbox)

# Very optional if want to output the table locally for inspection
portal2terra_transfer.save_pipeline_input_table(pipeline_input_table=portal_to_terra_input_table, output_dir='./')

# Upload this to Terra under the name 'DACC_Tester'
portal_to_terra_entity_type = 'DACC_single_cell_run_1'
api_tools.upload_portal_input_tsv_to_terra(terra_namespace=terra_namespace,
                                            terra_workspace=terra_workspace,
                                            terra_etype=portal_to_terra_entity_type, porta_input_table=portal_to_terra_input_table,
                                            verbose=True
                                            )
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
igvf_api_sandbox = api_tools.get_igvf_auth_and_api(igvf_site='sandbox')

# IGVF util (Sandbox or Production)
iu_conn_sandbox = Connection(igvf_mode='sandbox', submission=True)

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
                                                                    terra_etype=terra_to_portal_output_etype
                                                                   )

# Post all rows in the data table to IGVF portal
posting_all_data_report = terra2portal_transfer.add_post_status_summary_to_output_data_table(full_terra_data_table=terra_to_portal_post_datatable, igvf_api=igvf_api_sandbox, igvf_utils_api=iu_conn_sandbox, upload_file=True)

# Optional, if want to output the post results table locally
terra2portal_transfer.save_pipeline_postres_table(pipeline_postres_table=posting_all_data_report, output_dir='./')

# Upload the run results to Terra
api_tools.upload_output_post_res_to_terra(terra_namespace=terra_namespace,
                                          terra_workspace=terra_workspace,
                                          terra_etype=terra_to_portal_postres_etype,
                                          verbose=True
                                         )
```
