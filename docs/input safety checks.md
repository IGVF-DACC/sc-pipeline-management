# Pipeine input parameters safety checks

## IGVF data portal audits

* The IGVF data portal has built-in data audits. All the relevant audits will need to be cleared before starting a pipeline run.

## Seqspec index checks

1. If there are multiple input seqspec YAMLs for the same type of assay (i.e., RNAseq or ATACseq), there will be flags if these seqspec YAMLs don't generate the same read index upon running `seqspec index`.
2. Errors will be flagged if the command `seqspec index` itself failed to execute.

## Seqspec onlist checks

1. If there are multiple input seqspec YAMLs for the same type of assay (i.e., RNAseq or ATACseq), these seqspec YAMLs must have the same onlist files listed in the `region-type: barcode, onlist: !Onlist` sections.
2. A MeasurementSet `onlist_files` property values must match that of its linked seqspec YAMLs
3. If a MeasurementSet `onlist_method` is `no combination`, then the linked seqspec YAMLs must have only one onlist file in all sections listed as `region-type: barcode, onlist: !Onlist`
    * The IGVF data portal will flag errors if a MeasurementSet `onlist_method` is `no combination` but `onlist_files` have multiple files.

```bash
# How onlist files are parsed
>$ seqspec file -m <assay_type> -s region-type -i barcode -f index -k url seqspec.yaml.gz
```
