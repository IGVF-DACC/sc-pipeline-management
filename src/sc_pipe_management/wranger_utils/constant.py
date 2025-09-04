# Possible terms unsed in dataset descriptions and aliases
POSSIBLE_SINGLE_CELL_PIPELINE_PHRASES = ['scpipe',
                                         'single cell',
                                         'single cell uniform pipeline',
                                         'uniform pipeline',
                                         'sc-pipe',
                                         'scpipeline',
                                         'sc-uniform-pipeline',
                                         'sc-pipeline',
                                         'uniform-pipeline']

# Measurement sets will be excluded if they have any of these audits
EXCLUDED_NONCOMP_AUDITS = ['missing sequence specification',
                           'missing barcode replacement file',
                           'missing read names',
                           'missing barcode onlist']

EXCLUDED_ERROR_AUDITS = ['upload status not validated',
                         'unexpected barcode onlist',
                         'inconsistent sequence specifications',
                         'inconsistent preferred assay titles']
