#/bin/bash
set -exuo pipefail
gatk GenotypeGVCFs \
    -R $ref_fasta \
    -O $out_file  \
    -V gendb://${genomics_db_workspace}