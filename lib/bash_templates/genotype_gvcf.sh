#/bin/bash
set -exuo pipefail
gatk GenotypeGVCFs \
    -R $REF_FASTA \
    -O $out_file  \
    -V gendb://${genomics_db_workspace}