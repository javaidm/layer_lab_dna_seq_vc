#/bin/bash
set -exuo pipefail
gatk BaseRecalibrator \
            -I  $sample_cram\
            --known-sites $DBSNP \
            --known-sites $KNOWN_INDELS \
            -O ${sample}.recal.table \
            -R $REF_FASTA