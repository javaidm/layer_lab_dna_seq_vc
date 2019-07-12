#/bin/bash
set -exuo pipefail
gatk BaseRecalibrator \
            -I  $sample_cram\
            --known-sites $dbsnp \
            --known-sites $known_indels \
            -O ${sample}.recal.table \
            -R $ref_fasta