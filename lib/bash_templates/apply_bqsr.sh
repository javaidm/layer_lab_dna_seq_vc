#/bin/bash
set -exuo pipefail
gatk ApplyBQSR \
            -I $sample_cram \
            -R $REF_FASTA \
            -bqsr $recalibrated_file \
            -O ${sample}.bqsr.bam 