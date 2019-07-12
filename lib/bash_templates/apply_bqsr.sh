#/bin/bash
set -exuo pipefail
gatk ApplyBQSR \
            -I $sample_cram \
            -R $ref_fasta \
            -bqsr $recalibrated_file \
            -O ${sample}.bqsr.bam 