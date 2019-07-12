#/bin/bash
set -exuo pipefail
gatk HaplotypeCaller \
        -I $file_bqsr_bam \
        --dbsnp $dbsnp \
        -O ${sample}.gvcf \
        --emit-ref-confidence GVCF \
        -R $ref_fasta
    bgzip ${sample}.gvcf
    tabix ${sample}.gvcf.gz