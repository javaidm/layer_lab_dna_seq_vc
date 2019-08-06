#/bin/bash
set -exuo pipefail
gatk HaplotypeCaller \
        -I $file_bqsr_bam \
        --dbsnp $DBSNP \
        -O ${sample}.gvcf \
        --emit-ref-confidence GVCF \
        -R $REF_FASTA
    bgzip ${sample}.gvcf
    tabix ${sample}.gvcf.gz