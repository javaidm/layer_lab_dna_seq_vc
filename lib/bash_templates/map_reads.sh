#/bin/bash
set -exuo pipefail
bwa mem \
-t $THREADS -R "@RG\tID:$sample\tSM:$sample\tPL:$PLATFORM\tPU:$sample\tLB:$sample" $REF_FASTA $r1 $r2 \
| samblaster \
| samtools sort --output-fmt-option seqs_per_slice=4000 -O CRAM --reference $REF_FASTA -m 18G -@ 6 /dev/stdin -o ${sample}.cram