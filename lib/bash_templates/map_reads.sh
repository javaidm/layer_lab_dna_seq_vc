#/bin/bash
set -exuo pipefail
#echo "STARTING bwa cpus:$threads sample:$sample reference:$ref_fasta r1:$r1 r2:$r2"
bwa mem \
-t $threads -R "@RG\tID:$sample\tSM:$sample\tPL:$platform\tPU:$sample\tLB:$sample" $ref_fasta $r1 $r2 \
| samblaster \
| samtools sort --output-fmt-option seqs_per_slice=4000 -O CRAM --reference $ref_fasta -m 18G -@ 6 /dev/stdin -o ${sample}.cram