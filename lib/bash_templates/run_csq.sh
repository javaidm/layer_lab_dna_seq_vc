#/bin/bash
set -exuo pipefail
bcftools norm -m- -Ov  -f $ref_fasta -w 10000 $vcf |\
bcftools csq -s - --ncsq 40 -g $gff3 -l -f $ref_fasta -Ob -o $out_file