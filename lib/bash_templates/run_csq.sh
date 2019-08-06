#/bin/bash
set -exuo pipefail
bcftools norm -m- -Ov  -f $REF_FASTA -w 10000 $vcf |\
bcftools csq -s - --ncsq 40 -g $ENSEMBL_GENE_ANNOTATIONS -l -f $REF_FASTA -Ob -o $out_file