#/bin/bash
set -exuo pipefail
bcftools concat -o ${out_file} -Ov \
chr1.vcf \
chr2.vcf \
chr3.vcf \
chr4.vcf \
chr5.vcf \
chr6.vcf \
chr7.vcf \
chr8.vcf \
chr9.vcf \
chr10.vcf \
chr11.vcf \
chr12.vcf \
chr13.vcf \
chr14.vcf \
chr15.vcf \
chr16.vcf \
chr17.vcf \
chr18.vcf \
chr19.vcf \
chr20.vcf \
chr21.vcf \
chr22.vcf \
chrX.vcf \
chrY.vcf \
chrMT.vcf 
