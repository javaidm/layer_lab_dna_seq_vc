#/bin/bash
set -exuo pipefail
# Prepare sample map

for x in `ls *.gvcf.gz`
do
    sample=`echo \$x | cut -d '.' -f1`
    echo "\${sample}\t\${x}" >> ${sample_map}
done
gatk GenomicsDBImport \
--genomicsdb-workspace-path ${gDB} \
-L ${chr} \
--sample-name-map ${sample_map} \
--reader-threads 8