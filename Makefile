
all:
	nextflow run main.nf \
		--sample-tsv=/Shares/layer_shared/projects/cancer_center_merged/samples.tsv \
		--conditions-tsv /Shares/layer_shared/projects/cancer_center_merged/conditions.tsv \
		--run-name='cancer_center_hg38' -profile fiji -resume