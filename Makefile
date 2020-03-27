
all:
	nextflow run main.nf \
		--sample-tsv=/Shares/layer_shared/projects/chco/exome_bakeoff/downsample_tiny.tsv \
		--run-name='kristen_downsamle_tiny' -profile fiji_hg37 -resume