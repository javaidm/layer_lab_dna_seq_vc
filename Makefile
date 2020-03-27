
all:
	nextflow run main.nf \
		--sample-tsv=/Shares/layer_shared/projects/cancer_center_sampled/samples.tsv \
		--run-name='cancer_center_sampled_germline_hg37' -profile fiji_hg37 -resume