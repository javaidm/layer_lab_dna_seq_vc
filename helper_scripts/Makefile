all: run
run:
	nextflow run javaidm/layer_lab_dna_seq_vc --run-name kristen_toy_samples --reads "/scratch/Shares/layer/nextflow/kristen/three_sample_small_fastq_files/**/*_R{1,2}_*.fastq.gz" -profile fiji
	#nextflow -log logs/log run main.nf --run-name kristen_toy_samples --reads "/scratch/Shares/layer/nextflow/kristen/three_sample_small_fastq_files/**/*_R{1,2}_*.fastq.gz" -profile fiji
	#nextflow -log logs/log run main.nf --sample "../kristen_samples.tsv" -profile fiji
	
	#sbatch chco_v1_on_fiji.sh
clean:
	rm -fr *~
	rm -fr .nextflow
	rm -fr work
	rm -fr logs
	rm -fr *.out
	rm -fr *.err
	rm -fr nf_reports
cancel:
	scancel -u ${USER}	
