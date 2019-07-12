all: run
run:
	# nextflow -log logs/log run main.nf -profile fiji
	sbatch chco_v1_on_fiji.sh
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
