# Germline short indels Variants Analysis Pipeline
 `PLEASE NOTE THAT THIS IS VERY MUCH A WORK IN PROGRESS AND THE CONFIGURATION AND POSSIBLY THE CODE MIGHT NEED TWEAKING BEFORE RUNNING THE PIPELINE`

The pipeline is based on GATK best practices [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) and written using   [nextflow](https://www.nextflow.io).
## Pipeline Architecture
Below are the major compoents of the pipeline
- Nextflow is the main glue that wraps the individual commands and shell scripts (bwa mem for example), and orchestrate the workflow
- Configuration files (under the `conf` directory) carries the information that nextflow uses to work out which *executor* to use to run the pipeline
- Herlper utitlites (under `lib`) are some herlper funcions written to aid the main pipeline
- Individual scripts wrapped inside `nextflow processes` to carry out the actual tasks such as the *alighnemnt* or *marking duplicate reads*, or running a *structural variant caller*

## How to run the pipeline
#### Required bioinformatics packages
The current configuration on our Colorado Univeristy Fiji Cluster, we are using a singualarity container with all the required software installed in it. So if you plan to use this, you need to have `singuarlity` (current version we are using is 3.1.1, but should run with others too) on your *PATH* .

For a run on your machine (or any other infrastructure such as an HPC, or AWS), you need to add corrresponding configuration under the `conf`. The *configuration* defines a *profile* in Nextflow lingo and needs to be passed at the commandline when running the pipeline. See the `Makefile` in the top level directory for an example run. At the command prompt (while running thorugh Nextflow), you will need to pass a *samples.tsv* (tab delimited) carrying atleast three columns. The first column specifying the *sample name*, and the next two specifying the paths to the *first* and the *second reads*.
#### Example samples.tsv
The first column is the sample name, and the next two columns represents the first and the second read respectively. 


    KU1919GC	KU1919GC/KU1919GC_TD180603264_HT3LKCCXY_L3_1.fq.gz	KU1919GC/KU1919GC_TD180603264_HT3LKCCXY_L3_2.fq.gz
    KU1919P	KU1919P/KU1919P_TD180603261_HT3LKCCXY_L2_1.fq.gz	KU1919P/KU1919P_TD180603261_HT3LKCCXY_L2_2.fq.gz

#### A typical run of the pipeline
    nextflow run main.nf \
		--sample-tsv=/Shares/layer_shared/projects/cancer_center_sampled/samples.tsv \
		--results-dir=/Shares/layer_shared/projects/sequence_analysis
		--run-name='run_1' -profile fiji_hg37 -resume
	
The above run will create the results files in dir `/Shares/layer_shared/projects/sequence_analysis/run_1` *(by contcatenating the results-dir, and the run-name)*

#### Extras
#### Building singularity container
The `containers/dna_seq` subdirctory contains the required *envirnoment.yaml* and the singularity recipe file to generate the singuarlity image that we use. Our current working version of singularity is 3.1.1. Depending upon your singuarlity configuraion you will need something like the following to generate the image:

    sudo singularity build -F layer_lab_dna_seq_gatk_4.1.4.sif layer_lab_dna_seq.def




