# Germline short indels Variants Analysis Pipeline
 `PLEASE NOTE THAT THIS IS VERY MUCH A WORK IN PROGRESS AND THE CONFIGURATION AND POSSIBLY THE CODE MIGHT NEED TWEAKING BEFORE RUNNING THE PIPELINE`

The pipeline is based on GATK best practices [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) and written using   [nextflow](https://www.nextflow.io).
## Pipeline Architecture
Below are the major components of the pipeline
- Nextflow is the main glue that wraps the individual commands and shell scripts (bwa mem for example), and orchestrate the workflow
- Configuration files (under the `conf` directory) carries the information that nextflow uses to work out which *executor* to use to run the pipeline
- Helper utilities (under `lib`) are some helper functions written to aid the main pipeline
- Individual scripts wrapped inside `nextflow processes` to carry out the actual tasks such as the *alignment* or *marking duplicate reads*, or running a *structural variant caller*
- A Singularity container is our preferred way to run the pipeline. Currently we include the singularity configuration in the subdirectory *containers/dna_seq*

## How to run the pipeline
#### Required bioinformatics packages
In the current configuration on our Colorado University Boulder Fiji Cluster, we are using a singularity container with all the required software installed in it. So if you plan to use this, you need to have `singularity` (current version we are using is 3.1.1, but should run with others too) on your Unix *PATH* .

For a run on your machine (or any other infrastructure such as an HPC, or AWS), you need to add corresponding configuration under the `conf`. The *configuration* defines a *profile* in Nextflow lingo and needs to be passed at the commandline when running the pipeline. See the `Makefile` in the top level directory for an example run. At the command prompt (while running through Nextflow), you will need to pass a *samples.tsv* (tab delimited) carrying at least three columns. The first column specifying the *sample name*, and the next two specifying the paths to the *first* and the *second reads* respectively.
#### Example samples.tsv
The first column is the sample name, and the next two columns represents the first and the second read respectively. 


    KU1919GC	KU1919GC/KU1919GC_TD180603264_HT3LKCCXY_L3_1.fq.gz	KU1919GC/KU1919GC_TD180603264_HT3LKCCXY_L3_2.fq.gz
    KU1919P	KU1919P/KU1919P_TD180603261_HT3LKCCXY_L2_1.fq.gz	KU1919P/KU1919P_TD180603261_HT3LKCCXY_L2_2.fq.gz

#### A typical run of the pipeline
    nextflow run main.nf \
		--sample-tsv=/Shares/layer_shared/projects/cancer_center_sampled/samples.tsv \
		--results-dir=/Shares/layer_shared/projects/sequence_analysis
		--run-name='run_1' -profile fiji_hg37 -resume
	
The above run will create the results files in dir `/Shares/layer_shared/projects/sequence_analysis/run_1` *(by concatenating the results-dir, and the run-name)*
#### Passing parameters to the pipeline
Nextflow (and this pipeline) allows the following three ways to pass params to a pipeline run:
1. At the *commandline*, using syntex like *- -sample-tsv*, or *- -results-dir* (note the double dashes before param names)
2. Still at the commandline but putting params in a *yaml* file and by passing that file as *- -params my_params.yaml*. See below  an excerpt from a params.yaml file:

        sampleTsv: "/Shares/layer_shared/projects/cancer_center_tiny/samples.tsv"
        runName: "cancer_center_tiny"
        conditionsTsv: "/Shares/layer_shared/projects/cancer_center_tiny/conditions.tsv"
        sendEmail: False
        logFilesDur:    "logs"
        nfOptions: "-with-reports -with-trace -with-timeline -with-dag flowchart.png"
3. As part of the *profile config* file. See examples in the *conf* subdirectory.

#### Extras
#### Building singularity container
The `containers/dna_seq` subdirectory contains the required *environment.yaml* and the singularity recipe file to generate the singularity image that we use. Our current working version of singularity is 3.1.1. Depending upon your singularity configuration you will need something like the following to generate the image:

    sudo singularity build -F layer_lab_dna_seq_gatk_4.1.4.sif layer_lab_dna_seq.def




