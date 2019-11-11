#!/usr/bin/env nextflow
nextflow.preview.dsl=2
initParamsToDefaults()
if (params.help) exit 0, helpMessage()
platform = "ILLUMINA"
templateDir = "${workflow.projectDir}/lib/bash_templates"
// processParams()
/* Process the parameters and set the environemnt */
params.name = 'Layer Lab DNA Seq Analysis Pipeline'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>
params.contactMail ?: 'mahmood.javaid@colorado.edu'
if (!params.sample && !params.reads){
    exit 1, 'You must sepcify either a sample tsv file or the reads, see --help'
}
// Check parameters for profile fiji
// If they are not provided at the command line, set them for fiji from the env section of the config
if (workflow.profile == 'fiji'){
  resultsDir = params.resultsDir ?: "$fiji_results_dir"
  dbsnp = params.dbsnp ?: "$fiji_dbsnp"
  ensemblGeneAnnotation = params.ensemblGeneAnnotation ?: "$fiji_ensembl_gene_annotation"
  knownIndels = params.knownIndels ?: "$fiji_known_indels"
  refFasta = params.refFasta ?: "$fiji_ref_fasta"
}
// Check parameters for mendel
if (workflow.profile == 'mendel'){
  resultsDir = params.resultsDir ?: "$mendel_results_dir"
  dbsnp = params.dbsnp ?: "$mendel_dbsnp"
  ensemblGeneAnnotation = params.ensemblGeneAnnotation ?: "$mendel_ensembl_gene_annotation"
  knownIndels = params.knownIndels ?: "$mendel_known_indels"
  refFasta = params.refFasta ?: "$mendel_ref_fasta"
}
// Finally check if all the required parameters have either been provided via the commandline or
// by specifying a particular profile such as fiji or mendel
if (!resultsDir) exit 1, 'Specify the directory to store analysis results, see --help'  
if (!dbsnp) exit 1, 'Specify the dbsnp (*.vcf.gz) file, see --help'  
if (!ensemblGeneAnnotation) exit 1, 'Specify the ensemblGeneAnnotation (*.gff3.gz) file, see --help'  
if (!knownIndels) exit 1, 'Specify the known indels (*.vcf.gz) file, see --help'  
if (!refFasta) exit 1, 'Specify the refenrece fasta (*.fasta.gz) file, see --help'  

if (!params.runName) exit 1, 'Specify a run name, see --help'
runNameWithoutSpaces = params.runName.replaceAll('\\W','_')
outDir = "${resultsDir}/${runNameWithoutSpaces}"
log.info "Output dir: ${outDir}"

if (!workflow.commandLine.contains('-resume') && file(outDir).exists())
  exit 1, "${outDir} already exists, specify another results dir or run name"

/* End of parameter handling section */

// Handle input files
inputFiles = Channel.empty()
tsvPath = ''
if (params.sample) tsvPath = params.sample

if (tsvPath){
    tsvFile = file(tsvPath)
    if (!tsvFile.exists()) exit 1, "${tsvPath} does not exists!"
    inputFiles = LLabUtils.extractSample(tsvFile)
}else{
  // Define channel for reading file pairs
  Channel
      .fromFilePairs( params.reads, size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}"}
      .set {inputFiles}
}

if (params.verbose){
  println("Following samples will be processed")
  inputFiles.subscribe {println it}
}

// Generate a channel holding chromosome list
Channel.from(LLabUtils.getChrmList())
.set{chrmList}

// println params
include './lib/nf_modules/modules.nf' params(
                                outDir: outDir,
                                skipFastQc: params.skipFastQc,
                                skipQc: params.skipQc,
                                dbsnp: dbsnp,
                                ensemblGeneAnnotation: ensemblGeneAnnotation,
                                refFasta: refFasta,
                                knownIndels: knownIndels,
                                templateDir: templateDir,
                                platform: platform,
				                        customRunName:  params.runName,
                                threads: 32
                                )


workflow{
    RunFastQC(inputFiles)
    MapReads(inputFiles)
    CreateRecalibrationTable(MapReads.out[0])

    ApplyBQSR(
            CreateRecalibrationTable.out[0], 
            CreateRecalibrationTable.out[1]
            )
    
    CallVariants(
                    ApplyBQSR.out[0],
                    ApplyBQSR.out[1]
                    )
    
    RunGenomicsDBImport(
        CallVariants.out[0].collect(),
        CallVariants.out[1].collect(),
        CallVariants.out[2].collect(),
        chrmList
        )
    
    GenotypeGVCF(RunGenomicsDBImport.out)
    ConcatVCF(GenotypeGVCF.out.collect())
    RunCSQ(ConcatVCF.out[0])
    VariantEval(ConcatVCF.out[1], ConcatVCF.out[2], )
    
    RunMultiQC(
        RunFastQC.out[0].collect(),
        CreateRecalibrationTable.out[1].collect(),
        VariantEval.out
        )
} // end of workflow



workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

def initParamsToDefaults(){
  def null_path = ''
  params.reads = null_path
  params.runName = ''
  params.contactMail=''
  // params.resultsDir = null_path
  params.dbsnp = null_path
  params.ensemblGeneAnnotation = null_path
  params.refFasta = null_path
  params.knownIndels = null_path

  params.help = false  
  params.skipQc = false
  params.skipFastQc = false


}

/* Helper functions */

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}
def layerLabMessage() {
  log.info "Workflow For Somatic Variations ~ ${workflow.manifest.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def helpMessage() {
  // Display help message
  this.layerLabMessage()
  runStr = "nextflow run javaidm/layer_lab_dna_seq_vc "
  log.info "    Usage:"
  log.info "       $runStr --sample <file.tsv>"
  log.info "       $runStr --reads <shell glob pointing to the reads>"
  log.info "    --sample <file.csv>"
  log.info "       Specify a CSV file containing sample_id, path_to_read1, path_to_read2."
  log.info "    --reads <shell glob>"
  log.info "       Specify a shell glob pointing to the reads."
  log.info "    --results-dir <directory>"
  log.info "       Specify a directory to hold the analysis results."
  log.info "    --run-name"
  log.info "       Specify a run name, results will be stored under results-dir/run-name"
  log.info "    --manifest [optional <manifest.csv>]"
  log.info "       Specify an optional manifest file."  
  log.info "    --dbsnp [optional if already specified as part of a Nextflow Profile]"
  log.info "    --ensembl-gene-annotation [optional if already specified as part of a Nextflow Profile]"
  log.info "    --known-indels [optional if already specified as part of a Nextflow Profile]"
  log.info "    --ref-fasta [optional if already specified as part of a Nextflow Profile]"
  log.info "       Human Genome Reference File."  
  log.info "    --onlyQC"
  log.info "       Run only QC tools and gather reports."
  log.info "    --help"
  log.info "       you're reading it."
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow."
}