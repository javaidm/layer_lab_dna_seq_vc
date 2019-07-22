#!/usr/bin/env nextflow
nextflow.preview.dsl=2
initParamsToDefaults()
if (params.help) exit 0, helpMessage()
platform = "ILLUMINA"
templateDir = "${workflow.projectDir}/lib/bash_templates"
processParams()

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



// Generate a channel holding chromosome list
Channel.from(LLabUtils.getChrmList())
.set{chrmList}

// println params

include 'lib/nf_modules/modules.nf' params(
                                outDir: params.outDir,
                                skipFastQc: params.skipFastQc,
                                skipQc: params.skipQc,
                                dbsnp: params.dbsnp,
                                gff3: params.gff3,
                                refFasta: params.refFasta,
                                knownIndels: params.knownIndels,
                                templateDir: templateDir,
                                platform: platform,
				                        customRunName:  params.runName,
                                threads: 32
                                )


workflow{
    // Startup()
    
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
    RunCSQ(ConcatVCF.out)
    
    RunMultiQC(
        RunFastQC.out[0].collect(),
        CreateRecalibrationTable.out[1].collect()
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
  params.resultsDir = null_path
  params.dbsnp = null_path
  params.gff3 = null_path
  params.refFasta = null_path
  params.knownIndels = null_path

  params.help = false  
  params.skipQc = false
  params.skipFastQc = false


}

/* Helper functions */
def  processParams(){
  params.name = 'Layer Lab DNA Seq Analysis Pipeline'
  params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>
  params.contactMail ?: 'mahmood.javaid@colorado.edu'
  if (!params.sample && !params.reads){
      exit 1, 'You must sepcify either a sample tsv file or the reads, see --help'
  }

  if (!params.resultsDir) exit 1, 'Specify the directory to store analysis results, see --help'  
  if (!params.runName) exit 1, 'Specify a run name, see --help'
  runNameWithoutSpaces = params.runName.replaceAll('\\W','_')
  params.outDir = "${params.resultsDir}/${runNameWithoutSpaces}"
  log.info "Output directory for results is set to: ${params.outDir}"
  if (file(params.outDir).exists()) exit 1, "${params.outDir} already exists, specify another results dir or run name"
}

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
  log.info "    --onlyQC"
  log.info "       Run only QC tools and gather reports."
  log.info "    --help"
  log.info "       you're reading it."
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow."
}