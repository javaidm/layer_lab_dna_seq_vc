#!/usr/bin/env nextflow
nextflow.preview.dsl=2

platform = "ILLUMINA"
// templateDir = "${PWD}/lib/bash_templates"
templateDir = "${workflow.projectDir}/lib/bash_templates"
params.name = 'Layer Lab DNA Seq Analysis Pipeline'

params.help = false // Don't give help information
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>
params.verbose = false // Enable for more verbose information
params.contactMail='mahmood.javaid@colorado.edu'
params.reads = "/scratch/Shares/layer/nextflow/kristen/three_sample_small_fastq_files/**/*_R{1,2}_*.fastq.gz"

params.customRunName = params.name
params.skipFastqc = false
params.skipQc = false
params.outDir="analysis"

// Include modules from the module.nf
include 'lib/nf_modules/modules.nf' params(
                                outDir: params.outDir,
                                skipFastqc: params.skipFastqc,
                                skipQc: params.skipQc,
                                dbsnp: params.dbsnp,
                                gff3: params.gff3,
                                refFasta: params.refFasta,
                                knownIndels: params.knownIndels,
                                templateDir: templateDir,
                                platform: platform,
				customRunName:  params.customRunName,
                                threads: 32
                                )

// Define channel for reading file pairs
Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}"}
    .set {read_files}

// Generate a channel holding chromosome list
Channel.from(get_chrm_list())
.set{chrm_list}

workflow{
    RunFastQC(read_files)
    MapReads(read_files)
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
        chrm_list
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

