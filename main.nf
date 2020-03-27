#!/usr/bin/env nextflow
nextflow.preview.dsl=2
initParamsToDefaults()
if (params.help) exit 0, helpMessage()
PLATFORM = "ILLUMINA"
THREADS = 32

templateDir = "${workflow.projectDir}/lib/bash_templates"
// processParams()
/* Process the parameters and set the environemnt */
params.name = 'Layer Lab DNA Seq Analysis Pipeline'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>
params.contactMail ?: 'mahmood.javaid@colorado.edu'
if (!params.sampleTsv && !params.reads){
    exit 1, 'You must sepcify either a sample tsv file or the reads, see --help'
}
// Check parameters for profile fiji
// If they are not provided at the command line, set them for fiji from the env section of the config
if (workflow.profile == 'fiji' || workflow.profile == 'fiji_hg37'){
    println('setting from hg_37...')
  resultsDir = params.resultsDir ?: "$fiji_results_dir"
  dbsnp = params.dbsnp ?: "$fiji_dbsnp"
  ensemblGeneAnnotation = params.ensemblGeneAnnotation ?: "$fiji_ensembl_gene_annotation"
  knownIndels = params.knownIndels ?: "$fiji_known_indels"
  refFasta = params.refFasta ?: "$fiji_ref_fasta"
}

// if (workflow.profile == 'fiji_hg38'){
//   resultsDir = params.resultsDir ?: "$fiji_results_dir"
//   dbsnp = params.dbsnp ?: "$fiji_dbsnp"
//   ensemblGeneAnnotation = params.ensemblGeneAnnotation ?: "$fiji_ensembl_gene_annotation"
//   knownIndels = params.knownIndels ?: "$fiji_known_indels"
//   refFasta = params.refFasta ?: "$fiji_refFasta"
// }
// // Check parameters for mendel
// if (workflow.profile == 'mendel'){
//   resultsDir = params.resultsDir ?: "$mendel_results_dir"
//   dbsnp = params.dbsnp ?: "$mendel_dbsnp"
//   ensemblGeneAnnotation = params.ensemblGeneAnnotation ?: "$mendel_ensembl_gene_annotation"
//   knownIndels = params.knownIndels ?: "$mendel_known_indels"
//   refFasta = params.refFasta ?: "$mendel_refFasta"
// }
// Finally check if all the required parameters have either been provided via the commandline or
// by specifying a particular profile such as fiji or mendel
if (!resultsDir) exit 1, 'Specify the directory to store analysis results, see --help'  
if (!dbsnp) exit 1, 'Specify the dbsnp (*.vcf.gz) file, see --help'  
if (!ensemblGeneAnnotation) exit 1, 'Specify the ensemblGeneAnnotation (*.gff3.gz) file, see --help'  
if (!knownIndels) exit 1, 'Specify the known indels (*.vcf.gz) file, see --help'  
if (!refFasta) exit 1, 'Specify the refenrece fasta (*.fasta.gz) file, see --help'  

if (!params.runName) exit 1, 'Specify a run name, see --help'
runNameWithoutSpaces = params.runName.replaceAll('\\W','_')
OUT_DIR = "${resultsDir}/${runNameWithoutSpaces}"
log.info "Output dir: ${OUT_DIR}"

if (!workflow.commandLine.contains('-resume') && file(outDir).exists())
  exit 1, "${outDir} already exists, specify another results dir or run name"

/* End of parameter handling section */
// Handle input files
inputFiles = Channel.empty()
tsvPath = ''
if (params.sampleTsv) tsvPath = params.sampleTsv

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
// ch_conditions = LLabUtils.extractConditions(file(params.conditionsTsv))
// ch_conditions = LLabUtils.extractConditions(file(params.conditionsTsv))

if (params.verbose){
  println("Following samples will be processed")
  inputFiles.subscribe {println it}
}

// Generate a channel holding chromosome list
Channel.from(LLabUtils.getChrmList())
.set{chrmList}

Channel.from(LLabUtils.getChrmListHg38())
.set{chrmListHg38}

workflow{
    // RunFastQC(inputFiles)
    MapReads(inputFiles)
    CreateRecalibrationTable(MapReads.out[0])

    ApplyBQSR(
            MapReads.out[0], 
            CreateRecalibrationTable.out
            )
    
    CallVariants(
                    ApplyBQSR.out[0],
                    ApplyBQSR.out[1]
                    )
    
    RunGenomicsDBImport(
        CallVariants.out[0].collect(),
        CallVariants.out[1].collect(),
        CallVariants.out[2].collect(),
        chrmListHg38
        )
    
    GenotypeGVCF(RunGenomicsDBImport.out)
    ConcatVCF(GenotypeGVCF.out.collect())
    RunCSQ(ConcatVCF.out[0])
    VariantEval(ConcatVCF.out[1], ConcatVCF.out[2], )
    
    // RunMultiQC(
    //     RunFastQC.out[0].collect(),
    //     CreateRecalibrationTable.out[1].collect(),
    //     VariantEval.out
    //     )
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
/* Processes */




// process ExtractFilesFromSampleDir{
//     input:
//     val(srcPath)
//     output
// }
process RunFastQC {
    tag "$sample"
    publishDir "${OUT_DIR}/fastqc", mode: 'copy', overwrite: false
    when:
    !params.skipQc && !params.skipFastQc

    input:
    set val(sample), file(reads)

    output:
    file "*.{zip,html}"
    file '.command.out'

    script:
    """
    fastqc -q $reads
    """
}



process MapReads {
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/align/" , mode: 'copy', overwrite: false
    input:
    set val(sample), file(reads)

    output:
    file(out_file)
    file("${out_file}.crai")

    script:
    r1 = reads[0]
    r2 = reads[1]
    
    script:
    out_file = "${sample}.cram"
   """
   bwa mem \
    -t $THREADS -R "@RG\tID:$sample\tSM:$sample\tPL:$PLATFORM\tPU:$sample\tLB:$sample" $refFasta $r1 $r2 \
    | samblaster \
    | samtools sort --output-fmt-option seqs_per_slice=4000 -O CRAM --reference $refFasta -m 18G -@ 6 /dev/stdin -o $out_file \
    && samtools index $out_file
   """
}

process CreateRecalibrationTable{
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/misc/recal_table/", mode: 'copy', overwrite: false

    input:
    file (sample_cram)
    // set val(sample), file(sample_cram_index_pair)
    
    output:
    // file(sample_cram)
    file("${sample}.recal.table")
    // file '.command.log'
    
    script:
    // sample_cram = sample_cram_index_pair[0]
    sample = sample_cram.baseName 
    """
    gatk BaseRecalibrator \
            -I  $sample_cram\
            --known-sites $dbsnp \
            --known-sites $knownIndels \
            -O ${sample}.recal.table \
            -R $refFasta
    """
}

process ApplyBQSR {
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/misc/BQSR/", mode: 'copy', overwrite: false

    input:
    file(sample_cram) 
    file(recalibrated_file)

    
    output:
    file("${sample}.bqsr.bam")
    file("${sample}.bqsr.bai")

    script:
    sample = "${sample_cram.baseName}"
    """
    gatk ApplyBQSR \
            -I $sample_cram \
            -R $refFasta \
            -bqsr $recalibrated_file \
            -O ${sample}.bqsr.bam 
    """
}

process CallVariants {
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/variants/", mode: 'copy', overwrite: false

    input:
    file(file_bqsr_bam)
    file('*')
    
    output:
    file "${sample}.gvcf.gz" 
    file "${sample}.gvcf.idx" 
    file "${sample}.gvcf.gz.tbi"
    
    script:
    sample = file_bqsr_bam.getSimpleName()
    """
    gatk HaplotypeCaller \
        -I $file_bqsr_bam \
        --dbsnp $dbsnp \
        -O ${sample}.gvcf \
        --emit-ref-confidence GVCF \
        -R $refFasta
    bgzip ${sample}.gvcf
    tabix ${sample}.gvcf.gz
    """
}



process RunGenomicsDBImport {
    echo true
    publishDir "${OUT_DIR}/misc/genomicsdb/", mode: 'copy', overwrite: false

    input:
    file('*')
    file('*')
    file('*')
    val (chr)
    
    output:
    file (gDB)

    script:
    sample_map="cohort_samples.map"
    gDB = chr
    """
    for x in `ls *.gvcf.gz`
    do
        sample=`echo \$x | cut -d '.' -f1`
        echo "\${sample}\t\${x}" >> ${sample_map}
    done
    gatk GenomicsDBImport \
    --genomicsdb-workspace-path ${gDB} \
    -L ${chr} \
    --sample-name-map ${sample_map} \
    --reader-threads 8
    """
}


process GenotypeGVCF{
    echo true
    publishDir "${OUT_DIR}/misc/genotype_gvcf/", mode: 'copy', overwrite: false

    input:
    file(genomics_db_workspace)

    output:
    file "$out_file"

    script:
    out_file = "${genomics_db_workspace.baseName}.vcf"  
    """
    gatk GenotypeGVCFs \
    -R $refFasta \
    -O $out_file  \
    -V gendb://${genomics_db_workspace}
    """  
}

process ConcatVCF{
    echo true
    publishDir "${OUT_DIR}/results/", mode: 'copy', overwrite: false

    input:
    file('*')

    output:
    file(out_file)
    file "${out_file}.gz" 
    file "${out_file}.gz.tbi"

    script:
    out_file='cohort_joint.vcf'
    """
    bcftools concat -o ${out_file} -Ov \
    chr1.vcf \
    chr2.vcf \
    chr3.vcf \
    chr4.vcf \
    chr5.vcf \
    chr6.vcf \
    chr7.vcf \
    chr8.vcf \
    chr9.vcf \
    chr10.vcf \
    chr11.vcf \
    chr12.vcf \
    chr13.vcf \
    chr14.vcf \
    chr15.vcf \
    chr16.vcf \
    chr17.vcf \
    chr18.vcf \
    chr19.vcf \
    chr20.vcf \
    chr21.vcf \
    chr22.vcf \
    chrX.vcf \
    chrY.vcf

    bgzip -c ${out_file} > ${out_file}.gz
    tabix -p vcf ${out_file}.gz
"""
}

process RunCSQ{
    echo true
    publishDir "${OUT_DIR}/csq/", mode: 'copy', overwrite: false

    input:
    file(vcf)

    output:
    file "$out_file"

    script:
    out_file = "cohort.bcf"    
    """
    bcftools norm -m- -Ov  -f $refFasta -w 10000 $vcf |\
    bcftools csq -s - --ncsq 40 -g $ensemblGeneAnnotation -l -f $refFasta -Ob -o $out_file
    """
}

process VariantEval{
    echo true
    publishDir "${OUT_DIR}/variant_eval/", mode: 'copy', overwrite: false

    input:
    file(cohort_vcf)
    file(vcf_index)

    output:
    file "$out_file"

    script:
    out_file = "cohort.eval.grp"    
    script:
    """ 
    gatk VariantEval --eval $cohort_vcf --comp $dbsnp -R $refFasta --output $out_file
    """  
}

process RunMultiQC {
    publishDir "${OUT_DIR}/multiqc", mode: 'copy', overwrite: false

    input:
    file (fastqc:'fastqc/*')
    file ('gatk_base_recalibration/*')
    file ('gatk_variant_eval/*')
    
    output:
    file '*multiqc_report.html'
    file '*_data'
    file '.command.err'
    val prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename  . 2>&1
    
    """
}


