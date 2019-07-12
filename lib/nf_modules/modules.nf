
// set shorter versions of the  the passed parameter names
threads = params.threads
platform = params.platform 
outDir = params.outDir

dbsnp = params.dbsnp
known_indels = params.knownIndels
ref_fasta = params.refFasta
gff3 = params.gff3
templateDir=params.templateDir
customRunName = params.customRunName

process RunFastQC {
    tag "$sample"
    publishDir "${outDir}/fastqc", mode: 'copy'
    when:
    !params.skipQc && !params.skipFastqc

    input:
    set val(sample), file(reads)

    output:
    file "*.{zip,html}"
    file '.command.out'

    script:
    template "${templateDir}/run_fastqc.sh"
}



process MapReads {
    echo true
    tag "$sample"
    publishDir "${outDir}/align/"
    input:
    set val(sample), file(reads)

    output:
    file("${sample}.cram")
    file '.command.log'

    script:
    r1 = reads[0]
    r2 = reads[1]
    
    script:
    template "${templateDir}/map_reads.sh"
}

process CreateRecalibrationTable{
    echo true
    tag "$sample"
    publishDir "${outDir}/misc/recal_table/", mode: 'copy'

    input:
    file (sample_cram)
    
    output:
    file(sample_cram)
    file("${sample}.recal.table")
    file '.command.log'
    
    script:
    sample = "${sample_cram.baseName}" 
    template "${templateDir}/create_recalibration_table.sh"
}

process ApplyBQSR {
    echo true
    tag "$sample"
    publishDir "${outDir}/misc/BQSR/", mode: 'copy'

    input:
    file(sample_cram) 
    file(recalibrated_file)

    
    output:
    file("${sample}.bqsr.bam")
    file("${sample}.bqsr.bai")

    script:
    sample = "${sample_cram.baseName}"
    template "${templateDir}/apply_bqsr.sh"
}

process CallVariants {
    echo true
    tag "$sample"
    publishDir "${outDir}/variants/", mode: 'copy'

    input:
    file(file_bqsr_bam)
    file('*')
    
    output:
    file "${sample}.gvcf.gz" 
    file "${sample}.gvcf.idx" 
    file "${sample}.gvcf.gz.tbi"
    
    script:
    sample = file_bqsr_bam.getSimpleName()
    template "${templateDir}/call_variants.sh"
}



process RunGenomicsDBImport {
    echo true
    publishDir "${outDir}/misc/genomicsdb/"

    input:
    file('*')
    file('*')
    file('*')
    val (chr)
    
    output:
    file (gDB)

    script:
    sample_map="cohort_samples.map"
    gDB = "chr${chr}"
    template "${templateDir}/run_genomics_db_import.sh" 
}


process GenotypeGVCF{
    echo true
    publishDir "${outDir}/misc/genotype_gvcf/"

    input:
    file(genomics_db_workspace)

    output:
    file "$out_file"

    script:
    out_file = "${genomics_db_workspace.baseName}.vcf"  
    template "${templateDir}/genotype_gvcf.sh"  
}

process ConcatVCF{
    echo true
    publishDir "${outDir}/results/"

    input:
    file('*')

    output:
    file(out_file)

    script:
    out_file='cohort_joint.vcf'
    template "${templateDir}/concat_vcf.sh" 
}

process RunCSQ{
    echo true
    publishDir "${outDir}/csq/"

    input:
    file(vcf)

    output:
    file "$out_file"

    script:
    out_file = "cohort.bcf"    
    template "${templateDir}/run_csq.sh"  
}

process RunMultiQC {
    publishDir "${outDir}/multiqc", mode: 'copy'

    input:
    file (fastqc:'fastqc/*')
    file ('gatk_base_recalibration/*')
    
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


// Some Helper methods

def get_chrm_list(){
    chrs = (1..22).collect()
    chrs.addAll(['X', 'Y', 'MT'])
    return chrs
}
