import static nextflow.Nextflow.file
import nextflow.Channel

class LLabUtils {

static def extractSample(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: "\t")
  .map { row ->
    def idPatient  = row[0]
    def file1      = this.returnFile(row[1])
    def file2      = this.returnFile(row[2])
    [idPatient,[ file1, file2]]
  }
} // end of method extractSample()

static def runSanityChecks(dirPath){
  // this.checkIfBothReadsExist(dirPath)
}




// Layer Lab ascii art
  static def layerLabAscii() {
    // def ascii_str = 
    println \
    '''
    | |                         | |         | |    
    | |     __ _ _   _  ___ _ __| |     __ _| |__ 
    | |    / _` | | | |/ _ \\ '__| |    / _` | '_ \\ 
    | |___| (_| | |_| |  __/ |  | |___| (_| | |_) |
    |______\\__,_|\\__, |\\___|_|  |______\\__,_|_.__/ 
                __/ |                            
               |___/  
    '''
    // println ascii_str
  }

  // Check if a row has the expected number of item
  static def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  static def getChrmList(){
    def chrs = (1..22).collect()
    chrs.addAll(['X', 'Y', 'MT'])
    return chrs
  }
  
  // // Return element in list of allowed params
  static def checkParams(it) {
    return it in [
      'contact-mail',
      'contactMail',
      'container',
      'dbsnp',
      'gff3',
      'help',
      'known-indels',
      'knownIndels',
      'results-dir',
      'resultsDir',
      'publish-dir-mode',
      'publishDirMode',
      'ref-fasta',
      'refFasta',
      'run-name',
      'runName',
      'sample-dir',
      'sampleDir',
      'sample-tsv',
      'sampleTsv',
      'tag',
      'verbose',
      'version'
      ]
  }

  // // Check file extension
  static def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
  }

  static def exists (path){
    return file(path).exists()
  }
  // // Return file if it exists
  static def returnFile(it) {
    if (!this.exists(it)){
      println("Missing file listed  in TSV file: ${it}, see --help for more information")
      System.exit(1)
    }
    return file(it)
  }

  // // Return status [0,1]
  // // 0 == Normal, 1 == Tumor
  static def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
  }
}
