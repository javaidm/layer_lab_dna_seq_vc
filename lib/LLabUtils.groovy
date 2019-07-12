import static nextflow.Nextflow.file
import nextflow.Channel

class LLabUtils {

  // Check if a row has the expected number of item
  static def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  // Check parameter existence
  static def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
  }
  // Compare params to list of verified params
  static def isAllowedParams(params) {
    def test = true
    params.each{
      if (!checkParams(it.toString().split('=')[0])) {
        println "params ${it.toString().split('=')[0]} is unknown"
        test = false
      }
    }
    return test
  }

  // Compare each parameter with a list of parameters
  // static def checkParameterList(list, realList) {
  //   return list.every{ checkParameterExistence(it, realList) }
  // }

  // Return element in list of allowed params
  static def checkParams(it) {
    return it in [
      'contact-mail',
      'contactMail',
      'container',
      'dbsnp',
      'gff3'
      'help',
      'known-indels',
      'knownIndels',
      'out-dir',
      'outDir',
      'publish-dir-mode',
      'publishDirMode',
      'ref-fasta',
      'refFasta',
      'sample-dir',
      'sampleDir',
      'sample-tsv',
      'sampleTsv',
      'tag',
      'verbose',
      'version']
  }

  // Check file extension
  static def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
  }


  // Return file if it exists
  static def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
  }

  // Return status [0,1]
  // 0 == Normal, 1 == Tumor
  static def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
  }

  // Layer Lab ascii art
  static def layerLabAscii() {

  println "| |                         | |         | |    "
  println "| |     __ _ _   _  ___ _ __| |     __ _| |__  "
  println "| |    / _` | | | |/ _ \ '__| |    / _` | '_ \ "
  println "| |___| (_| | |_| |  __/ |  | |___| (_| | |_) |"
  println "|______\__,_|\__, |\___|_|  |______\__,_|_.__/ "
  println "              __/ |                            "
  println "             |___/                             "
  }

}
