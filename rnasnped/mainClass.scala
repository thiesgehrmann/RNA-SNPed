package rnasnped{

object mainClass{

  val operations = Map( "assignOrigin"  -> assignOrigin.main _,
                        "filterVCF"     -> filterVCF.main _,
                        "variantCaller" -> variantCaller.main _,
                        "slidingWindowSNPs2" -> slidingWindowSNPs2.main _,
                        "phaseSNPs"          -> phaseSNPs.main _,
                        "help"          -> help _)

  /////////////////////////////////////////////////////////////////////////////

  def main(args: Array[String]) = {

   if (args.length < 1 || !operations.contains(args(0)) ) help(Array.empty[String]) else operations(args(0))(args.drop(1))

  }

  /////////////////////////////////////////////////////////////////////////////

  def help(args: Array[String]) = {
    println("rnaSNPed toolkit")
    println("Usage: java -jar rnasnped.jar <task> <options>")
    println("")
    println(" <task> is one of the following:")
    println("   variantCaller      - to call variants from BAM files")
    println("   assignOrigin       - to call the origin of the SNPs in a sample tree")
    println("   filterVCF          - to perform filtering operations on the resulting VCF files")
    println("   slidingWindowSNPs2 - to count SNPs in a sliding window")
    println("   phaseSNPs          - to phase SNPs within a gene based on VAFs")
    println("   help               - to see this help message")
    println("")
    println("P.S. you can see the help for each of those options by having the second argument be 'help', e.g. 'assignOrigin help'")
  }

  /////////////////////////////////////////////////////////////////////////////

}

}
