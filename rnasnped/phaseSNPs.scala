import scala.collection.JavaConversions._

import scala.io.Source
import java.io._

package rnasnped {
object phaseSNPs {

  def main(args: Array[String]) = {

    if (args.length < 3 || args(0) == "help") {
      help()
    } else {

      val vcfFile   = args(0);
      val outFile   = args(1);
      val fieldName = args(2);

      //print(vcfFile + "\n")

      val vcf  = VCF.readVCF(vcfFile);
      val iter = new groupVCIteratorGene(vcf.buffered, fieldName);

      val phasedSNPs = iter.filter(g => g.length > 0).map(g => phase(g, fieldName)).flatMap(x => x)
      //iter.foreach( g => print("--%d--\n" format g.length))
      //vcf.foreach(v => print(v))

      VCF.write(phasedSNPs, outFile)
    }

  }

  /////////////////////////////////////////////////////////////////////////////

  def help() = {
    println("phaseSNPs: Phase SNPs in VCF file within a gene based on VAFs")
    println("Usage: phaseSNPs <vcfFile> <outFile> <fieldName>")
    println("")
    println(" vcfFile:   The input VCF")
    println(" outFile:   The output VCF location")
    println(" fieldName: The name of the field to use as a gene identifier e.e. GCR")
  }

  /////////////////////////////////////////////////////////////////////////////

  class groupVCIteratorGene(iter: BufferedIterator[VCF.VC], fieldName: String) extends VCF.groupVCIterator(iter) {
    override def elemCompare(vcA: VCF.VC, vcB: VCF.VC) : Boolean = {
      val aField = if(vcA.info contains fieldName) vcA.info(fieldName) else "a";
      val bField = if(vcB.info contains fieldName) vcB.info(fieldName) else "b";
      //print (fieldName, aField, bField)
      //print(vcA)
      //print(vcB)
      //print("\n")
      (aField == bField)
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def phase(VCs: Array[VCF.VC], fieldName: String) = {
    val geneID  = VCs(0).info(fieldName)
    val VAFs    = VCs.map(x => x.info("AD").toFloat / x.info("DP").toInt)
    //if(geneID == "1343077"){
    //  print(geneID + "\n");
    //}
    (VCs zip VAFs).map{ case (vc, vaf) => VCF.VC(vc.contig, vc.pos, vc.ident, vc.ref, vc.alt, vc.qual, vc.filter, vc.info + ("VP" -> {if(vaf >= 0.5) "A" else "B"}))}
  }

  /////////////////////////////////////////////////////////////////////////////

}
}
