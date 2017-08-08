import scala.io.Source
import java.io._

import scala.collection.JavaConversions._

import java.util.Date;
import java.text.SimpleDateFormat;
import java.util.Calendar;

package rnasnped {
object VCF{

  final var emptyString = ".";
  final var emptyChar   = '.';
  final var pass        = "PASS";


  case class VC(contig: String = "",
                pos: Int = -1,
                ident: String = "",
                ref: Char = VCF.emptyChar,
		alt: Array[Char] = Array.empty,
                qual: Int = -1,
                filter: Array[String] = Array.empty,
                info: Map[String,String] = Map.empty);

  /////////////////////////////////////////////////////////////////////////////


  def readVCF(vcfFile: String) = {
    val stream = if(vcfFile == "-"){ io.Source.stdin } else {  io.Source.fromFile(vcfFile) }
    stream.getLines.filter(x => x.length > 0 && x(0) != '#').map( x => parseVCFLine(x))
  }

  /////////////////////////////////////////////////////////////////////////////

  def groupVCs(VCs: Iterator[VC]) = {
    new groupVCIterator(VCs.buffered)
  }

  /////////////////////////////////////////////////////////////////////////////

  class groupVCIterator(iter: BufferedIterator[VCF.VC]) extends Iterator[Array[VCF.VC]] {
    var prevVC = VCF.VC()
    def hasNext = iter.hasNext
    def next = {
      if (!iter.hasNext) {
        Iterator.empty.next
      } else {
        @scala.annotation.tailrec def untilNext(arr: Array[VCF.VC]): Array[VCF.VC] = {
          val vc = if(iter.hasNext) iter.head else VCF.VC()
          if (!this.elemCompare(vc, prevVC) || !iter.hasNext) {
            //print("%s: %d \n" format (prevVC.contig, arr.length))
            prevVC = vc
            arr
          } else{
            untilNext(iter.next +: arr)
          }
        }
        untilNext(Array.empty[VCF.VC])
      }
    }
    def elemCompare(vcA: VCF.VC, vcB: VCF.VC): Boolean = {
      (vcA.contig == vcB.contig && vcA.pos == vcB.pos)
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def parseVCFLine(line: String) = {

    val s = line.split('\t');

    val contig = s(0);
    val pos    = s(1).toInt;
    val ident  = s(2);
    val ref    = s(3)(0);
    val alt    = if(s(4).length == 0) Array(VCF.emptyChar) else s(4).split(',').map(x => x(0))
    val qual   = if (s(5) == ".") -1 else s(5).toFloat.toInt;
    val filter = if(s(6) contains ",") s(6).split(",") else s(6).split(";") // Because I fucked up earlier, split on comma and semi-colon
    val info   = if(s.length < 8) Map.empty[String,String] else s(7).split(";").map( x => if(x.contains('=') && x.last != '=') x.split('=') else Array(x, "") ).map{ case Array(f1,f2) => (f1,f2)}.toMap;

    VC(contig, pos, ident, ref, alt, qual, filter, info)
  }

  /////////////////////////////////////////////////////////////////////////////

  def write(VCs: TraversableOnce[VC], outFile: String) = {
    val outfd = if(outFile == "-") new BufferedWriter(new OutputStreamWriter(System.out, "utf-8")) else new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "utf-8"))
    val today = Calendar.getInstance.getTime()
    val datafmt = new SimpleDateFormat("yyyyMMdd")
    outfd.write("##fileformat=VCFv4.1\n");
    outfd.write("##fileDate=%s\n" format datafmt.format(today))
    outfd.write("##source=RNASNEP\n")
    outfd.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples\">\n");
    outfd.write("##INFO=<ID=SO,Number=.,Type=Integer,Description=\"SNP origin nodes\">\n");
    outfd.write("##INFO=<ID=SP,Number=.,Type=Integer,Description=\"Support for each origin node\">\n")
    outfd.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Number of reads observed for alternative base at each origin node\">\n");
    outfd.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
    outfd.write("##INFO=<ID=CA,Number=1,Type=Float,Description=\"Count of A\">\n")
    outfd.write("##INFO=<ID=CC,Number=1,Type=Float,Description=\"Count of C\">\n")
    outfd.write("##INFO=<ID=CG,Number=1,Type=Float,Description=\"Count of G\">\n")
    outfd.write("##INFO=<ID=CT,Number=1,Type=Float,Description=\"Count of T\">\n")
    outfd.write("##INFO=<ID=SAD,Number=.,Type=Integer,Description=\"Sample AD\">\n")
    outfd.write("##INFO=<ID=SCA,Number=.,Type=Integer,Description=\"Sample CA\">\n")
    outfd.write("##INFO=<ID=SCC,Number=.,Type=Integer,Description=\"Sample CC\">\n")
    outfd.write("##INFO=<ID=SCG,Number=.,Type=Integer,Description=\"Sample CG\">\n")
    outfd.write("##INFO=<ID=SCT,Number=.,Type=Integer,Description=\"Sample CT\">\n")
    outfd.write("##INFO=<ID=SDP,Number=.,Type=Integer,Description=\"Sample DP\">\n")
    outfd.write("##INFO=<ID=NN,Number=.,Type=String,Description=\"Origin node names\">\n")
    outfd.write("##INFO=<ID=NI,Number=.,Type=Integer,Description=\"Origin node IDS in tree\">\n")
    outfd.write("##INFO=<ID=NP,Number=1,Type=String,Description=\"Tree possibility description\">\n")
    outfd.write("##INFO=<ID=GN,Number=1,Type=String,Description=\"Gene Name\">\n")
    outfd.write("##INFO=<ID=GID,Number=1,Type=String,Description=\"Gene ID\">\n")
    outfd.write("##INFO=<ID=GCR,Number=1,Type=String,Description=\"Gene coding region, gene ID\">\n")
    outfd.write("##INFO=<ID=MTL,Number=1,Type=String,Description=\"Mating type locus ID\">\n")
    outfd.write("##INFO=<ID=ISO,Number=1,Type=String,Description=\"Isoform transcript ID\">\n")
    outfd.write("##INFO=<ID=GI,Number=1,Type=String,Description=\"Gene annotation\">\n")
    outfd.write("##INFO=<ID=DL,Number=1,Type=String,Description=\"Deleteriousness of SNP\">\n");
    outfd.write("##INFO=<ID=VP,Number=1,Type=String,Description=\"Allele of SNP\">\n");
    outfd.write("##INFO=<ID=DM,Number=1,Type=String,Description=\"Domain annotation\">\n");
    outfd.write("##INFO=<ID=CS,Number=1,Type=Float,Description=\"Support for the call (number of leaves with the snp / number of leaves)\">\n")
    outfd.write("##FILTER=<ID=NX,Description=\"No Expression\">\n")
    outfd.write("##FILTER=<ID=LX,Description=\"Low Expression\">\n")
    outfd.write("##FILTER=<ID=MSNP,Description=\"SNP multiple alternative alleles\">\n")
    outfd.write("##FILTER=<ID=NALT,Description=\"No SNP at this location\">\n")
    outfd.write("##FILTER=<ID=NO,Description=\"No suitable origin could be found\">\n")
    outfd.write("##FILTER=<ID=MO,Description=\"Multiple origins were found :S\">\n")

    outfd.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    VCs.foreach{
      x => outfd.write("%s\t%d\t%s\t%c\t%s\t%s\t%s\t%s\n" format (x.contig,
                                                                  x.pos,
                                                                  x.ident,
                                                                  x.ref,
                                                                  if(x.alt.length == 0) "." else x.alt.mkString(","),
                                                                  if(x.qual == -1) "." else x.qual.toString,
                                                                  x.filter.mkString(";"),
                                                                  x.info.map{ case (x,y) => "%s%s%s" format (x, if(y.length==0) "" else "=" ,y)}.mkString(";")));
    }
    outfd.close();
  }

  /////////////////////////////////////////////////////////////////////////////

  def toGFF(vc: VC) = {

    GFF.GFFEntry(vc.contig,
                 vc.ident,
                 "SNP",
                 vc.pos,
                 vc.pos+1,
                 emptyString,
                 emptyString,
                 emptyString,
                 "REF=%s;ALT=%s;QUAL=%d;FILT=%s;%s" format (vc.ref.toString,
                                                            vc.alt.mkString(","),
                                                            vc.qual,
                                                            vc.filter.mkString(","),
                                                            vc.info.map{ case (x,y) => "%s%s%s" format (x, if(y.length==0) "" else "=" ,y)}.mkString(";") ))

  }

}

}
