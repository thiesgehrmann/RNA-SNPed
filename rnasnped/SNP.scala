import scala.io.Source
import java.io._
import java.util.Date;
import java.text.SimpleDateFormat;
import java.util.Calendar;

package rnasnped {
object SNP {

  val noExpression       = "NX";
  val lowExpression      = "LX";
  val possibleDeletion   = "PD";
  val possibleInsertion  = "PI";
  val nearSpliceJunction = "SJ";
  val referenceIsN       = "RN";
  val referenceMax       = "RM";
  val alternateMax       = "AM";

  

  case class SNPcall(contig: String,
                     pos: Int,
                     call: Boolean,
                     refBase: Char,
                     altBase: List[Char],
                     totalCount: Int,
                     countA: Int,
                     countC: Int,
                     countG: Int,
                     countT: Int,
                     countN: Int,
                     callFlags: List[String],
                     qual: List[Int])

  /////////////////////////////////////////////////////////////////////////////

  def addFlags(snp: SNPcall, flags: List[String]): SNPcall = {
    new SNPcall(snp.contig, snp.pos, snp.call, snp.refBase, snp.altBase, snp.totalCount, snp.countA, snp.countC, snp.countG, snp.countT, snp.countN, snp.callFlags ++ flags, snp.qual)
  }

  /////////////////////////////////////////////////////////////////////////////

  def writeVCF(calls: TraversableOnce[SNP.SNPcall], outfile: String, sample_id: String) = {
    val outfd = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfile), "utf-8"))

    println("Writing output to %s\n" format outfile);

    outfd.write("##fileformat=VCFv4.2\n");
    val today = Calendar.getInstance.getTime()
    val datafmt = new SimpleDateFormat("yyyyMMdd")
    outfd.write("##fileDate=%s\n" format datafmt.format(today))
    outfd.write("##source=RNASNEP\n")
    outfd.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
    outfd.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
    outfd.write("##INFO=<ID=CA,Number=1,Type=Float,Description=\"Count of A\">\n")
    outfd.write("##INFO=<ID=CC,Number=1,Type=Float,Description=\"Count of C\">\n")
    outfd.write("##INFO=<ID=CG,Number=1,Type=Float,Description=\"Count of G\">\n")
    outfd.write("##INFO=<ID=CT,Number=1,Type=Float,Description=\"Count of T\">\n")
    outfd.write("##INFO=<ID=CN,Number=1,Type=Float,Description=\"Count of N\">\n")
    outfd.write("##INFO=<ID=QA,Number=A,Type=Float,Description=\"Quality of Alternates\">\n")
    outfd.write("##INFO=<ID=LX,Number=0,Type=Flag,Description=\"Alternative allele has low expression\">\n")
    outfd.write("##FILTER=<ID=NX,Description=\"No Expression\">\n")
    outfd.write("##FILTER=<ID=PI,Description=\"Near/In Possible Insertion\">\n")
    outfd.write("##FILTER=<ID=PD,Description=\"Near/In Possible Deletion\">\n")
    outfd.write("##FILTER=<ID=SJ,Description=\"Near Splice Junction\">\n")
    outfd.write("##FILTER=<ID=RN,Description=\"Reference is N\">\n")
    outfd.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    calls.foreach{
        x => outfd.write("%s\t%d\t%s\t%c\t%s\t%s\t%s\tDP=%d;CA=%d;CC=%d;CG=%d;CT=%d;CN=%d%s%s\n" format ( x.contig,
                                                                                      x.pos+1,
                                                                                      "%s;%d;%s" format ( x.contig, x.pos+1, sample_id),
                                                                                      x.refBase,
                                                                                      if (x.call) x.altBase.mkString(",") else ".",
                                                                                      if (x.call) "%d" format x.qual.reduceLeft(_ min _) else ".",
                                                                                      if (x.call) "PASS" else if(x.callFlags.length == 0) "." else x.callFlags.mkString(";"),
                                                                                      x.totalCount,
                                                                                      x.countA,
                                                                                      x.countC,
                                                                                      x.countG,
                                                                                      x.countT,
                                                                                      x.countN,
                                                                                      if (x.call) ";QA=" + x.qual.map( y => "%d" format y).mkString(",") else "",
                                                                                      if (x.callFlags.contains(SNP.lowExpression)){ ";LX" } else { "" }))
      }
      outfd.close();

  }

  /////////////////////////////////////////////////////////////////////////////

  def writeCalls(calls: List[SNP.SNPcall], outfile: String) = {
    val outfd = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfile), "utf-8"))
    
    outfd.write("#contig\tpos\tcall\trefbase\taltbase\ttotalcount\ta\tc\tg\tt\tn\tcallflags\n");
    
    calls.foreach{
        x => outfd.write( "%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" format ( x.contig,
                                                                                      x.pos+1,
                                                                                      if(x.call){ "T" }else{ "F"},
                                                                                      x.refBase,
                                                                                      x.altBase.mkString(","),
                                                                                      x.totalCount,
                                                                                      x.countA,
                                                                                      x.countC,
                                                                                      x.countG,
                                                                                      x.countT,
                                                                                      x.countN,
                                                                                      x.callFlags.mkString(",")))
      }
      outfd.close();
  }

  /////////////////////////////////////////////////////////////////////////////

  def readCalls(fileName: String): Iterator[SNP.SNPcall] = {
    Source.fromFile(fileName).getLines().filter(x => x.head != '#' || x.length == 0).map( line => readLine(line) )
  }

  /////////////////////////////////////////////////////////////////////////////

  def readLine(line: String): SNP.SNPcall = {
    val lineSplit = line.split('\t');
    if(lineSplit.length != 12){
      println(line)
    }
    new SNP.SNPcall(lineSplit(0),
                    lineSplit(1).toInt,
                    if(lineSplit(2) == "T") { true } else { false },
                    lineSplit(3).head,
                    lineSplit(4).split(',').map(x => x.head).toList,
                    lineSplit(5).toInt,
                    lineSplit(6).toInt,
                    lineSplit(7).toInt,
                    lineSplit(8).toInt,
                    lineSplit(9).toInt,
                    lineSplit(10).toInt,
                    lineSplit(11).split(',').toList,
                    lineSplit(12).split(',').toList.map( x => x.toInt))
  }

  /////////////////////////////////////////////////////////////////////////////

  def toGFF(call: SNP.SNPcall): GFF.GFFEntry = {
    new GFF.GFFEntry(call.contig,
                     "RNAseq",
                     "SNP",
                     call.pos,
                     call.pos,
                     "1000",
                     ".",
                     ".",
                     "baseCounts=%d,%d,%d,%d,%d;flags=%s" format (call.countA, call.countC, call.countG, call.countT, call.countN,call.callFlags.mkString(",")))
  }

  /////////////////////////////////////////////////////////////////////////////

}

}
