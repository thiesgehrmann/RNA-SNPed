
import scala.io.Source
import java.io._
import scala.collection.JavaConversions._

object cmpCalls {

  case class callCmp(contig: String,
                    pos: Int,
                    cmp: String,
                    refBase: Char,
                    altBase: List[List[Char]],
                    totalCount: List[Int],
                    countA: List[Int],
                    countC: List[Int],
                    countG: List[Int],
                    countT: List[Int],
                    countN: List[Int],
                    callFlags: List[String])


  /////////////////////////////////////////////////////////////////////////////

  def main(args: Array[String]) = {

    val fileA = args(0);
    val fileB = args(1);
    val out   = args(2);
    //val outA  = args(2);
    //val outB  = args(3);

    println("A: %s\nB: %s\nO: %s" format (fileA, fileB, out))

    val varAStream = SNP.readCalls(fileA);
    val varBStream = SNP.readCalls(fileB);

    println("Comparing");
    val cmps = (varAStream zip varBStream).map( x => cmp(x._1, x._2)).filter(x => x.cmp contains "P")

    writeCmps(cmps.toList, out);

    //val varPairs = (varAStream zip varBStream).map( x => (x._1, x._2, cmp(x._1, x._2))).map( x => (SNP.addFlags(x._1, prepareFlags(x._3._1, x._3._2, true)), SNP.addFlags(x._2, prepareFlags(x._3._1, x._3._2, false))))
    //val varAs    = varPairs.map( x => x._1);
    //val varBs    = varPairs.map( x => x._2);

    //SNP.writeCalls(List(varAs.toList).toList, outA);
    //SNP.writeCalls(List(varBs.toList).toList, outB);
  }

  /////////////////////////////////////////////////////////////////////////////

  def writeCmps(c : List[callCmp], out: String) = {
   
    println("Writing");
 
    val outfd = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(out), "utf-8"))
    
    outfd.write("#contig\tpos\tcall\trefbase\taltbase\ttotalcount\ta\tc\tg\tt\tn\tcallflags\n");
    
    c.foreach{
     x => outfd.write( "%s\t%d\t%s\t%s\t%s\t%d;%d\t%d;%d\t%d;%d\t%d;%d\t%d;%d\t%d;%d\t%s\n" format ( x.contig,
                                                                                   x.pos+1,
                                                                                   x.cmp,
                                                                                   x.refBase,
                                                                                   x.altBase.map( x => x.mkString(",")).mkString(";"),
                                                                                   x.totalCount(0), x.totalCount(1),
                                                                                   x.countA(0), x.countA(1),
                                                                                   x.countC(0), x.countC(1),
                                                                                   x.countG(0), x.countG(1),
                                                                                   x.countT(0), x.countT(1),
                                                                                   x.countN(0), x.countN(1),
                                                                                   x.callFlags.mkString(",")))
      }
      outfd.close();

  }


  /////////////////////////////////////////////////////////////////////////////

  def prepareFlags(cmptype: String, ea: Boolean, v: Boolean): List[String] = {
    // Invert the cmp types if it is NOT the first var (i.e. if it is var B)
    val cmptypeflag = if(v){cmptype} else{invTypes(cmptype)}
    val known       = if(isKnown(cmptype)){ cmpCalls.knownCmp } else{ cmpCalls.unknownCmp }
    if(ea){List(cmptypeflag, known)} else {List(cmptypeflag, known, cmpCalls.differentAlternates)}
  }

  /////////////////////////////////////////////////////////////////////////////

  val invTypes = Map( "PAPB" -> "PAPB",
                      "PAUB" -> "UAPB",
                      "PANB" -> "NAPB",
                      "UAPB" -> "PAUB",
                      "UAUB" -> "UAUB",
                      "UANB" -> "NAUB",
                      "NAPB" -> "PANB",
                      "NAUB" -> "UANB",
                      "NANB" -> "NANB");
  val isKnown = Map( "PAPB" ->  true,
                     "PAPB" ->  true,
                     "PAUB" ->  false,
                     "UAPB" ->  false,
                     "PANB" ->  true,
                     "NAPB" ->  true,
                     "UAPB" ->  false,
                     "PAUB" ->  false,
                     "UAUB" ->  false,
                     "UAUB" ->  false,
                     "UANB" ->  false,
                     "NAUB" ->  false,
                     "NAPB" ->  true,
                     "PANB" ->  true,
                     "NAUB" ->  false,
                     "UANB" ->  false,
                     "NANB" ->  true,
                     "NANB" ->  true);

  val presentA    = "PA";
  val unknownA    = "UA";
  val notPresentA = "NA";
  val presentB    = "PB";
  val unknownB    = "UB";
  val notPresentB = "NB";
  val knownCmp    = "KC";
  val unknownCmp  = "UC";
  val differentAlternates = "DA";

  /////////////////////////////////////////////////////////////////////////////

  def cmp(varInA: SNP.SNPcall, varInB: SNP.SNPcall): callCmp = {
    val ta = if(varInA.call){ cmpCalls.presentA  } else if( varInA.callFlags.contains(SNP.lowExpression)){ cmpCalls.unknownA } else { cmpCalls.notPresentA }
    val tb = if(varInB.call){ cmpCalls.presentB } else if( varInB.callFlags.contains(SNP.lowExpression)){ cmpCalls.unknownB } else { cmpCalls.notPresentB }

    val ea = eqAlt(varInA, varInB)
    new callCmp(varInA.contig,
                   varInA.pos,
                   ta+tb,
                   varInA.refBase,
                   if(ea){ List(varInA.altBase) } else { List( varInA.altBase, varInB.altBase ) },
                   List(varInA.totalCount, varInB.totalCount),
                   List(varInA.countA, varInB.countA),
                   List(varInA.countC, varInB.countC),
                   List(varInA.countG, varInB.countG),
                   List(varInA.countT, varInB.countT),
                   List(varInA.countN, varInB.countN),
                   if(ea) { varInA.callFlags.toSet.union(varInB.callFlags.toSet).toList } else {cmpCalls.differentAlternates :: varInA.callFlags.toSet.union(varInB.callFlags.toSet).toList })
  }
  /////////////////////////////////////////////////////////////////////////////

  //def cmp(varInA: SNP.SNPcall, varInB: SNP.SNPcall): (String, Boolean) = {
  //  val ta = if(varInA.call){ cmpCalls.presentA  } else if( varInA.callFlags.contains(SNP.lowExpression)){ cmpCalls.unknownA } else { cmpCalls.notPresentA }
  //  val tb = if(varInB.call){ cmpCalls.presentB } else if( varInB.callFlags.contains(SNP.lowExpression)){ cmpCalls.unknownB } else { cmpCalls.notPresentB }
  //  //println("###########################");
  //  //println(varInA)
  //  //println(varInB)
  //  //println(ta + tb);
  //  (ta + tb, eqAlt(varInA, varInB));
  //}

  /////////////////////////////////////////////////////////////////////////////

  def eqAlt(varInA: SNP.SNPcall, varInB: SNP.SNPcall): Boolean = {
    varInA.altBase == varInB.altBase
  }

}
