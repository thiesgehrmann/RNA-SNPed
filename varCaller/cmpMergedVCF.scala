import scala.collection.JavaConversions._

import scala.io.Source
import java.io._

import scala.collection.mutable.ArrayBuffer

//import VCF;

object assignOrigin{

  final val minExpr        =  3
  final val possible_yes   =  1
  final val possible_maybe =  0
  final val possible_no    = -1
  final val infoFeats      = Array("CA", "CC", "CG", "CT", "DP")

  /////////////////////////////////////////////////////////////////////////////

  def main(args: Array[String]) = {

    val inputFile  = args(0);
    val treeFile   = args(1);
    val outputFile = args(2);
    val gatkFlag   = if(args.length > 3) true else false

    val VCGroups   = new groupVCIterator(VCF.readVCF(inputFile));
    val sTree      = sampleTree.readTree(treeFile)

    val cmpVCs = if(!gatkFlag) {
      VCGroups.filter(x => x.length > 0).map( x => assignOriginToVCGroup(x, sTree))
    } else {
      VCGroups.filter(x => x.length > 0).map( x => assignOriginToVCGroupGATK(x, sTree))
    }
    VCF.write(cmpVCs, outputFile);

  }

  /////////////////////////////////////////////////////////////////////////////

  class groupVCIterator(iter: Iterator[VCF.VC]) extends Iterator[Array[VCF.VC]] {
    var prevContig = "";
    var prevPos    = 0;
    def hasNext = iter.hasNext
    def next = {
      if (!iter.hasNext)
        Iterator.empty.next
      @scala.annotation.tailrec def untilNext(arr: Array[VCF.VC]): Array[VCF.VC] = {
        val vc = iter.next()
        if (vc.contig != prevContig || vc.pos != prevPos || !iter.hasNext) {
          if(iter.hasNext){
            val nextVC = (iter take 1).next
            prevContig = nextVC.contig
            prevPos    = nextVC.pos
          }
          arr
        } else {
          untilNext(vc +: arr)
        }
      }
      untilNext(Array.empty[VCF.VC])
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def assignOriginToVCGroup(VCs: Array[VCF.VC], sTree: Array[sampleTree.treeNode]) = {
    val firstVC       = VCs(0)
    val callFlags     = VCs.map( x => x.filter).flatten.distinct
    val legalBase     = !(callFlags.contains(SNP.possibleDeletion) || callFlags.contains(SNP.possibleInsertion) || callFlags.contains(SNP.nearSpliceJunction) || firstVC.ref == 'N')
    val observedBases = VCs.map( x => x.alt).flatten.distinct
    val observedAlts  = observedBases.filter( b => b != firstVC.ref && b != VCF.emptyChar )

    //print(VCs)

    if (!legalBase){
      VCF.VC(firstVC.contig,
             firstVC.pos,
             VCF.emptyString,
             firstVC.ref,
             observedAlts,
             -1,
             callFlags.filter( x => x != VCF.emptyString && x != VCF.pass),
             Map.empty[String,String])
    } else if(observedAlts.length == 0){
      VCF.VC(firstVC.contig,
             firstVC.pos,
             VCF.emptyString,
             firstVC.ref,
             Array(VCF.emptyChar),
             -1,
             Array("NOSNP"),
             Map.empty[String,String])
    } else if(observedAlts.length > 1){
      VCF.VC(firstVC.contig,
             firstVC.pos,
             VCF.emptyString,
             firstVC.ref,
             observedAlts,
             -1,
             Array("MSNP"),
             Map.empty[String,String])
    } else {
   
      val alt = observedAlts(0);

      // Get the observed counts for the alternative base.
      // If it is 0, but there is depth >3 at the node, then set the count to -1
      val infoFeatCounts = infoFeats.map( f => (f, VCs.map( y => (y.ident.split('.')(1).toInt - 1, y.info.getOrElse(f, "0").toInt) ).toMap) ).toMap
      val observedAltExpr = VCs.map{
        y =>
         val exprA = y.info.getOrElse("C%c" format alt, "0").toInt
         (y.ident.split('.')(1).toInt - 1, if(exprA > 0 || y.info.getOrElse("DP", "0").toInt < minExpr) exprA else -1 )
      }
      val observedAltExprMap = observedAltExpr.toMap

      val (nodeExpr, nodeSupport, nodePossible, nodeInfoFeatCounts) = calcNodeFeatures(sTree, observedAltExprMap, infoFeatCounts);

      if(nodeSupport.last == 1) {
        // This is a somatic mutation on a leaf node ONLY if we have at least three reads supporting it!
        val ((id, expr), index) = observedAltExpr.zipWithIndex.filter(x => x._1._2 > 0).last
        VCF.VC(firstVC.contig,
               firstVC.pos,
               "%s;%s;%d" format (firstVC.contig, firstVC.pos, id),
               firstVC.ref,
               VCs(index).alt,
               (100.0 / VCs.length).toInt,
               if(expr < minExpr) Array("LX") else Array(VCF.pass),
               Map("SO" -> "1",
                   "SP" -> "1",
                   "NS" -> VCs.length.toString,
                   "AD" -> expr.toString,
                   "DP" -> nodeInfoFeatCounts(id)("DP").toString,
                   "CA" -> nodeInfoFeatCounts(id)("CA").toString,
                   "CC" -> nodeInfoFeatCounts(id)("CC").toString,
                   "CG" -> nodeInfoFeatCounts(id)("CG").toString,
                   "CT" -> nodeInfoFeatCounts(id)("CT").toString,
                   "NN" -> sTree(id).name,
                   "NI" -> id.toString))

      } else{

        // Determine which node is best!, and filter those origins
        val origins = determineOrigin(sTree, nodePossible).filter( o => !(nodeExpr(o) < minExpr) || nodeSupport(o) > 1)

        // This is a debug message
        //print("###################\n%s:%d: ALT=%c VCs=%d\n%s\n--\n" format (firstVC.contig, firstVC.pos, alt, VCs.length, formatNodeFeatures(sTree, nodeExpr, nodeSupport, nodePossible, origins)))

        // If there was no supported origin, then there is no SNP here!
        if(origins.length == 0) {
          VCF.VC(firstVC.contig,
                 firstVC.pos,
                 ".",
                 firstVC.ref,
                 Array.empty[Char],
                 -1,
                 Array("NO"),
                 Map.empty[String,String])
        }
        else{ //Otherwise
          VCF.VC(firstVC.contig,
                 firstVC.pos,
                 "%s;%d;%s.%s" format (firstVC.contig, firstVC.pos, origins.map(x => sTree(x).name).mkString(","), origins.map( x => (x + 1).toString).mkString(",")),
                 firstVC.ref,
                 if(observedBases contains firstVC.ref) Array(firstVC.ref, alt) else Array(alt),
                 (100.0*(origins.map(o => nodeSupport(o)).foldLeft(0)((a,b) => math.max(a,b) ).toFloat / VCs.length)).toInt,
                 if(origins.length == 1) Array("PASS") else Array("MO"),
                 Map("SO" -> origins.length.toString,
                     "SP" -> origins.map( o => nodeSupport(o).toString).mkString(","),
                     "NS" -> VCs.length.toString,
                     "AD" -> origins.map( o => nodeExpr(o)).foldLeft(0)(_ + _).toString,
                     "DP" -> origins.map( id => nodeInfoFeatCounts(id)("DP")).foldLeft(0)(_ + _).toString,
                     "CA" -> origins.map( id => nodeInfoFeatCounts(id)("CA")).foldLeft(0)(_ + _).toString,
                     "CC" -> origins.map( id => nodeInfoFeatCounts(id)("CC")).foldLeft(0)(_ + _).toString,
                     "CG" -> origins.map( id => nodeInfoFeatCounts(id)("CG")).foldLeft(0)(_ + _).toString,
                     "CT" -> origins.map( id => nodeInfoFeatCounts(id)("CT")).foldLeft(0)(_ + _).toString,
                     "NN" -> origins.map( o => sTree(o).name).mkString(","),
                     "NI" -> origins.map( _.toString).mkString(","),
                     "NP" -> formatNodePossible(sTree, nodePossible)))
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def calcNodeFeatures(sTree: Array[sampleTree.treeNode], observedAltExprMap: Map[Int,Int], infoFeatCounts: Map[String,Map[Int,Int]]) = {
    // Construct two arrays describing 
    // expr: The total number of reads observed beneath this node
    // support: the total number of samples with reads observed beneath this node
    // possible: Whether or not the alternative allele is POSSIBLE at this node
    //  Leaf nodes with NX/DP=0 (which have been filtered out) are considered as being a possibility
    //  However, internal nodes (nodes with children) may not have conflicting nodes
    var expr     = collection.mutable.ArrayBuffer.fill(sTree.length)(0)
    var support  = collection.mutable.ArrayBuffer.fill(sTree.length)(0)
    var possible = collection.mutable.ArrayBuffer.fill(sTree.length)(possible_maybe)
    var infos    = collection.mutable.ArrayBuffer.fill(sTree.length)(collection.mutable.Map.empty[String,Int])

    sTree.foreach{
      case sampleTree.treeNode(id, name, children, isLeaf, height) =>
        expr(id)     = if(isLeaf) observedAltExprMap.getOrElse(id, 0) else children.map(c => math.max(0, expr(c))).foldLeft(0)(_ + _);
        support(id)  = if(isLeaf){ if(expr(id) > 0) 1 else 0 } else children.map(c => support(c)).foldLeft(0)(_ + _);

        possible(id) = if(isLeaf) { if(expr(id) > 0) possible_yes else if(expr(id) == 0) possible_maybe else possible_no }
                       else{ children.map(c => possible(c)).foldLeft(possible_yes)( (a,b) => math.min(a,b))  }

        infoFeats.foreach{ f =>
          if(isLeaf)  infos(id)(f) = infoFeatCounts.getOrElse(f, Map.empty[Int,Int]).getOrElse(id, 0)
          else        infos(id)(f) = children.map(c => infos(c)(f)).foldLeft(0)(_ + _)
        }
    }

    (expr.toArray, support.toArray, possible.toArray, infos.map(x => x.toMap).toArray)
  }

  /////////////////////////////////////////////////////////////////////////////

  def formatNodeFeatures(sTree: Array[sampleTree.treeNode], nodeExpr: Array[Int], nodeSupport: Array[Int], nodePossible: Array[Int], origins: Array[Int]) = {
    ("%5s\t%12s\t%6s\t%6s\t%6s\t%6s\n" format ("ID", "Name","Height", "Expr", "Sup", "Pos")) + 
      sTree.map( s => "%5d\t%12s\t%6d\t%6d\t%6d\t%6s%s" format (s.id+1, s.name, s.height, nodeExpr(s.id), nodeSupport(s.id), if(nodePossible(s.id) == 1) "YES" else if(nodePossible(s.id) == 0) "MAYBE" else "NO", if(origins contains s.id) "*" else "")).mkString("\n");
  }

  /////////////////////////////////////////////////////////////////////////////

  def formatNodePossible(sTree: Array[sampleTree.treeNode], nodePossible: Array[Int]) = {
    nodePossible.map( x => if(x == possible_yes) 'y' else if(x == possible_maybe) 'm' else 'n').mkString
  }

  /////////////////////////////////////////////////////////////////////////////

  def assignOriginToVCGroupGATK(VCs: Array[VCF.VC], sTree: Array[sampleTree.treeNode]) = {
    //GATK doesn't output all the data we need, so we need to write a seperate  algorithm for it.
    // Many parts are the same, however
    val firstVC      = VCs(0)
    val callFlags    = VCs.map( x => x.filter).flatten.distinct
    val legalBase    = callFlags.length == 1 &&  callFlags.head == VCF.pass && firstVC.ref != 'N'
    val observedAlts = VCs.map( x => x.alt.filter( b => b != firstVC.ref && b != VCF.emptyChar )).flatten.distinct

    if (!legalBase){
      VCF.VC(firstVC.contig,
             firstVC.pos,
             VCF.emptyString,
             firstVC.ref,
             observedAlts,
             -1,
             callFlags.filter( x => x != VCF.emptyString && x != VCF.pass),
             Map.empty[String,String])
    } else if(observedAlts.length == 0){
      VCF.VC(firstVC.contig,
             firstVC.pos,
             VCF.emptyString,
             firstVC.ref,
             Array(VCF.emptyChar),
             -1,
             Array("NOSNP"),
             Map.empty[String,String])
    } else if(observedAlts.length > 1){
      VCF.VC(firstVC.contig,
             firstVC.pos,
             VCF.emptyString,
             firstVC.ref,
             observedAlts,
             -1,
             Array("MSNP"),
             Map.empty[String,String])
    } else {

      val alt = observedAlts.head

      // We must assume that if we don't observe the base in the GATK call, then there was NO expression of it.
      // However, we can assume that if a call is MISSING, then there MAY be expression of this base
      val observedAltExprMap = VCs.map( y => (y.ident.split('.')(1).toInt - 1, if(y.alt contains alt) 1 else -1 ) ).toMap
      val filledAltExprMap   = sTree.map( x => (x.id, observedAltExprMap.getOrElse(x.id, -1))).toMap
      val (nodeExpr, nodeSupport, nodePossible, nodeInfoFeatCounts) = calcNodeFeatures(sTree, filledAltExprMap, Map.empty[String, Map[Int, Int]]);

      val origins = determineOrigin(sTree, nodePossible)

      // This is a debug message
      print("###################\n%s:%d: ALT=%c VCs=%d\n%s\n--\n" format (firstVC.contig, firstVC.pos, alt, VCs.length, formatNodeFeatures(sTree, nodeExpr, nodeSupport, nodePossible, origins)))

      val originExpr = origins.map( o => nodeExpr(o))
      VCF.VC(firstVC.contig,
             firstVC.pos,
             "%s;%d;%s" format (firstVC.contig, firstVC.pos, origins.map( x => (x + 1).toString).mkString(",")),
             firstVC.ref,
             observedAlts,
             (100.0*(origins.map(o => nodeSupport(o)).foldLeft(0)((a,b) => math.max(a,b) ).toFloat / VCs.length)).toInt,
             if(origins.length == 1) Array("PASS") else Array("MO"),
             Map("SO" -> origins.length.toString,
                 "SP" -> origins.map( o => nodeSupport(o).toString).mkString(","),
                 "NS" -> VCs.length.toString,
                 "NN" -> origins.map( o => sTree(o).name).mkString(","),
                 "NI" -> origins.map( _.toString).mkString(",")))

    } 
   
  }

  /////////////////////////////////////////////////////////////////////////////

  def determineOrigin(sTree: Array[sampleTree.treeNode], possible: Array[Int]) = {
    // Find the highest node in the tree?
    // Assumption of one GAIN event, and multiple LOSS events

    def traverseTree(node: Int) : Array[Int] = {
      val sampleTree.treeNode(id, name, children, isLeaf, height) = sTree(node);
      val p = possible(node)

      if(p >= possible_maybe){
        Array(id)
      } else if(children.length == 0){
          Array.empty[Int]
      } else {
        children.map(c => traverseTree(c)).flatten
      }
    }

    traverseTree(sTree.length - 1)
  }

  /////////////////////////////////////////////////////////////////////////////

//  def cmpVCGroup(VCs: Array[VCF.VC], groupName: String, nSamples: Int, thresholds: Array[Int]) = {
//
//    val firstVC   = VCs(0)
//    val ref       = firstVC.ref
//    val validVCs  = VCs.filter(x => (x.filter contains "PASS") or (x.filter contains "."));
//
//    val observedAlts    = validVCs.map( _.alt.filter( _ != ref)).flatten.distinct
//    val MSNP = observedAlts.length > 1
//    val NALT = observedAlts.length == 0
//
//    val (alts,passThresholds,callFlags) = if(!MSNP && NALT) {
//      val alt = observedAlts(0)
//      val observedAltExpr = validVCs.map(y => y.info.getOrElse("C%c" format alt, "0").toInt)
//      (Array(ref, alt), thresholds.map(x => observedAltExpr.map( y => y >= x ).foldLeft(0)( (a,b) => a + (if(b) 1 else 0))), Array.empty[String])
//    }
//    else if(MSNP) {
//      (Array('.'), thresholds.map( x => 0), Array("MSNP"))
//    }
//    else{
//      (Array('.'), thresholds.map( x => 0), Array("NALT"))
//    }
//
//    (thresholds zip passThresholds).map( x => VCF.VC(firstVC.contig,
//                                                     firstVC.pos,
//                                                     "%s;%d;%s;%d" format (firstVC.contig, firstVC.pos, groupName, x._1),
//                                                     ref,
//                                                     alts,
//                                                     (x._2.toFloat / nSamples * 100).toInt,
//                                                     if(nSamples == 0) callFlags :+ "NX" else callFlags ,
//                                                     Map("NS" -> validVCs.length.toString, "NP" -> x._2.toString)))
//  }
//
//  /////////////////////////////////////////////////////////////////////////////
//
//  def cmpVCGroupGATK(VCs: Array[VCF.VC], groupName: String, nSamples: Int, thresholds: Array[Int]) = {
//    val firstVC = VCs(0)
//    val nObs    = VCs.filter(x => x.filter contains "PASS");
//    print("%s\t%d\t%d\t%d\n" format (firstVC.contig, firstVC.pos, nObs, nSamples))
//    thresholds.map{ x =>
//      VCF.VC(firstVC.contig,
//             firstVC.pos,
//             "%s;%d;%s;%d" format (firstVC.contig, firstVC.pos, groupName, x),
//             firstVC.ref,
//             firstVC.alt,
//             (nObs.toFloat / nSamples * 100).toInt,
//             Array("PASS"),
//             Map("NS" -> VCs.length.toString, "NP" -> nObs.toString))
//    }
//  }
}
