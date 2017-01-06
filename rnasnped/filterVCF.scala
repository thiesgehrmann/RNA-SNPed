import scala.collection.JavaConversions._

import scala.io.Source
import java.io._

import scala.collection.mutable.ArrayBuffer
import org.arabidopsis.interval._

package rnasnped{
object filterVCF {

  val operations = Map( "removeRegions"  -> remove _,
                        "isolateRegions" -> isolate _,
                        "annotate"       -> annotate _,
                        "infoMin"        -> infoMin _,
                        "infoMax"        -> infoMax _,
                        "infoStringEq"   -> infoStringEq _,
                        "pass"           -> pass _,
                        "toGFF3"         -> toGFF _,
                        "toDOT"          -> toDOT _,
                        "hetSNP"         -> hetSNP _,
                        "homoSNP"        -> homoSNP _,
                        "rmKOSampleSNP"  -> rmKOSampleSNP _);

  /////////////////////////////////////////////////////////////////////////////

  def main(args: Array[String]) = {

    if( args.length < 3 || args(0) == "help") {
      help()
    } else {

    val vcfFile   = args(0);
    val outFile   = args(1);
    val operation = args(2);
    val operationParams = args.slice(3, args.length);

    val vcf      = VCF.readVCF(vcfFile);
    val vcfFilt = operations(operation)(operationParams, vcf)
    VCF.write(vcfFilt, outFile)

    }

  }

  /////////////////////////////////////////////////////////////////////////////

  def help() = {

    println("filterVCF: perform filtering operations on VCF files from the rnaSNPed toolkit")
    println("Usage: filterVCF <inputVCF> <outputVCF> <operation> <operation_options>")
    println("")
    println(" inputVCF: location of input VCF (- for stdin)")
    println(" outputVCF: location of output VCF (- for stdout)")
    println(" operation:")
    println("")
    println("  help")
    println("  Show this help message")
    println("")
    println("  removeRegions <regionsFile>")
    println("  Remove items from the VCF file that lie within regions specified in <regionsFile>")
    println("    regionsFile: A tab delimited file with three columns, chromosome_name, start, and end")
    println("")
    println("  isolateRegions <regionsFile>")
    println("  Remove items from the VCF file that DO NOT lie within regions specified in <regionsFile>")
    println("    regionsFile: A tab delimited file with three columns, chromosome_name, start, and end")
    println("")
    println("  annotate <regionsFile> <fieldName>")
    println("  Add an extra field in the info column of the VCF with name <fieldName>, and a value specified in <regionsFile>")
    println("    regionsFile: A tab delimited file with FOUR columns, chromosome_name, start, end, label")
    println("    fieldName:   The fieldname to give this annotation")
    println("")
    println("  infoMin <fieldName> <value>")
    println("  Threshold the VCF file based on <fieldName> in the info column")
    println("    fieldName: The fieldname to use")
    println("    value:     The minimum value to allow in outputVCF")
    println("")
    println("  infoMax <fieldName> <value>")
    println("  Threshold the VCF file based on <fieldName> in the info column")
    println("    fieldName: The fieldname to use")
    println("    value:     The maximum value to allow in outputVCF")
    println("")
    println("  infoStringEq <fieldName> <value> [complement]")
    println("  Remove elements in the VCF file whose info field value in <fieldName> is not equal to <value>")
    println("    fieldName:  The fieldname to use")
    println("    value:      The value to compare against")
    println("    complement: any non-empty value given here will produce the complementary output")
    println("")
    println("  pass")
    println("  Remove all variants which did not pass a previous filtering step")
    println("")
    println("  toGFF3 <outputFile>")
    println("  Output a GFF3 version of this VCF file")
    println("    outputFile: location of output GFF3 file (- for stdout)")
    println("")
    println("  toDOT <treeFile> <graphName> <outFile>")
    println("  Given a tree structure (the same used for assignOrigin), produce a DOT formatted graph with SNP counts")
    println("    treeFile: the location of the tree file")
    println("    graphName: The name to give the graph (e.g. SNP)")
    println("    outFile:   The location of the output (- for stdout)")
    println("")
    println("  hetSNP")
    println("  Return only heterozygous SNPs")
    println("")
    println("  homoSNP")
    println("  Return only homozygous SNPs")
    println("")
    println("  rmKOSampleSNP <regionsFile>")
    println("  If some of your samples are knockouts, SNPs identified in those regions could be a result of the knockout")
    println("  This command removes them if the SNP lies within the region of a knocked out gene AND is present in the knockout sample.")
    println("   regionsFile: a tab separated file with FOUR fields: chromosome_name, start, end and (comma separated list of sample ids that have this gene knocked out.)")
  }

  /////////////////////////////////////////////////////////////////////////////
    
  def readRegions(regionsFile: String) = {
    var regionsMap     = scala.collection.mutable.Map.empty[String,IntervalTree]
    var regionNamesMap = scala.collection.mutable.Map.empty[(String,Int,Int),String]
    val stream = if(regionsFile == "-"){ io.Source.stdin } else {  io.Source.fromFile(regionsFile) }

    stream.getLines.filter(x => x.length != 0 && x(0) != '#').foreach{
      line =>
        val lsplit = line.split('\t');
        val contig = lsplit(0);
        val low    = lsplit(1).toInt;
        val high   = lsplit(2).toInt;
        val name   = lsplit(3);
        
        if ( !(regionsMap contains contig)) {
          regionsMap(contig) = new IntervalTree();
        }
        regionsMap(contig).insert(new Interval(low, high))
        regionNamesMap((contig,low,high)) = name;
    }
    (regionsMap.toMap, regionNamesMap.toMap)
  }

  /////////////////////////////////////////////////////////////////////////////

  def remove(params: Array[String], vcf: Iterator[VCF.VC]) = {

    val (regionsMap, regionNamesMap) = readRegions(params(0))

    vcf.filter{
      vc =>
        if(regionsMap contains vc.contig) {
          val presentInRegions = regionsMap(vc.contig).searchAll(new Interval(vc.pos, vc.pos))
          presentInRegions.size() == 0
        } else {
          true
        }
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def isolate(params: Array[String], vcf: Iterator[VCF.VC]) = {
    val (regionsMap, regionNamesMap) = readRegions(params(0))
    val fieldName = if(params.length > 1) params(1) else "RG";
    
    vcf.map{
      vc =>
        if(regionsMap contains vc.contig) {
          val presentInRegions = regionsMap(vc.contig).searchAll(new Interval(vc.pos, vc.pos))
          if(presentInRegions.size() > 0) {
           (true, VCF.VC(vc.contig, vc.pos, vc.ident, vc.ref, vc.alt, vc.qual, vc.filter, vc.info + (fieldName -> presentInRegions.map( r => regionNamesMap((vc.contig,r.getLow(),r.getHigh()))).mkString(",") )))
          } else{
           (false, vc)
          }
        } else {
          (false, vc)
        }
    }.filter( x => x._1).map( x => x._2)
  }

  /////////////////////////////////////////////////////////////////////////////

  def annotate(params: Array[String], vcf: Iterator[VCF.VC]) = {
    val (regionsMap, regionNamesMap) = readRegions(params(0))
    val fieldName = params(1)
    
    vcf.map{
      vc =>
        if(regionsMap contains vc.contig) {
          val presentInRegions = regionsMap(vc.contig).searchAll(new Interval(vc.pos, vc.pos))
          if(presentInRegions.size() > 0) {
           VCF.VC(vc.contig, vc.pos, vc.ident, vc.ref, vc.alt, vc.qual, vc.filter, vc.info + (fieldName -> presentInRegions.map( r => regionNamesMap((vc.contig,r.getLow(),r.getHigh()))).mkString(",") ))
          } else{
           vc
          }    
        } else {
          vc
        }                   
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////

  def infoMin(params: Array[String], vcf: Iterator[VCF.VC]) = {
    val field = params(0);
    val value = params(1).toFloat;
    vcf.filter{ vc => (vc.info contains field) && (vc.info(field).toFloat >= value)}
  }

  /////////////////////////////////////////////////////////////////////////////

  def infoMax(params: Array[String], vcf: Iterator[VCF.VC]) = {
    val field = params(0);
    val value = params(1).toFloat;
    vcf.filter{ vc => (vc.info contains field) && (vc.info(field).toFloat <= value)}
  }

  /////////////////////////////////////////////////////////////////////////////

  def infoStringEq(params: Array[String], vcf: Iterator[VCF.VC]) = {
    val field = params(0);
    val value = params(1);
    val not   = if (params.length > 2) true else false
    vcf.filter{ vc => (vc.info contains field) && (if(!not) (vc.info(field) == value) else (vc.info(field) != value))}
  }

  /////////////////////////////////////////////////////////////////////////////

  def pass(params: Array[String], vcf: Iterator[VCF.VC]) = { vcf.filter(vc => vc.filter contains VCF.pass) }

  /////////////////////////////////////////////////////////////////////////////

  def toGFF(params: Array[String], vcf: Iterator[VCF.VC]) = {
    val outGFF = params(0);
    GFF.write(vcf.map(x => VCF.toGFF(x)), outGFF)
    Iterator.empty
  }

  /////////////////////////////////////////////////////////////////////////////

  def toDOT(params: Array[String], vcf: Iterator[VCF.VC]) = {

    val tree_file = params(0);
    val graphName = params(1);
    val outFile   = params(2);

    val sTree  = sampleTree.readTree(tree_file)

    var nodeYesCounts        = collection.mutable.ArrayBuffer.fill(sTree.length)(0)
    var nodeMaybeCounts      = collection.mutable.ArrayBuffer.fill(sTree.length)(0)
    var totalNodeYesCounts   = collection.mutable.ArrayBuffer.fill(sTree.length)(0)
    var totalNodeMaybeCounts = collection.mutable.ArrayBuffer.fill(sTree.length)(0)


    def traverseTreeAdd(nodeID: Int, vc: VCF.VC, root: Boolean): Unit = {
      val nodePossible = if(vc.info contains "NP") vc.info("NP")(nodeID) else 'y'
      if(root) {
        if(nodePossible == 'y') {
          nodeYesCounts(nodeID)   += 1
        } else if(nodePossible == 'm'){
          nodeMaybeCounts(nodeID) += 1
        }
      } else {
        if(nodePossible == 'y') {
          totalNodeYesCounts(nodeID)   += 1
        } else if(nodePossible == 'm'){
          totalNodeMaybeCounts(nodeID) += 1
        }   
      }
      sTree(nodeID).children.foreach(c => traverseTreeAdd(c, vc, false))
    }

    vcf.foreach{ vc =>
      val nodeID = vc.info("NI").toInt - 1;
      traverseTreeAdd(nodeID, vc, true)
    }

    val outfd = if(outFile == "-") new BufferedWriter(new OutputStreamWriter(System.out, "utf-8")) else new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "utf-8"))

    outfd.write("digraph %s{\n" format graphName);
    outfd.write("  rankdir=LR;\n");
    outfd.write("  ranksep=1.0;\n");

    sTree.foreach{ node =>
      val nodeID = node.id
      outfd.write("  n%d [ shape=record, label=\"{ <enter> %s | { Unique: %d | Yes: %d | Maybe: %d } | { Total: %d | <exit> Yes: %d | Maybe: %d}}\"]\n" format
        (nodeID+1,
         node.name,
         nodeYesCounts(nodeID) + nodeMaybeCounts(nodeID),
         nodeYesCounts(nodeID),
         nodeMaybeCounts(nodeID),
         nodeYesCounts(nodeID) + nodeMaybeCounts(nodeID) + totalNodeYesCounts(nodeID) + totalNodeMaybeCounts(nodeID),
         nodeYesCounts(nodeID) + totalNodeYesCounts(nodeID),
         nodeMaybeCounts(nodeID) + totalNodeMaybeCounts(nodeID)))
      outfd.write("/*%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d*/\n" format
        (nodeID+1,
         node.name,
         nodeYesCounts(nodeID) + nodeMaybeCounts(nodeID),
         nodeYesCounts(nodeID),
         nodeMaybeCounts(nodeID),
         nodeYesCounts(nodeID) + nodeMaybeCounts(nodeID) + totalNodeYesCounts(nodeID) + totalNodeMaybeCounts(nodeID),
         nodeYesCounts(nodeID) + totalNodeYesCounts(nodeID),
         nodeMaybeCounts(nodeID) + totalNodeMaybeCounts(nodeID)))
     }
     
     sTree.foreach( node => node.children.foreach( c => outfd.write("  n%d:exit -> n%d:enter\n" format (node.id+1, c+1)) ) )

     outfd.write("}");
     outfd.close();

     Iterator.empty
  }

  /////////////////////////////////////////////////////////////////////////////

  def hetSNP(params: Array[String], vcf: Iterator[VCF.VC]) = {
    vcf.filter(vc => vc.alt.length > 1)
  }

  /////////////////////////////////////////////////////////////////////////////

  def homoSNP(params: Array[String], vcf: Iterator[VCF.VC]) = {
    vcf.filter(vc => (vc.alt.length == 1) && !(vc.alt contains vc.ref))
  }  
  /////////////////////////////////////////////////////////////////////////////

  def rmKOSampleSNP(params: Array[String], vcf: Iterator[VCF.VC]) = {
    val (regionsMap, regionNamesMap) = readRegions(params(0))
    val sampleIDs = regionNamesMap.map{ case (k,v) => (k, v.split(",").map(x => x.toInt))}
    vcf.filter{ vc =>
      if (!regionsMap.contains(vc.contig)){
        true
      } else {
        val presentInRegions = regionsMap(vc.contig).searchAll(new Interval(vc.pos, vc.pos))
        if(presentInRegions.length == 0){
          true
        } else {
          val infoPossible = vc.info("NP")
          val rSampleIDs = presentInRegions.map( r => sampleIDs(vc.contig, r.getLow(), r.getHigh())).flatten.toSet
          val pass = presentInRegions.zipWithIndex.map{ case (p, i) => if(p == 'y' && !rSampleIDs.contains(i)) true else false}.foldLeft(false)(_ | _)
          pass
        }                                                                                                                                                                                                                    
      }
    }
  }

}
}
