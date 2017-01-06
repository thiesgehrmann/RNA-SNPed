import scala.collection.JavaConversions._

import scala.io.Source
import java.io._

import org.arabidopsis.interval._

/* NOTE: DOES NOT WORK FOR MULTIPLE WINDOWSIZES PARAMETERS
 * However, it might be sufficient for your downstream analyses
 * It produces lists of overlapping regions
 * If you want to merge these regions later, it is ok, but they should not be
 * compared individually*/

package rnasnped{
object slidingWindowSNPs2 {

  def main(args: Array[String]) = {

    if (args.length < 3 || args(0) == "help") {
      help()
    } else {
      perform(args)
    }

  }

  /////////////////////////////////////////////////////////////////////////////

  def help() = {
    println("slidingWindowSNPs2: Count SNPs at each base in a sliding window")
    println("Usage: slidingwindowSNPs2 <vcfFile> <outFile> <windowSize1>[,windowsize2[,windowsize3...]]")
    println("")
    println(" vcfFile: The input VCF file")
    println(" outFile: The output file")
    println(" windowsize: Integer, the size of the window you want to consider.")
    println("              ATTN: a provided window size of 5 will count 5 bases upstream and 5 bases downstream")
  }

  /////////////////////////////////////////////////////////////////////////////

  def perform(args: Array[String]) = {
    
    val vcfFile     = args(0);
    val outFile     = args(1);
    val windowSizes = args(2).split(",").map(x => x.toInt);
    
    val vcf    = VCF.readVCF(vcfFile);
    val groups = new groupVCIteratorContig(vcf.buffered).map( x => x.toArray);

    var SNPdict = scala.collection.mutable.Map.empty[String,IntervalTree]

    val windowBoundaries  = groups.filter(_.length > 0).map{ contigGroup =>
      var contig = contigGroup(0).contig;
      var minp = contigGroup(0).pos;
      var maxp = -1;
      var nvc  = 0
      SNPdict(contig) = new IntervalTree()
      contigGroup.foreach{ vc =>
        maxp = math.max(maxp, vc.pos);
        minp = math.min(minp, vc.pos);
        nvc  += 1
        SNPdict(contig).insert(new Interval(vc.pos, vc.pos))
      }
      (contig, minp, maxp, nvc)
    }

    val outfd = if(outFile == "-") new BufferedWriter(new OutputStreamWriter(System.out, "utf-8")) else new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile), "utf-8"))

    outfd.write("#contig\tstart\tend\t%s\n" format windowSizes.map(x => "w%d" format x).mkString("\t"))
    outfd.write("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");

    windowBoundaries.foreach{ case (contig, min, max, n) =>
      (min to max).foldLeft( List(Tuple4("", 0, 0, windowSizes.map(x => 0)))){ case (pList, pos) =>
        val ws = windowSizes.map( w => SNPdict(contig).searchAll(new Interval(pos-w, pos+w)).length )
        val pElem = pList.head
        /* This reduction doesn't work properly if there is more than one window size */
        if( (pElem._4 zip ws).map( x => x._1 == x._2).foldLeft(true)(_ && _)){
          //The elements are the same
          Tuple4(pElem._1, pElem._2, pos, ws) :: pList.drop(1)
        }
        else{
          Tuple4(contig, pos, pos, ws) :: pList
        }
      }.filter( x => x._4.foldLeft(0)(math.max(_, _)) > 0).foreach{
        case (contig, s, e, ws) =>
          outfd.write("%s\t%d\t%d\t%s\n" format (contig, s, e, ws.map(x => x.toString).mkString("\t")))
        }
    }

    outfd.close()
    
  }

  /////////////////////////////////////////////////////////////////////////////

  class groupVCIteratorContig(iter: BufferedIterator[VCF.VC]) extends VCF.groupVCIterator(iter) {
    override def elemCompare(vcA: VCF.VC, vcB: VCF.VC) : Boolean = {
      vcA.contig == vcB.contig
    }
  }

  /////////////////////////////////////////////////////////////////////////////

}
}
