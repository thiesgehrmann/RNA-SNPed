

//import htsjdk.samtools.SAMFileReader
//import htsjdk.samtools.SAMRecordIterator
//import htsjdk.samtools.SAMSequenceRecord
import htsjdk.samtools._
import scala.io.Source
import java.io._

import org.arabidopsis.interval._
import scala.collection.JavaConversions._

import org.apache.commons.math3.distribution.BinomialDistribution

//import Fasta;

object variantCaller {

  def main(args: Array[String]) = {

    val sjdb_file      = args(0);
    val reference_file = args(1);
    val bam_file       = args(2);
    val sample_id      = args(3);
    val out_file       = args(4);

    println(sjdb_file)
    println(bam_file);
    println(out_file);

    val reference = Fasta.read(reference_file).map(x => (x.description, x.sequence)).toMap
    val BAM       = new SAMFileReader(new File(bam_file));
    val sjdb      = Source.fromFile(sjdb_file).getLines();
    

    val counts = new PositionBaseCounts(BAM, reference, sjdb)
    counts.processReads();
    //SNP.writeCalls(counts.callSNPs().flatten, out_file);
    SNP.writeVCF(counts.callSNPs().flatten, out_file, sample_id);

  }

}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

class PositionBaseCounts(val BAM: SAMFileReader, val reference : Map[String,String], sjdb: Iterator[String]) {

  val MINEXPR = 2;
  val MINQUAL = 30;
  val MINFRAC = 0.01;
  val PERROR  = 0.01;
  val MINPVAL = 0.01;
  val TBOUND  = 4; // Boundary around SJDB, INSERT & DELETION

  val bases   = List('A', 'C', 'G', 'T', 'N');
  val baseMap = Map( 'A' -> 0, 'a' -> 0,
                     'C' -> 1, 'c' -> 1,
                     'G' -> 2, 'g' -> 2,
                     'T' -> 3, 't' -> 3,
                     'N' -> 4, 'n' -> 4);

  val seqIDMap = this.BAM.getFileHeader().getSequenceDictionary().getSequences().toList.zipWithIndex.map( (x: (SAMSequenceRecord, Int)) => (x._1.getSequenceName(), x._2) ).toMap
  val IDSeqMap = this.seqIDMap.map(_.swap)

  var insertionTree = this.BAM.getFileHeader().getSequenceDictionary().getSequences().toList.map( x => new IntervalTree());
  var deletionTree  = this.BAM.getFileHeader().getSequenceDictionary().getSequences().toList.map( x => new IntervalTree());
  var sjdbTree      = this.BAM.getFileHeader().getSequenceDictionary().getSequences().toList.map( x => new IntervalTree());
  var baseCountMap  = this.BAM.getFileHeader().getSequenceDictionary().getSequences().toList.map( x => Array.ofDim[Int](x.getSequenceLength(), 5))

  sjdb.foreach(x => insertSJDBLine(x))

  /////////////////////////////////////////////////////////////////////////////

  def insertSJDBLine(sjdbLine: String) = {
    val lineSplit = sjdbLine.split('\t');
    this.sjdbTree(this.seqIDMap(lineSplit(0))).insert(new Interval(lineSplit(1).toInt - this.TBOUND, lineSplit(2).toInt + this.TBOUND))
  }

  /////////////////////////////////////////////////////////////////////////////

  def processReads() = {
    println("Processing reads!");
    // Here I deviate from the normal SCALA style, because I use a mutable type
    for( read <- this.BAM.iterator()){
      this.processRead(read)
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def callSNPs() : List[List[SNP.SNPcall]] = {
    println("Calling SNPs");
    this.BAM.getFileHeader().getSequenceDictionary().getSequences().toList.map( x => (0 to x.getSequenceLength()-1).map( y => this.testBase(x.getSequenceName(), y)).toList )
  }

  /////////////////////////////////////////////////////////////////////////////

  def testBase(contig: String, pos: Int) : SNP.SNPcall = {
    val contigID   = this.seqIDMap(contig)
    //println("%s: %d" format (contig, pos))
    val refBase    = this.reference(contig)(pos)
    val counts     = this.baseCountMap(contigID)(pos)
    val total      = counts.fold(0){(a,b) => a + b}
    val totalACGT  = total - counts(4); //all minus counts of N
    //val fractions  = counts.map( x => x.toFloat / total)
    //val fractionsACGT = counts.dropRight(1).map( x => x.toFloat / totalACGT)
    //val nonRef     = totalACGT - this.baseCountMap(contigID)(pos)(this.baseMap(refBase));


      // Determine flags
    val noExpr = totalACGT == 0;
    val inSJDB = this.sjdbTree(contigID).searchAll(new Interval(pos, pos)).size() > 0;
    val inDEL  = this.deletionTree(contigID).searchAll(new Interval(pos, pos)).size() > 0;
    val inINS  = this.insertionTree(contigID).searchAll(new Interval(pos, pos)).size() > 0;
    val refMAX = this.baseCountMap(contigID)(pos)(this.baseMap(refBase)) == counts.max;
    val refN   = refBase == 'N';

    var callFlags = scala.collection.mutable.ListBuffer.empty[String];

    if(noExpr){ callFlags += SNP.noExpression}
    if(inSJDB){ callFlags += SNP.nearSpliceJunction }
    if(inDEL){  callFlags += SNP.possibleDeletion }
    if(inINS){  callFlags += SNP.possibleInsertion }
    if(refN){   callFlags += SNP.referenceIsN }

    if( noExpr | inSJDB || inDEL || inINS || refN ) {
      new SNP.SNPcall(contig, pos, false, refBase, '.':: Nil, totalACGT, counts(0), counts(1), counts(2), counts(3), counts(4), callFlags.toList, -1 :: Nil);
    }
    else {
      //val alt = (this.bases zip fractions).dropRight(1).filter( x => x._2 >= this.MINFRAC).map(x => x._1);
      val binom = new BinomialDistribution(totalACGT, this.PERROR);
      val tests = (this.bases zip counts.dropRight(1)).filter(x => x._2 > 0).map(x => (x._1, x._2, 1 - binom.cumulativeProbability(x._2)))
      val alt   = tests.filter( x => x._3  < this.MINPVAL);
      val altb  = alt.map(x => x._1);
      
      val (call, qual) = if(alt.length > 1 || !refMAX){
        // I don't know what is happening, but it seems that the qual function produces -1 sometimes, and I don' know why!!
        val qualcalc = alt.map(x => math.round(-10*math.log10(math.max(x._3, 1e-100f))).toInt ).map( x => if(x == -1) 200  else x)
        (true, qualcalc)
      }
      else {
        (false, -1 :: Nil)
      }

      if(call && !alt.map(x => x._2 > this.MINEXPR).reduceLeft(_ && _)){
        callFlags += SNP.lowExpression
      }

      new SNP.SNPcall(contig, pos, call, refBase, altb, totalACGT, counts(0), counts(1), counts(2), counts(3), counts(4), callFlags.toList, qual)
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def processRead(read: SAMRecord) = {
    val alignmentStart = read.getAlignmentStart();
    val quality = read.getBaseQualityString().toCharArray.map( x => x.toInt - 33)
    val qualitys = read.getBaseQualityString().toCharArray
    val readSeq = read.getReadString()
    val contig  = read.getContig();
    val contigLength = this.reference(contig).size;
    var seq_offset = 0;
    var ref_offset = 0;
    //if ( read.getCigarString() contains "I") {
    //  println("***********************************************");
    //}
    if(!read.getDuplicateReadFlag()){
      for( element <- read.getCigar().getCigarElements() ) {
        element.getOperator() match {
          case CigarOperator.M | CigarOperator.EQ => {
            for (base <- readSeq slice(seq_offset, seq_offset+element.getLength())) {
              val base_pos = alignmentStart+ref_offset-1;
              if(base_pos >= contigLength) {
                 println("SEQ: %s" format (readSeq));
                 println("CIG: %s" format (read.getCigarString()));
                 println("alignmentStart: %d" format (alignmentStart));
                 println("seq: %d, ref: %d" format (seq_offset, ref_offset));
                 println("basePos: %d" format (base_pos));
                 println("contigLength: %d" format (contigLength));
              }
              if(quality(seq_offset) >= this.MINQUAL && base_pos < contigLength){
                this.baseCountMap(this.seqIDMap(contig))(base_pos)(this.baseMap(base)) += 1;
              }
              //println("%s, %d, %d, %s" format (base, seq_offset, ref_offset, this.reference(contig) slice(alignmentStart+ref_offset-1, alignmentStart+ref_offset)))
              seq_offset += 1;
              ref_offset += 1;
            }
          }
          case CigarOperator.I => {
            // Add a four base window around the inserted region on the genome, just to be safe.
            this.insertionTree(this.seqIDMap(contig)).insert(new Interval(alignmentStart+ref_offset-1-this.TBOUND, alignmentStart+ref_offset+this.TBOUND));
            seq_offset += element.getLength();
          }
          case CigarOperator.D | CigarOperator.N => {
            this.deletionTree(this.seqIDMap(contig)).insert(new Interval(alignmentStart+ref_offset-1-this.TBOUND, alignmentStart+ref_offset-1+element.getLength()+this.TBOUND));
            ref_offset += element.getLength();
          }
          case CigarOperator.X => { // DO NOT COUNT THESE!
            ref_offset += element.getLength();
            seq_offset += element.getLength();
          }
          case CigarOperator.S => seq_offset += element.getLength();
          case CigarOperator.H => Nil
          case CigarOperator.P => seq_offset += element.getLength();
          //case _ => DO NOTHING
        }
      }
    }    
  }

  /////////////////////////////////////////////////////////////////////////////

}

///////////////////////////////////////////////////////////////////////////////
