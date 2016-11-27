import java.io.{PrintWriter, FileWriter}

object GFF{

  case class GFFEntry(seqname: String,
                      source: String,
                      feature: String,
                      start: Int,
                      end: Int,
                      score: String,
                      strand: String,
                      frame: String,
                      attr: String)

  /////////////////////////////////////////////////////////////////////////////

  def write(entries: TraversableOnce[GFFEntry], outfile: String) = {

    val outfd = new PrintWriter(new FileWriter(outfile, false));

    for(entry <- entries){
      outfd.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" format (entry.seqname, entry.source, entry.feature, entry.start, entry.end, entry.score, entry.strand, entry.frame, entry.attr))
    }
   outfd.close();
  }

 ///////////////////////////////////////////////////////////////////////////// 

}
