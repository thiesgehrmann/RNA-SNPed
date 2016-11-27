object call2GFF {

  def main(args: Array[String]) = {

    // Usage: call2GFF input_calls.tsv output_gff.gff

    val callFile = args(0);
    val gffFile = args(1);
    val varAStream = SNP.readCalls(callFile);

    GFF.write(varAStream.filter( _.call).map( x => SNP.call2gff(x)).toList, gffFile)
  }
}
