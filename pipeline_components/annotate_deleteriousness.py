#!/usr/bin/env python

import numpy as np

import biu

import re
import sys

def apply_mutation(gff, fasta, vcf_record):
  oSeq = fasta[v.CHROM]
  oFasta = biu.formats.Fasta([ oSeq ])
  vAlt = [ a.sequence for a in vcf_record.ALT if a.sequence != v.REF][0]
  vSeq = oSeq[:vcf_record.POS] + vAlt + oSeq[vcf_record.POS+1:]
  vFasta = biu.formats.Fasta([ vSeq ])
  
  oTrans = str(gff.seq(vcf_record.INFO["GCR"], oFasta).translate())
  vTrans = str(gff.seq(vcf_record.INFO["GCR"], vFasta).translate())
  
  #for i,(ob,vb) in enumerate(zip(oTrans,vTrans)):
  #  if ob != vb:
  #    print(i, vcf_record.POS, ob, '->', vb)
  
  if len(oTrans) > len(vTrans):
    return 'N'
  elif oTrans != vTrans:
    return 'M'
  else:
    return 'S'
  #fi
#edef

###############################################################################

def help(arg0):
  print("%s: Annotate a deleteriousness to SNPs in a VCF file" % arg0)
  print(" NOTE: IT ONLY ANNOTATES VARIANTS WITH A 'GCR=TRANSCRIPT_ID' tag in them. This can be added with rnasnped filterVCF annotate.")
  print("usage: %s <input_vcf> <output_vcf> <gff_file> <fa_file>")
  print("")
  print(" input_vcf:  The input VCF file")
  print(" output_vcf: The location of the output VCF")
  print(" gff_file:   A GFF file containing gene descriptions")
  print(" fa_file:    The genome sequence in multifasta format")
  print("")
  print(" Each SNP in a coding region will be annotated with an extra field in the info field, dl=X, where X can be either")
  print("  S: Synonymous")
  print("  M: Missense")
  print("  N: Nonsense")

###############################################################################

if __name__ == '__main__':

  if (len(sys.argv) < 5 or sys.argv[1] == 'help'):
    help(sys.argv[0])
    sys.exit(1)
  #fi

  input_vcf  = sys.argv[1];
  output_vcf = sys.argv[2];
  gff_file   = sys.argv[3];
  fa_file    = sys.argv[4];

  biu.config.settings.setDebugState(False)

  vcf   = biu.formats.VCF(input_vcf)
  gff   = biu.formats.GFF3(gff_file)
  fasta = biu.formats.Fasta(fa_file)

  with open(output_vcf, 'w') as ofd:
    ofd.write('##fileformat=VCFv4.1\n');
    ofd.write('##fileDate=20160818\n');
    ofd.write('##source=RNASNEP\n');
    ofd.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples">\n');
    ofd.write('##INFO=<ID=SO,Number=.,Type=Integer,Description="SNP origin nodes">\n');
    ofd.write('##INFO=<ID=SP,Number=.,Type=Integer,Description="Support for each origin node">\n');
    ofd.write('##INFO=<ID=AD,Number=.,Type=Integer,Description="Number of reads observed for alternative base at each origin node">\n');
    ofd.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n');
    ofd.write('##INFO=<ID=CA,Number=1,Type=Float,Description="Count of A">\n');
    ofd.write('##INFO=<ID=CC,Number=1,Type=Float,Description="Count of C">\n');
    ofd.write('##INFO=<ID=CG,Number=1,Type=Float,Description="Count of G">\n');
    ofd.write('##INFO=<ID=CT,Number=1,Type=Float,Description="Count of T">\n');
    ofd.write("##INFO=<ID=SAD,Number=.,Type=Integer,Description=\"Sample AD\">\n")
    ofd.write("##INFO=<ID=SCA,Number=.,Type=Integer,Description=\"Sample CA\">\n")
    ofd.write("##INFO=<ID=SCC,Number=.,Type=Integer,Description=\"Sample CC\">\n")
    ofd.write("##INFO=<ID=SCG,Number=.,Type=Integer,Description=\"Sample CG\">\n")
    ofd.write("##INFO=<ID=SCT,Number=.,Type=Integer,Description=\"Sample CT\">\n")
    ofd.write("##INFO=<ID=SDP,Number=.,Type=Integer,Description=\"Sample DP\">\n")
    ofd.write('##INFO=<ID=NN,Number=.,Type=String,Description="Origin node names">\n');
    ofd.write('##INFO=<ID=NI,Number=.,Type=Integer,Description="Origin node IDS in tree">\n');
    ofd.write('##INFO=<ID=NP,Number=1,Type=String,Description="Tree possibility description">\n');
    ofd.write('##INFO=<ID=GN,Number=1,Type=String,Description="Gene Name">\n');
    ofd.write('##INFO=<ID=GID,Number=1,Type=String,Description="Gene ID">\n');
    ofd.write('##INFO=<ID=GCR,Number=1,Type=String,Description="Gene coding region, gene ID">\n');
    ofd.write('##INFO=<ID=MTL,Number=1,Type=String,Description="Mating type locus ID">\n');
    ofd.write('##INFO=<ID=ISO,Number=1,Type=String,Description="Isoform transcript ID">\n');
    ofd.write('##INFO=<ID=GI,Number=1,Type=String,Description="Gene annotation">\n');
    ofd.write('##INFO=<ID=DL,Number=1,Type=String,Description="Deleteriousness of the SNP">\n');
    ofd.write('##INFO=<ID=VP,Number=1,Type=String,Description="Allele of SNP">\n');
    ofd.write("##INFO=<ID=DM,Number=1,Type=String,Description=\"Domain annotation\">\n");
    ofd.write('##FILTER=<ID=NX,Description="No Expression">\n');
    ofd.write('##FILTER=<ID=LX,Description="Low Expression">\n');
    ofd.write('##FILTER=<ID=MSNP,Description="SNP multiple alternative alleles">\n');
    ofd.write('##FILTER=<ID=NALT,Description="No SNP at this location">\n');
    ofd.write('##FILTER=<ID=NO,Description="No suitable origin could be found">\n');
    ofd.write('##FILTER=<ID=MO,Description="Multiple origins were found :S">\n');
    ofd.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    for v in vcf:
      chrom = v.CHROM
      pos   = v.POS
      ident = v.ID
      ref   = v.REF
      alt   = ','.join([ a.sequence if hasattr(a, 'sequence') else '.' for a in v.ALT ])
      qual  = v.QUAL
      filt  = "PASS" if len(v.FILTER) == 0 else ';'.join(v.FILTER)
      info  = v.INFO
      if 'GCR' in v.INFO:
        info['DL'] = apply_mutation(gff, fasta, v)
      #fi
      info  = ';'.join([ f.upper() + '=' + ','.join([ str(x) for x in np.array([v.INFO[f]]).flatten().tolist() ]) for f in v.INFO])
      info = re.sub('\[|\]|-', '', info)

      ofd.write('\t'.join([ str(f) for f in [chrom, pos, ident, ref, alt, qual, filt, info] ]) + '\n')
    #efor
  #ewith
#fi
