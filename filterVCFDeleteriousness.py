from ibidas import *

import pandas as pd
import matplotlib.pylab as plt;
import seaborn as sns;
import numpy as np;
import matplotlib as mpl

os.sys.argv = [ os.sys.argv[0], 'input_vcf', 'output_vcf' 'gff_file', 'fa_file']

input_vcf  = os.sys.argv[1];
output_vcf = os.sys.argv[2];
gff_file   = os.sys.argv[3];
fa_file    = os.sys.argv[4];

###############################################################################


bases = ['t', 'c', 'a', 'g'];
codons = [a+b+c for a in bases for b in bases for c in bases];
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
codon_table = dict(zip(codons, amino_acids));

###############################################################################

def revcomp(seq):
  dict = { 'A' : 'T', 'T' : 'A', 'a' : 't', 't' : 'a',
           'C' : 'G', 'G' : 'C', 'c' : 'g', 'g' : 'c',
           'R' : 'Y', 'Y' : 'R', 'r' : 'y', 'y' : 'r',
           'K' : 'M', 'M' : 'K', 'k' : 'm', 'm' : 'k',
           'S' : 'W', 'W' : 'S', 's' : 'w', 'w' : 's',
           'B' : 'V', 'V' : 'B', 'b' : 'v', 'v' : 'b',
           'D' : 'H', 'H' : 'D', 'd' : 'h', 'h' : 'd',
           'N' : 'N', 'X' : 'X', 'n' : 'n', 'x' : 'x' }
  
  return ''.join([ dict[b] for b in seq[::-1] ]);

#edef

###############################################################################

def translate(seq, starts, ends, strand):

  transcript = ''.join([ seq[start-1:end] for (start, end) in zip(starts, ends) ])
  transcript = transcript if strand == '+' else revcomp(transcript)

  aa = '';
  for i in xrange(0, len(transcript), 3):
    cod = transcript[i:i+3].lower();
    aa += codon_table[cod] if (cod in codon_table) else '*';
    if(aa[-1] == '*'):
      break;
    #fi
  #efor
  return aa;

#edef

###############################################################################

def apply_mut_translate(pos, ref, alts, seq, starts, ends, strand, orig_transcript):
  alt = [ a for a in alts if a != ref ][0]
  seq = seq[0:pos] + a + seq[pos+1:]
  new_transcript = translate(seq, starts, ends, strand)
  if(len(new_transcript) < len(orig_transcript)):
    return "N" #Nonsence
  elif(new_transcript != orig_transcript):
    return "M" # Missence
  else:
    return "S" # Synonymous
  #fi
#edef

###############################################################################

def flatten(items, seqtypes=(list, tuple)):
  for i, x in enumerate(items):
    while i < len(items) and isinstance(items[i], seqtypes):
       items[i:i+1] = items[i]
    #ewhile
  #efor
  return items
#edef

###############################################################################

def write_vcf(vcf, outfile):
  outfd = open(outfile, 'w')
  outfd.write('##fileformat=VCFv4.1\n');
  outfd.write('##fileDate=20160818\n');
  outfd.write('##source=RNASNEP\n');
  outfd.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples">\n');
  outfd.write('##INFO=<ID=SO,Number=.,Type=Integer,Description="SNP origin nodes">\n');
  outfd.write('##INFO=<ID=SP,Number=.,Type=Integer,Description="Support for each origin node">\n');
  outfd.write('##INFO=<ID=AD,Number=.,Type=Integer,Description="Number of reads observed for alternative base at each origin node">\n');
  outfd.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n');
  outfd.write('##INFO=<ID=CA,Number=1,Type=Float,Description="Count of A">\n');
  outfd.write('##INFO=<ID=CC,Number=1,Type=Float,Description="Count of C">\n');
  outfd.write('##INFO=<ID=CG,Number=1,Type=Float,Description="Count of G">\n');
  outfd.write('##INFO=<ID=CT,Number=1,Type=Float,Description="Count of T">\n');
  outfd.write('##INFO=<ID=NN,Number=.,Type=String,Description="Origin node names">\n');
  outfd.write('##INFO=<ID=NI,Number=.,Type=Integer,Description="Origin node IDS in tree">\n');
  outfd.write('##INFO=<ID=NP,Number=1,Type=String,Description="Tree possibility description">\n');
  outfd.write('##INFO=<ID=GN,Number=1,Type=String,Description="Gene Name">\n');
  outfd.write('##INFO=<ID=GID,Number=1,Type=String,Description="Gene ID">\n');
  outfd.write('##INFO=<ID=GCR,Number=1,Type=String,Description="Gene coding region, gene ID">\n');
  outfd.write('##INFO=<ID=MTL,Number=1,Type=String,Description="Mating type locus ID">\n');
  outfd.write('##INFO=<ID=ISO,Number=1,Type=String,Description="Isoform transcript ID">\n');
  outfd.write('##INFO=<ID=GI,Number=1,Type=String,Description="Gene annotation">\n');
  outfd.write('##INFO=<ID=DL,Number=1,Type=String,Description="Deleteriousness of the annotation">\n');
  outfd.write('##FILTER=<ID=NX,Description="No Expression">\n');
  outfd.write('##FILTER=<ID=LX,Description="Low Expression">\n');
  outfd.write('##FILTER=<ID=MSNP,Description="SNP multiple alternative alleles">\n');
  outfd.write('##FILTER=<ID=NALT,Description="No SNP at this location">\n');
  outfd.write('##FILTER=<ID=NO,Description="No suitable origin could be found">\n');
  outfd.write('##FILTER=<ID=MO,Description="Multiple origins were found :S">\n');

  reserved_fields = [ 'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter' ]
  info_fields     = [ n for n in vcf.Names if n not in reserved_fields ]
  for vc in zip(*vcf()):
    vcd = dict(zip(vcf.Names, vc));
    outfd.write('\t'.join([vcd['chrom'],
                    str(vcd['pos']),
                    '.'.join(vcd["id"]),
                    vcd["ref"],
                    ','.join(vcd["alt"]),
                    str(vcd["qual"]),
                    "PASS" if len(vcd["filter"]) == 0 else ';'.join(vcd["filter"]),
                    ';'.join([ f.upper() + '=' + ','.join([ str(x) for x in flatten([vcd[f]]) ]) for f in info_fields]).translate(None, "'[]-") ]) + '\n');
  #efor
  outfd.close()
#edef

###############################################################################

VCF = Read(input_vcf)
G   = Read(gff_file);
S   = Read(fa_file)

G_exons = G[_.feature == 'exon'].GroupBy(_.parent).Sort(_.start).Get(_.seqname[0], _.parent, _.start.Array(), _.end.Array(), _.strand[0]);
G_seq   = (G_exons | Match(_.seqname, _.f0) | S).Get(_.seqname, _.parent, _.start, _.end, _.strand, _.Get(_.seq, _.start, _.end, _.strand).Each(lambda seq, starts, ends, strand: translate(seq, starts, ends, strand)) / 'transcript')


VCF_genes_seq = ((VCF.To(_.gcr, Do=_.ReplaceMissing("")) | Match(_.gcr, _.parent) | G_seq[_.parent.In(VCF.gcr.Unique())]).Without(_.seqname, _.parent) | Match(_.chrom, _.f0) | S).Copy()

VCF_dl = VCF_genes_seq.Get( _.Get(_.pos, _.ref, _.alt.Array(), _.seq, _.start, _.end, _.strand, _.transcript).Each(lambda pos, ref, alts, seq, starts, ends, strand, t: apply_mut_translate(pos, ref, alts, seq, starts, ends, strand, t)) / "dl", *VCF.Names).Copy()


write_vcf(VCF_dl.Without(_.samples).Sort((_.chrom, _.pos)), output_vcf)

###############################################################################

