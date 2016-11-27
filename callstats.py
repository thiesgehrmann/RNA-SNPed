#!/usr/bin/env python

import ibidas;
from matplotlib import pylab as plt;
import seaborn as sns

os.sys.argv = [ os.sys.argv[0], '/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3/varcalls.binom.1.tsv', '/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3/callstats.1' ]
os.sys.argv = [ os.sys.argv[0], '/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3/varcmp.1.2.tsv', '/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3/callstats.1,2' ]

SNPcmp = os.sys.argv[1];
out_prefix = os.sys.argv[2];

###############################################################################

def readSNPcmp(file):
  calls = Read(SNPcalls, delimiter='\t')
  for c in [ "altbase", "totalcount", "a", "c", "g", "t", "n" ]:
    calls = calls.To(c, Do=_.Each(lambda x: x.split(';')))
  #efor
  #calls = calls.Cast(str, int, str, str, str, int, int, int, int, int, int, str)

  return calls.Detect().To(_.contig, Do=_.Cast(str)).Copy();
#edef

###############################################################################

def callsIntersection(cmp):
  #PAPB and no DA
  return cmp[_.call == 'PAPB'][~_.callflags.HasPattern('DA')]
#edef

def callsdifference(cmp)

###############################################################################

def getVAF(altbases, A, C, G, T):
  bmap = { 'A': A, 'C': C, 'G' : G, 'T': T };
  
  ab = altbases.split(',');
  a  = float(bmap[ab[0]]);
  b  = float(bmap[ab[1]]);
  
  apb = a + b;
  
  return min( a/apb, b/apb);
#edef

###############################################################################

def getVAFcmp(altbases, A, C, G, T):
  bmap = { 'A': A, 'C': C, 'G' : G, 'T': T };

  ab = altbases.split(',');
  a  = float(bmap[ab[0]]);
  b  = float(bmap[ab[1]]);

  apb = a + b;

  return min( a/apb, b/apb);
#edef

###############################################################################

cmp = readSNPcmp(SNPcmp)

cmp[

calls     = Read(SNPcalls, delimiter='\t')[_.Get(2) == 'T'].Copy();
calls     = calls.Cast(str, int, str, str, str, int, int, int, int, int, int, str)
nolxcalls = calls[~_.callflags.HasPattern('LX')];


homo_calls = calls[_.altbase.Each(lambda x: len(x)   == 1).Cast(bool)]
hetero_calls = calls[_.altbase.Each(lambda x: len(x) == 3).Cast(bool)]

nolx_homo_calls   = nolxcalls[_.altbase.Each(lambda x: len(x) == 1).Cast(bool)]
nolx_hetero_calls = nolxcalls[_.altbase.Each(lambda x: len(x) == 3).Cast(bool)]

def getVAF(altbases, A, C, G, T):
  bmap = { 'A': A, 'C': C, 'G' : G, 'T': T };

  ab = altbases.split(',');
  a  = float(bmap[ab[0]]);
  b  = float(bmap[ab[1]]);

  apb = a + b;

  return min( a/apb, b/apb);
#edef

vaf      = hetero_calls.Get(_.contig, _.pos, Tuple(_.altbase, _.a, _.c, _.g, _.t).Each( lambda x: getVAF(*x)).Cast(float) / 'vaf')
nolx_vaf = nolx_hetero_calls.Get(_.contig, _.pos, Tuple(_.altbase, _.a, _.c, _.g, _.t).Each( lambda x: getVAF(*x)).Cast(float) / 'vaf')

plt.cla(); plt.clf();
sns.distplot(vaf.vaf());
plt.savefig(out_prefix + '.vaf.png');
plt.savefig(out_prefix + '.vaf.svg');

plt.cla(); plt.clf();
sns.distplot(nolx_vaf.vaf());
plt.savefig(out_prefix + '.vaf.nolx.png');
plt.savefig(out_prefix + '.vaf.nolx.svg');


