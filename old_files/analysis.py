#!/usr/bin/env python

from ibidas import *

import pandas as pd
import matplotlib.pylab as plt;
import seaborn as sns;
import numpy as np;
import matplotlib as mpl
import scipy.stats

os.sys.argv = [ os.sys.argv[0], '/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3_unique' ];

output_dir = os.sys.argv[1];

###############################################################################
###############################################################################
###############################################################################
   # Correlation between number of reads vs number of mutations found

sample_snp_read_data = Read("%s/sample_SNPs_reads.tsv" % output_dir).Detect() / ("id", "name", "unique", "unique_yes", "unique_maybe", "total", "total_yes", "total_maybe", "nreads")
read_counts = sample_snp_read_data.nreads()
snp_counts  = sample_snp_read_data.total_yes()
second_sequencing_indexes = np.array([0,1,2,3,6,7])


plt.clf(); plt.cla();
ax = sns.regplot(read_counts, snp_counts)
plt.xlabel("Reads mapped to reference");
plt.ylabel("SNPs detected in sample");
plt.scatter(read_counts[second_sequencing_indexes], snp_counts[second_sequencing_indexes], c='r')
plt.ylim([0,3500])
plt.savefig("%s/analysis_figures/correlation_reads_snps.png" % output_dir);
plt.savefig("%s/analysis_figures/correlation_reads_snps.svg" % output_dir);

scipy.stats.spearmanr(sample_snp_read_data.nreads(), sample_snp_read_data.total_yes())

np.corrcoef(sample_snp_read_data.nreads(), sample_snp_read_data.total_yes())


###############################################################################
###############################################################################
###############################################################################
  # correlation to Normalization of number of SNPs

norm_read_counts = read_counts
norm_snp_counts  = sample_snp_read_data.Get(_.total_yes * _.nreads.Max().Cast(float) / _.nreads)()

plt.clf(); plt.cla();
ax = sns.regplot(norm_read_counts, norm_snp_counts)
plt.xlabel("Reads mapped to reference");
plt.ylabel("SNPs in sample (normalized)");
plt.scatter(norm_read_counts[second_sequencing_indexes], norm_snp_counts[second_sequencing_indexes], c='r')
plt.ylim([0,3500])
plt.savefig("%s/analysis_figures/correlation_reads_snps_norm.png" % output_dir);
plt.savefig("%s/analysis_figures/correlation_reads_snps_norm.svg" % output_dir);

np.corrcoef(norm_read_counts, norm_snp_counts)

###############################################################################
###############################################################################
###############################################################################
  # Correlation to the number of covered reads...
sample_min5  = Read('%s/covered_regions.min5.csv' % output_dir).Detect() / ('id', 'nreads')

sample_snp_covb_data = sample_snp_read_data.Get(_.id, _.total_yes) | Match(_.id) | sample_min5

covb_counts = sample_snp_covb_data.nreads()
covb_snp_counts  = sample_snp_covb_data.total_yes()

plt.clf(); plt.cla();
ax = sns.regplot(read_counts, snp_counts)
plt.xlabel("Number of bases with at least 5 reads");
plt.ylabel("SNPs detected in sample");
plt.scatter(covb_counts[second_sequencing_indexes], covb_snp_counts[second_sequencing_indexes], c='r')
plt.ylim([0,3500])
plt.savefig("%s/analysis_figures/correlation_covb_snps.png" % output_dir);
plt.savefig("%s/analysis_figures/correlation_covb_snps.svg" % output_dir);

np.corrcoef(covb_counts, covb_snp_counts)

###############################################################################
###############################################################################
###############################################################################
  # correlation to Normalization with covered reads

covb_norm_read_counts = sample_snp_covb_data.nreads()
covb_norm_snp_counts  = sample_snp_covb_data.Get(_.total_yes * _.nreads.Max().Cast(float) / _.nreads)()

plt.clf(); plt.cla();
ax = sns.regplot(norm_read_counts, norm_snp_counts)
plt.xlabel("Number of bases with at least 5 reads");
plt.ylabel("SNPs in sample (normalized)");
plt.scatter(covb_norm_read_counts[second_sequencing_indexes], covb_norm_snp_counts[second_sequencing_indexes], c='r')
plt.ylim([0,3500])
plt.savefig("%s/analysis_figures/correlation_reads_snps_norm_covb.png" % output_dir);
plt.savefig("%s/analysis_figures/correlation_reads_snps_norm_covb.svg" % output_dir);

np.corrcoef(covb_norm_read_counts, covb_norm_snp_counts)

###############################################################################
###############################################################################
###############################################################################
  # Examine the VAF of the SNPs

VCF = Read('%s/varcall.binom.origins.annotated.vcf' % output_dir)
VCF_NMTL = Read('%s/varcall.binom.origins.nomatingtype.vcf' % output_dir)
VAF = VCF.Get(_.ad[0].Cast(float) / _.dp).Detect()()

plt.cla(); plt.clf();
sns.distplot(VAF)
plt.xlim([0,1])
plt.xlabel("Variant Allele Frequency")
plt.ylabel("Frequency")
plt.savefig("%s/analysis_figures/VAF_hist.png" % output_dir);
plt.savefig("%s/analysis_figures/VAF_hist.svg" % output_dir);

###############################################################################
###############################################################################
###############################################################################
  # Average number of SNPs present in each sample

np.mean([ sum(m) for m in zip(*[ [True if x == 'y' else False for x in y[:46]]  for y in VCF.np() ])])

###############################################################################
###############################################################################
###############################################################################
  # Print vaf as scatterplot

plt.cla(); plt.clf();
ad = VCF.Get(_.ad[0])();
dp = VCF.Get(_.dp)()
(x,y,s,c) = Rep(zip(ad,dp)).GroupBy((_.f0, _.f1)).Get(_.f0[0], _.f1[0] - _.f0[0], (_.f1.Count()+1).Log10()*20 / 's', (_.f0[0].Cast(float) / _.f1[0]) / 'c')()
plt.scatter(x, y, c=c, cmap=mpl.cm.Reds, s=s);
ax = plt.gca()
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')
ax.set_xlim([1, x.max()])
ax.set_ylim([1, dp.max()])
plt.xlabel("Alternative base count");
plt.ylabel("Reference base count");
plt.savefig("%s/analysis_figures/VAF_scatter.png" % output_dir);
plt.savefig("%s/analysis_figures/VAF_scatter.svg" % output_dir);

###############################################################################
###############################################################################
###############################################################################
  #Support and alternative base expression

plt.cla(); plt.clf();
ad = VCF.Get(_.ad[0])();
sp = VCF.Get(_.sp[0])();
(x,y,s) = Rep(zip(ad,sp)).GroupBy((_.f0, _.f1)).Get(_.f0[0], _.f1[0], _.f1.Count())()
plt.scatter(x, y, s=s);
ax = plt.gca()
ax.set_xscale("log", nonposx='clip')
#ax.set_yscale("log", nonposy='clip')
ax.set_xlim([1, x.max()])
ax.set_ylim([1, y.max()])
plt.xlabel("Alternative base count");
plt.ylabel("Number of samplues with expression of SNP");
plt.savefig("%s/analysis_figures/Expression_support_group_scatter.png" % output_dir);
plt.savefig("%s/analysis_figures/Expression_support_group_scatter.svg" % output_dir);

###############################################################################
###############################################################################
###############################################################################
  # Average expression and support
plt.cla(); plt.clf();
ad = VCF.Get(_.ad[0])();
sp = VCF.Get(_.sp[0])();
(x,y,yerr) = Rep(zip(ad,sp)).GroupBy(_.f1).Get(_.f1, (_.f0.Cast(float) / _.f1).Mean(), (_.f0.Cast(float) / _.f1).Std())()
plt.errorbar(x,y, yerr=yerr);
#sns.regplot(x, y);
ax = plt.gca()
#ax.set_xscale("log", nonposx='clip')
#ax.set_yscale("log", nonposy='clip')
ax.set_xlim([0, x.max()+1])
ax.set_ylim([0, y.max() + yerr.max()])
plt.ylabel("Average alternative base count per sample");
plt.xlabel("Number of samples with expression of SNP");
plt.savefig("%s/analysis_figures/Expression_support_scatter.png" % output_dir);
plt.savefig("%s/analysis_figures/Expression_support_scatter.svg" % output_dir);

###############################################################################
###############################################################################
###############################################################################
  # Look at SNP conversion rates

hetero_conv = VCF[_.alt.Count() > 1][_.ref != _.alt].Get(_.ref, _.alt, _.nn[0]).Flat()[_.nn != '70s_wildtype'].GroupBy((_.ref, _.alt)).Get(_.ref[0], _.alt[0], _.alt.Count() / 'count').Sort(_.ref, _.alt)
homo_conv   = VCF[_.alt.Count() == 1][_.ref != _.alt].Get(_.ref, _.alt, _.nn[0]).Flat()[_.nn != '70s_wildtype'].GroupBy((_.ref, _.alt)).Get(_.ref[0], _.alt[0], _.alt.Count() / 'count').Sort(_.ref, _.alt)

###############################################################################

from scipy.stats import chi2_contingency

  # All genes vs mating type genes
#mat = [ [ 25, 11782 ], [71568, 39091968] ]
#chi2_contingency(mat)

  # Only the hydrophobins
#mat = [ [ 25, 11782 ], [1, 7379] ]
#chi2_contingency(mat)

###############################################################################
###############################################################################
###############################################################################
  # Mutation rate estimateion

#def calculate_mut_rate(D, n_bases, generations):
#  nodes_to_exclude = [ 'second_experiment',
#                       'ku80_knockout',
#                       'vegetative_induced.2',
#                       'primordia.1',
#                       'vegetative_mycelium.2',
#                       'primordia.2',
#                       'vegetative_induced.1',
#                       'seed_vegetative_mycelium',
#                       'vegetative_mycelium.1',
#                       'seed_primordia',
#                       'seed_vegetative_induced',
#                       'double_bri1_knockout',
#                       'double_c2h2_knockout',
#                       'double_fst3_knockout',
#                       'double_fst4_knockout',
#                       'double_gat1_knockout',
#                       'double_hom1_knockout',
#                       'double_hom2_knockout',
#                       'double_wc1_knockout',
#                       'double_wc2_knockout',
#                       '70s_wildtype' ];
#  D = D.Detect() / ('node_id', 'node_name', 'unique', 'uy', 'um', 'total', 'ty', 'tm')
#  D = D[~_.node_name.In(nodes_to_exclude)]
#  return D.Get(_.node_name, _.unique, _.unique.Cast(float).Each(lambda x: x/(n_bases*generations)).Cast(float) /'mrate').Get(_.mrate.Mean(), _.mrate.Std())
##edef

#n_bases = 39163626 - 71658

  # background mutation rate
#calculate_mut_rate(Read('%s/mutation_rate_samples.tsv' % output_dir), n_bases, 200)

  # Sliding window
#calculate_mut_rate(Read('%s/mutation_rate_slidingWindow.tsv' % output_dir), 20000*27, 200)

  # Mating type loci
#calculate_mut_rate(Read('%s/mutation_rate_mtl.tsv' % output_dir), 71568, 200)

###############################################################################
###############################################################################
###############################################################################
  # Normalized SNP counts
nodes_to_exclude = [ 'ku80_knockout',
                     'double_bri1_knockout',
                     'double_c2h2_knockout',
                     'double_fst3_knockout',
                     'double_fst4_knockout',
                     'double_gat1_knockout',
                     'double_hom1_knockout',
                     'double_hom2_knockout',
                     'double_wc1_knockout',
                     'double_wc2_knockout',
                     '70s_wildtype' ];

def fill_tree_reads(D, T, R, op=np.max):
  Rtree       = T.Detect()
  tree        = [ (x, y, [] if z == '-' else [ int(c) for c in z.split(',')]) for (x,y,z) in zip(*Rtree()) ]
  sample_snp_read_data = R.Detect()
  SNPs_reads  = sample_snp_read_data.Get(_.id, _.nreads).Cast(int, int)
  SNPs_counts = D.Get(_.f0, _.f1, _.f2).Cast(int, str, int) / ('id', 'name', 'nsnps')

  SNP_counts_reads      = (SNPs_counts | Match(_.id, jointype='full', merge_same='equi') | SNPs_reads)
  SNP_counts_reads_list = zip(*SNP_counts_reads())
  SNP_counts_reads_fill = []
  for (i,name,nsnps,nreads) in SNP_counts_reads_list:
    SNP_counts_reads_fill += [( i, name, nsnps, nreads if len(tree[i-1][2]) == 0 else op([SNP_counts_reads_fill[c-1][3] for c in tree[i-1][2]]))]
  #efor

  SNP_counts_reads = Rep(SNP_counts_reads_fill).Cast(int,str,int,int) / ('id', 'name', 'nsnps', 'nreads')
  
  return SNP_counts_reads;
#edef

def calculate_mut_rate_norm(D, T, R, n_bases, generations, op=np.max, nodes_to_exclude=nodes_to_exclude):
  SNP_counts_reads = fill_tree_reads(D,T,R, op)
  SNP_counts_reads = SNP_counts_reads[~_.name.In(nodes_to_exclude)];
  SNP_counts_norm  = SNP_counts_reads.Get(_.id, _.name, _.nreads,
                                          _.nsnps,
                                          (_.nsnps.Cast(float) * (float(SNP_counts_reads.nreads.Min()())/_.nreads) / 'min'), 
                                          (_.nsnps.Cast(float) * (float(SNP_counts_reads.nreads.Max()())/_.nreads) / 'max'),
                                          (_.nsnps.Cast(float) * (float(SNP_counts_reads.nreads.Mean()())/_.nreads)) / 'mean',
                                          (_.nsnps.Cast(float) * (float(SNP_counts_reads.nreads.Median()())/_.nreads) / 'median'))
  SNP_counts_norm = SNP_counts_norm.Get(_.id, _.name, _.nreads,
                                        _.min, (_.min.Cast(float) / (n_bases*generations)) / 'min_mrate',
                                        _.max, (_.max.Cast(float) / (n_bases*generations)) / 'max_mrate',
                                        _.mean, (_.mean.Cast(float) / (n_bases*generations)) / 'mean_mrate',
                                        _.median, (_.median.Cast(float) / (n_bases*generations)) / 'median_mrate')


  return SNP_counts_norm;
#edef

sample_snps  = Read('%s/mutation_rate_samples.tsv' % output_dir)
mtl_snps     = Read('%s/mutation_rate_mtl.tsv' % output_dir)
tree         = Read('%s/../tree_with_wildtype.tsv' % output_dir)
sample_reads = Read("%s/sample_SNPs_reads.tsv" % output_dir).Get(_.f0,_.f8).Detect() / ('id', 'nreads')
sample_min5  = Read('%s/covered_regions.min5.csv' % output_dir).Detect() / ('id', 'nreads')

import numpy as np
import scipy as sp
import scipy.stats

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h
#edef

n_bases = 39163626 - 71658

#snp_reads_norm = calculate_mut_rate_norm(sample_snps, tree, sample_reads, n_bases, 200, op=np.max)
#snp_min5_norm  = calculate_mut_rate_norm(sample_snps, tree, sample_min5, n_bases, 200, op=np.max)

# Per haploid genome!!
min5_n_bases = (sample_min5.nreads.Max()()-71658) * 2
snp_min5_norm  = calculate_mut_rate_norm(sample_snps, tree, sample_min5, min5_n_bases, 200, op=np.max)
mean_confidence_interval(snp_min5_norm.max_mrate())

snp_min5_norm_all = calculate_mut_rate_norm(sample_snps, tree, sample_min5, min5_n_bases, 200, op=np.max, nodes_to_exclude=[])

plt.cla(); plt.clf();
sns.distplot(snp_min5_norm_all[~_.name.In(nodes_to_exclude)].max(), bins=40, kde=False)
plt.xlabel('Number of SNPs')
plt.ylabel('Number of samples')
plt.xlim([0, 1000]);
plt.savefig('%s/analysis_figures/snps_dist_ku80.svg' % output_dir)
plt.savefig('%s/analysis_figures/snps_dist_ku80.png' % output_dir)

#sns.distplot(snp_min5_norm_all[~_.name.HasPattern('ouble')].max())
#sns.distplot(snp_min5_norm_all[_.name
#plt.xlabel('Number of SNPs')
#plt.ylabel('Number of samples')
#plt.xlim([0, 1000]);
#plt.savefig('%s/analysis_figures/snps_dist_ku80.svg' % output_dir)
#plt.savefig('%s/analysis_figures/snps_dist_ku80.png' % output_dir)


min5_n_bases_mtl   = 71658 * 2
snp_min5_norm_mtl = calculate_mut_rate_norm(mtl_snps, tree, sample_min5, min5_n_bases_mtl, 200, op=np.max)
mean_confidence_interval(snp_min5_norm_mtl.max_mrate())

plt.clf(); plt.cla();
ax = sns.regplot(snp_min5_norm.nreads(), snp_min5_norm.mean())
plt.xlabel("Reads mapped to reference");
plt.ylabel("SNPs detected in sample");
plt.savefig("%s/analysis_figures/correlation_reads_snps_norm.png" % output_dir);
plt.savefig("%s/analysis_figures/correlation_reads_snps_norm.svg" % output_dir);

np.corrcoef(sample_snp_read_data.nreads(), sample_snp_read_data.total_yes())
mean_confidence_interval(snp_min5_norm.max_mrate())

from scipy.stats import ttest_ind
from scipy.stats import ttest_1samp
from scipy.stats import ks_2samp
from scipy.stats import ranksums
from scipy.stats import mannwhitneyu
# Just the ku80 knockout
x = calculate_mut_rate_norm(sample_snps, tree, sample_min5, n_bases, 200, op=np.max, nodes_to_exclude=nodes_to_exclude[1:])
ku80_r  = x[_.name == 'ku80_knockout'].max.Max()()
other_r = x[_.name != 'ku80_knockout'].max()
ttest_1samp(other_r, ku80_r)

# All double knockouts

x = calculate_mut_rate_norm(sample_snps, tree, sample_min5, n_bases, 200, op=np.max, nodes_to_exclude=nodes_to_exclude[:1] + nodes_to_exclude[-1:])
double_r = x[_.name.HasPattern('double')].max()
other_r  = x[~_.name.HasPattern('double')].max()
ttest_ind(double_r, other_r, equal_var=False)
ks_2samp(double_r, other_r)



sample_snp_read_norm = (snp_min5_norm.Get(_.id, _.name, _.max) | Match(_.id, merge_same='equi', jointype='right') | sample_snp_read_data.Get(_.id, _.nreads))



plt.clf(); plt.cla();
ax = sns.regplot(read_counts, snp_counts)
plt.xlabel("Reads mapped to reference");
plt.ylabel("SNPs detected in sample");
plt.scatter(read_counts[second_sequencing_indexes], snp_counts[second_sequencing_indexes], c='r')
plt.savefig("%s/analysis_figures/correlation_reads_snps_corrected.png" % output_dir);
plt.savefig("%s/analysis_figures/correlation_reads_snps_corrected.svg" % output_dir);


###############################################################################
###############################################################################
###############################################################################
  # SNP Windows
chr_sizes = [ ('scaffold_1', 4463208),
              ('scaffold_2', 3875761),
              ('scaffold_3', 3537448),
              ('scaffold_4', 3492607),
              ('scaffold_5', 3342965),
              ('scaffold_6', 2560162),
              ('scaffold_7', 2428502),
              ('scaffold_8', 2413352),
              ('scaffold_9', 2127050),
              ('scaffold_10', 1964478),
              ('scaffold_11', 1833194),
              ('scaffold_12', 1497036),
              ('scaffold_13', 1212600),
              ('scaffold_14', 1149858),
              ('scaffold_15', 710001),
              ('scaffold_16', 609140),
              ('scaffold_17', 608737),
              ('scaffold_18', 527346),
              ('scaffold_19', 94025),
              ('scaffold_20', 70536),
              ('scaffold_21', 65238),
              ('scaffold_22', 41516),
              ('scaffold_23', 39147),
              ('scaffold_24', 3667),
              ('scaffold_25', 2805) ]
chr_sizes = Rep(chr_sizes) / ('seqname', 'length')
windows   = Read("%s/windowSizes.nowildtype.het.bedgraph" % output_dir, skiprows=2).Detect() / ('seqname', 'start', 'end', 'w10000')
wtwindows = Read("%s/windowSizes.wildtype.het.bedgraph" % output_dir, skiprows=2).Detect() / ('seqname', 'start', 'end', 'w10000')

import matplotlib
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from matplotlib import gridspec

def draw_windows(windows, chr_sizes, field, out_prefix, spacing=1, width=.75):

  print field

  cmap = matplotlib.cm.get_cmap('Reds')
  norm = matplotlib.colors.Normalize(vmin=0, vmax=windows.Get(field).Max()())
  cmf  = cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba

  color_boxes = (windows | Match(_.seqname) | chr_sizes).Get(_.seqname.Each(lambda x: x.split('_')[1]).Cast(int), _.start.Cast(float) / _.length, _.end.Cast(float) / _.length, field)

  max_chr = color_boxes.seqname.Max()() + 1
  width   = width / max_chr;

  plt.cla(); plt.clf();
  fig  = plt.figure(figsize=(12,5))
  gs   = gridspec.GridSpec(1, 2, width_ratios=[19, 1]) 
  ax   = plt.subplot(gs[0])
  cbax = plt.subplot(gs[1])

  # Now adding the colorbar
  ax.grid(False)

  cb = mpl.colorbar.ColorbarBase(cbax, cmap=cmap, norm=norm)

  for (h, s, e, v) in zip(*color_boxes[_.Get(field) > 0].Sort(field, descend=False)()):
    ax.add_patch(Rectangle((s, float((max_chr-h)*spacing)/max_chr), e-s, width, facecolor=cmf(v), edgecolor="none", linewidth=0))
  #efor

  for (c,l) in zip(*chr_sizes.Get(_.seqname.Each(lambda y: y.split('_')[1]).Cast(int), _.length)()):
    ax.add_patch(Rectangle((0, float((max_chr-c)*spacing)/max_chr ), 1.0, width, fill=False, edgecolor="black"))
  #efor

  plt.xlim([0,1])
  plt.ylim([0,1])
  plt.savefig("%s.png" % out_prefix, transparent=True);
  plt.savefig("%s.svg" % out_prefix, transparent=True);
  plt.figure(figsize=(8,5))
#edef

#draw_windows(windows, chr_sizes, "w1000", '%s/chr_window_1000' % output_dir)
#draw_windows(windows, chr_sizes, "w5000", '%s/chr_window_5000' % output_dir)
draw_windows(windows, chr_sizes, "w10000", '%s/chr_window_10000' % output_dir)


#draw_windows(wtwindows, chr_sizes, "w1000", '%s/chr_wt_window_1000' % output_dir)
#draw_windows(wtwindows, chr_sizes, "w5000", '%s/chr_wt_window_5000' % output_dir)
draw_windows(wtwindows, chr_sizes, "w10000", '%s/chr_wt_window_10000' % output_dir)


###############################################################################
###############################################################################
###############################################################################
  # Mutation rates across the genome
size = 10000


nsamples = 83 - 1

genome_mutrates = windows.Get(_.seqname, _.start, _.end, _.Get('w%s' % size), _.Get('w%d' % size).Cast(float) / (nsamples*200*(size + _.end - _.start)), 2*(2*size + _.end - _.start)) / ('seqname', 'start', 'end', 'size', 'mutrate', 'length')
genome_mutrates = genome_mutrates.Get(_.seqname, _.start, _.end, _.size, _.mutrate, _.length, (_.length.Cast(float) / genome_mutrates.length.Sum()()) / 'contrib')


def merge_ranges(ranges):
  # From http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
  sorted_intervals = sorted(ranges, key=lambda x: x[0])

  if not sorted_intervals:  # no intervals to merge
    return

# low and high represent the bounds of the current run of merges
  low, high = sorted_intervals[0]

  for iv in sorted_intervals[1:]:
    if iv[0] <= high:  # new interval overlaps current run
      high = max(high, iv[1])  # merge with the current run
    else:  # current run is over
      yield low, high  # yield accumulated interval
      low, high = iv  # start new run

  yield low, high  # end the final run
#edef

def highly_variable_regions(VCF, windows, size, min, nodes_to_exclude=['70s_wildtype']):

  mut_windows = windows.Get(_.seqname, _.start, _.end, _.Get('w%s' % size))
  chr_ranges = zip(*mut_windows[_.Get('w%d' % size) >= min].GroupBy(_.seqname)())

  high_mut_ranges = [ item for sublist in [ [(chr, x[0], x[1], starts, ends) for x in merge_ranges(zip(starts-size, ends+size))] for (chr, starts, ends, counts) in chr_ranges ] for item in sublist ]

  X = Rep([ (chr, start, end, VCF.To(_.nn, Do=_.Get(_.nn[0]).Cast(str))[~_.nn.In(nodes_to_exclude)][_.chrom == chr][_.pos >= start][_.pos <= end].chrom.Shape()()) for (chr, start, end, starts, ends) in high_mut_ranges ])
  return X / ('seqname', 'start', 'end', 'count')
#edef

VCF_nomtl_het = Read('%s/varcall.binom.origins.nomatingtype.het.vcf' % output_dir).Copy()

hvr_nwt = highly_variable_regions(VCF_nomtl_het, windows, 10000, min=20)
hvr_wt  = highly_variable_regions(VCF_nomtl_het, wtwindows, 10000, min=10, nodes_to_exclude=list(set(VCF.nn.Flat().Unique()()) - set(['70s_wildtype'])))
###############################################################################
###############################################################################
###############################################################################
  # Look at phasing
#varcall.binom.origins.phased.vcf

VCF_phased = Read('%s/varcall.binom.origins.nomatingtype.GCR.deleteriousness.phased.vcf' % output_dir)

###
  # Plot the different mutation types

y, c, x, s = VCF_phased.Get(_.dl, _.ad[0].Cast(float) / _.dp).ReplaceMissing().GroupBy((_.dl, _.result)).Get(_.dl[0], _.result[0], _.result.Count()).Get(_.dl.TakeFrom({'S' : 1, 'M' : 2, 'N' : 3}), _.dl.TakeFrom({'S' : 'blue', 'M' : 'green', 'N' : 'red'}), _.result)()

plt.cla(); plt.clf();

plt.scatter(x,y, s=s, c=c.tolist())
plt.xlabel("Variant Allele Frequency")
plt.ylabel("Mutation type");
plt.savefig('%s/analysis_figures/VAF_deleteriousness.svg' % output_dir)
plt.savefig('%s/analysis_figures/VAF_deleteriousness.png' % output_dir)

###############################################################################
###############################################################################
###############################################################################
  #Functional aspects


altsplice = Read("/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/alt_splicing_annots.gff")
iprannots = Read("/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/Schco3_GeneCatalog_proteins_20130812_IPR.tab")
goannots  = Read("/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/Schco3_GeneCatalog_proteins_20130812_GO_CUT_NAMES.tab")
tfdomains = Read("/home/nfs/thiesgehrmann/groups/w/phd/data/fungal_transfac_domains")
domlocs   = Read("/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/ipr_domains.tsv");

asprots    = altsplice[_.feature == 'transcript'].GroupBy(_.parent)[_.feature.Count() >= 2].parent
tfprots    = iprannots[_.iprid.In(tfdomains.f1)].proteinid.Unique()
cp450prots = iprannots[_.iprid == 'IPR001128'].proteinid.Unique()
metabprots = goannots[_.goacc == 'GO:0008152'].proteinid.Unique()
cazyprots  = Read("/home/nfs/thiesgehrmann/groups/w/phd/tasks/cazyme_predicted/schco3/predicted_cazymes.txt") / 'proteinid'


dlvcf  = Read("%s/varcall.binom.origins.nomatingtype.GCR.deleteriousness.vcf" % output_dir)
dlfunc = dlvcf.Get(_.gcr, _.pos, _.nn, _.dl, _.ad[0].Cast(float)/_.dp, ~_.dm.IsMissing(), _.np).Flat()[_.nn != '70s_wildtype'].Copy()
dlfuncnames = dlvcf.Get(_.gcr, _.nn, _.gn, _.dl, _.ad[0].Cast(float)/_.dp).Flat()[_.nn != '70s_wildtype'][_.gn != ''].Copy()

def sig_test(all, a, b):

  import scipy.stats as sstats #chi2_contingency as c2

  all = set(all)
  a   = set(a) & all
  b   = set(b) & all

  na_nb = len(all - (a | b))
  na_yb = len(b - a)
  ya_nb = len(a - b)
  ya_yb = len(a & b)

  print len(all), len(a), len(b)
  print [ [ ya_yb, ya_nb], [na_yb, na_nb]]
  return sstats.chi2_contingency([ [ ya_yb, ya_nb], [na_yb, na_nb] ])

#edef

[ [ dlfunc[_.dl == e][_.gcr.In(x)].gcr.Unique().Shape()() for x in [tfprots, cp450prots, metabprots, cazyprots] ] for e in ['S', 'M', 'N']]

[ [ dlfunc[_.dm][_.dl == e][_.gcr.In(x)].gcr.Unique().Shape()() for x in [tfprots, cp450prots, metabprots, cazyprots] ] for e in ['S', 'M', 'N']]

[ [ sig_test(iprannots.proteinid.Unique()(), x(), dlfunc[_.dl == e].gcr.Unique()()) for x in [tfprots, cp450prots, metabprots, cazyprots] ] for e in ['S', 'M', 'N']]

dlfunc[_.gcr.In(tfprots)].Shape()
dlfunc[_.gcr.In(metabprots)].Shape()
dlfunc[_.gcr.In(cp450prots)].Shape()
dlfunc[_.gcr.In(cazyprots)].Shape()

###############################################################################
###############################################################################
###############################################################################
  # How many of these deleterious SNPs occur before the last exon ends?

last_exon = domlocs.Detect().GroupBy(_.proteinid).Get(_.proteinid.Cast(str), _.Get(_.strand[0], _.starts.Array(), _.ends.Array()).Each(lambda s,b,e: max(e) if s == '+' else min(b)).Cast(int) /'loc', _.strand[0])

dlfunc_exon = (dlfunc | Match(_.gcr, _.proteinid, jointype='inner', merge_same='equi') | last_exon)
dlfunc_exon = dlfunc_exon[_.Get(_.pos, _.loc, _.strand).Each(lambda p, l, s: p < l if s == '+' else p > l).Cast(bool)]

[ [ dlfunc_exon[_.dl == e][_.gcr.In(x)].gcr.Unique().Shape()() for x in [tfprots, cp450prots, metabprots, cazyprots] ] for e in ['S', 'M', 'N']]

[ [ sig_test(iprannots.proteinid.Unique()(), x(), dlfunc_exon[_.dl == e].gcr.Unique()()) for x in [tfprots, cp450prots, metabprots, cazyprots] ] for e in ['S', 'M', 'N']]


###############################################################################
###############################################################################
###############################################################################
  # Examine the probabilities of neutral, missense and nonsense mutation

codon_usage_table = Read('/home/nfs/thiesgehrmann/groups/w/phd/tasks/codon_usage/schco3.tsv').Detect()

bases = ['t', 'c', 'a', 'g'];
codons = [a+b+c for a in bases for b in bases for c in bases];
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
codon_table = dict(zip(codons, amino_acids));

d1_effect = Rep([ (codons[y],[ codons[x] for x in range(len(codons)) if len([a for (a,b) in zip(codons[x], codons[y]) if a != b]) == 1]) for y in range(len(codons))]) / ('codon', 'mut_codon');
d1_effect = (d1_effect | Match(_.codon) | codon_usage_table).Get(_.codon, _.mut_codon, _.aa, _.mut_codon.Each(lambda x: codon_table[x]).Cast(str) / 'mut_aa',  _.counts, _.percentage)
d1_effect = d1_effect.Get(_.codon, _.mut_codon, _.aa, _.mut_aa, _.Get(_.aa, _.mut_aa).Each(lambda x,y: 'S' if x == y else 'N' if y == '*' else 'M').Cast(str) / 'effect', _.counts, _.percentage).Flat()

p_s = d1_effect[_.effect == 'S'].percentage.Sum()() / d1_effect.percentage.Sum()()
p_m = d1_effect[_.effect == 'M'].percentage.Sum()() / d1_effect.percentage.Sum()()
p_n = d1_effect[_.effect == 'N'].percentage.Sum()() / d1_effect.percentage.Sum()()
m, n, s = dlfunc.Get(_.dl, _.gcr).GroupBy(_.dl).To(_.gcr, Do=_.Count()).Sort(_.dl).gcr()
observed = np.array([ s, m, n ]);
total_mut = observed.sum()
expected = np.array([ p_s, p_m, p_n]) * total_mut

import scipy.stats as sstats

sstats.chisquare(observed, expected)


###############################################################################
###############################################################################
###############################################################################
  # Look at UV signature mutations

#VP = VCF.Get(_.chrom, _.pos, _.ref, _.alt, _.nn[0])

#(VP | Stack | VP.To(_.pos, Do=_.Each(lambda x: x-1).Cast(int))).GroupBy(_.pos)[_.ref.In(['C', 'T'])][_.chrom.Count() > 1][_.ref[0] == _.ref[1]][_.nn[0] == _.nn[1]][_.alt.In(['C', 'T'])][_.alt != _.ref[0]][_.alt.Flat().Count() > 1]

###############################################################################
###############################################################################
###############################################################################
  # PCA of samples
from sklearn.decomposition import PCA

D = VCF.np.Each(lambda s: [ 1 if x == 'y' else 0 if x == 'm' else -1 for x in s ]).Detect()().transpose()

pca_w = PCA(n_components=2).fit(D)
pca_m = pca_w.transform(D)

plt.cla(); plt.clf()
cm = plt.cm.get_cmap('RdYlBu')
plt.scatter(pca_m[:,0], pca_m[:,1], c=range(pca_m.shape[0]), vmin=0, vmax=pca_m.shape[0], s=35, cmap=cm)

for (i, (x,y)) in enumerate(pca_m):
  plt.text(x, y, '%d' % (i+1), horizontalalignment='center', verticalalignment='center', fontsize=5);
#efor

plt.xlabel("First component (%d%%)"  % (pca_w.explained_variance_ratio_[0] * 100));
plt.ylabel("Second component (%d%%)" % (pca_w.explained_variance_ratio_[1] * 100));

plt.savefig("%s/PCA_samples.png" % output_dir);
plt.savefig("%s/PCA_samples.svg" % output_dir);


###############################################################################
###############################################################################
###############################################################################

  # Format the expression data labels as necessary


def fdr_gen(p, alpha, type='bh'):

  nt     = len(p);  # Number of tests
  ps   = sorted(p);
  indx = [ i[0] for i in sorted(enumerate(p), key=lambda x:x[1]) ];

  if type == 'bhy':
    cm    = sum([ 1.0/float(i) for i in xrange(nt)] );         
    klist = [ (float(i+1)/(float(nt) * cm)) for i in xrange(nt) ];
  else:
    klist = [ (float(i+1)/float(nt)) for i in xrange(nt) ];
  #fi

    # Adjust pvalues to qvalues
  q = [ ps[i] / klist[i] for i in xrange(nt)];
    # Fix pvalues larger than 1
  q = [ qi if qi < 1.0 else 1.0 for qi in q ];

    # Monotonicity
  qm = [];
  prev_v = q[0];
  for v in q:
    qm.append(max(prev_v, v));
    prev_v = qm[-1];
  #efor

    # get back to original sorting
  qrs = [0] * nt;
  for i in xrange(nt):
    qrs[indx[i]] = qm[i];
  #efor

  return qrs;
#edef

###############################################################################

expr = Read('/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/rnaseq_second_sample_expression.tsv').Detect().To(_.test_id, Do=_.Cast(str)).Copy()

expr_id_mappings = { 'veg_myc_1_normalized'   : 1,
                     'veg_myc_2_normalized'   : 2,
                     'veg_ind_1_normalized'   : 3,
                     'veg_ind_2_normalized'   : 4,
                     '_4_8_1_normalized'      : 5,
                     '_4_8_2_normalized'      : 6,
                     'primordia_1_normalized' : 7,
                     'primordia_2_normalized' : 8,
                     '_4_8_3_normalized'      : 9,
                     '_4_8_4_normalized'      : 10,
                     'bri1_1_normalized'      : 11,
                     'bri1_2_normalized'      : 12,
                     'bri1_3_normalized'      : 13,
                     'bri1_4_normalized'      : 14,
                     'c2h2_1_normalized'      : 15,
                     'c2h2_2_normalized'      : 16,
                     'c2h2_3_normalized'      : 17,
                     'c2h2_4_normalized'      : 18,
                     'fst3_1_normalized'      : 19,
                     'fst3_2_normalized'      : 20,
                     'fst3_3_normalized'      : 21,
                     'fst3_4_normalized'      : 22,
                     'fst4_1_normalized'      : 23,
                     'fst4_2_normalized'      : 24,
                     'fst4_3_normalized'      : 25,
                     'fst4_4_normalized'      : 26,
                     'gat1_1_normalized'      : 27,
                     'gat1_2_normalized'      : 28,
                     'gat1_3_normalized'      : 29,
                     'gat1_4_normalized'      : 30,
                     'hom1_1_normalized'      : 31,
                     'hom1_2_normalized'      : 32,
                     'hom1_3_normalized'      : 33,
                     'hom1_4_normalized'      : 34,
                     'hom2_1_normalized'      : 35,
                     'hom2_2_normalized'      : 36,
                     'hom2_3_normalized'      : 37,
                     'hom2_4_normalized'      : 38,
                     'wc1_1_normalized'       : 39,
                     'wc1_2_normalized'       : 40,
                     'wc1_3_normalized'       : 41,
                     'wc1_4_normalized'       : 42,
                     'wc2_1_normalized'       : 43,
                     'wc2_2_normalized'       : 44,
                     'wc2_3_normalized'       : 45,
                     'wc2_4_normalized'       : 46}

Ed = dict([ (p[0], np.array(p[1:])) for p in zip(*expr.Get(_.test_id, *[ a[0] for a in sorted(expr_id_mappings.items(), key=lambda x: x[1])])())])

dlfunc_expr = []
for (gcr, pos, nn, dl, result, dm, node_p) in zip(*dlfunc()):
  gene_expr = Ed[gcr]
  if gene_expr.sum() == 0:
    continue;
  #fi
  node_p = np.array(list(node_p[:len(expr_id_mappings)]))
  gene_expr_y = gene_expr[np.where(node_p == 'y')]
  gene_expr_n = gene_expr[np.where(node_p == 'n')]
  #stat, pval = ks_2samp(gene_expr_y, gene_expr_n);
  stat, pval = ranksums(gene_expr_y, gene_expr_n);
  #stat, pval = ttest_ind(gene_expr_y, gene_expr_n, equal_var=False)
  dlfunc_expr.append((gcr, dl, gene_expr_y, gene_expr_n, stat, pval))
#efor

import statsmodels
from statsmodels.sandbox.stats.multicomp import multipletests
dlfunc_expr     = zip(*dlfunc_expr)
#dlfunc_expr[-1] = fdr_gen(np.array(dlfunc_expr[-1]), 0.05, 'bhy')
pval = np.array(dlfunc_expr[-1])
pval[np.isnan(pval)] = 1.0
discard_1, qval, discard_2, discard_3 = multipletests(pval, 0.05, 'fdr_bh')
dlfunc_expr     = dlfunc_expr + [ qval ]
dlfunc_expr     = zip(*dlfunc_expr)
dlfunc_expr     = Rep(dlfunc_expr) / ('gcr', 'dl', 'y', 'n', 'stat', 'pval', 'qval')

###############################################################################
###############################################################################
###############################################################################
  # Do SNPs next to neighboring genes influence expression?
  
gff_file = '/tudelft.net/staff-groups/ewi/insy/DBL/thiesgehrmann/w/phd/data/schco3/cleaned_gene_model_proteinid.gff3.gff'
G        = Read(gff_file)

from bx.intervals.intersection import IntervalTree

itrees = {};
for (seqname, id, start, end, strand) in zip(*G[_.feature == 'mRNA'].Get(_.seqname, _.id, _.start, _.end, _.strand)()):
  if seqname not in itrees:
    itrees[seqname] = IntervalTree()
  #fi
  if strand == '+':
    r_s = start - 500
    r_e = start
  else:
    r_s = end
    r_e = end + 500
  itrees[seqname].add(r_s, r_e, id)
#efor

intergenic_vcf = VCF[_.gid.IsMissing()].Get(_.chrom, _.pos, _.nn, _.np, _.ad.Cast(float) / _.dp).Copy()
assoc_expr = []
for (seqname, pos, nn, node_p, vaf) in zip(*intergenic_vcf()):
  assoc_genes = itrees[seqname].find(pos, pos);
  print seqname, pos
  if(len(assoc_genes) == 0):
    continue;
  #fi
  for gene in assoc_genes:
    gene_expr = Ed[gene];
    node_p = np.array(list(node_p[:len(expr_id_mappings)]))
    gene_expr_y = gene_expr[np.where(node_p == 'y')]
    gene_expr_n = gene_expr[np.where(node_p == 'n')]
    #stat, pval = ks_2samp(gene_expr_y, gene_expr_n);
    #stat, pval = ranksums(gene_expr_y, gene_expr_n);
    stat, pval = ttest_ind(gene_expr_y, gene_expr_n, equal_var=False)
    assoc_expr.append((seqname, pos, nn, gene, vaf, gene_expr_y, gene_expr_n, stat, pval))
  #efor
#efor

assoc_expr     = zip(*assoc_expr)
#dlfunc_expr[-1] = fdr_gen(np.array(dlfunc_expr[-1]), 0.05, 'bhy')
pval = np.array(assoc_expr[-1])
pval[np.isnan(pval)] = 1.0
discard_1, qval, discard_2, discard_3 = multipletests(pval, 0.05, 'fdr_bh')
assoc_expr     = assoc_expr + [ qval ]
assoc_expr     = zip(*assoc_expr)
assoc_expr     = Rep(assoc_expr) / ('chrom', 'pos', 'nn', 'gid', 'vaf', 'y', 'n', 'stat', 'pval', 'qval')

  # Significant values:
sig_assoc = assoc_expr[_.qval < 0.01].Show().Flat(_.nn).To(_.y, Do=_.Mean()).To(_.n, Do=_.Mean())

CDS = Read('/tmp/cuffdiff_sig.tsv').Detect().Copy()
difftest_names = dict(enumerate([ "wc1-1", "wc1-2", "wc2-1", "wc2-2", "hom1-1", "hom1-2", "hom2-1", "hom2-2","fst3-1", "fst3-2", "fst4-1", "fst4-2", "bri1-1", "bri1-2", "gat1-1", "gat1-2", "c2h2-1", "c2h2-2", "4.8-1", "4.8-2", "4.8-4", "4.8-5", "4.8-6", ]))
CDS = CDS.To(_.f4, Do=_.TakeFrom(difftest_names) / 'cond1').To(_.f5, Do=_.TakeFrom(difftest_names) / 'cond2').To(_.f0, Do=_.Cast(str) / 'gid').Copy()

(sig_assoc | Match(_.gid, _.gid) | CDS).GroupBy(_.chrom, _.pos)

sig_assoc_explained = (sig_assoc | Match(_.gid, _.gid, jointype='left') | CDS).GroupBy((_.chrom, _.pos))

sig_assoc_explained = sig_assoc_explained.Get(_.chrom[0], _.pos[0], _.nn[0], _.gid[0], _.vaf[0], _.y[0], _.n[0], _.qval[0], _.Get(_.cond1, _.cond2).Array().Each(lambda x,y: ['%s/%s' % (a,b) for (a,b) in zip(x,y)] ))
