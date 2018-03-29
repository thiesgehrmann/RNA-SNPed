#!/usr/bin/env python
import biu

import sys 

if __name__ == '__main__':

  if len(sys.argv) < 8:
    print("Generate gene feature files")
    print(' usage ./gene_annots.py "{input_gff}" "{output_names}" "{output_regions}" "{output_coding_regions}" "{params_nameAttr}" "{params_geneFeature}" "{params_cdsFeature}"')
    sys.exit(1)
  #fi

  input_gff             = sys.argv[1]
  output_names          = sys.argv[2]
  output_regions        = sys.argv[3]
  output_coding_regions = sys.argv[4]
  params_nameAttr       = sys.argv[5]
  params_geneFeature    = sys.argv[6]
  params_cdsFeature     = sys.argv[7]

  gff = biu.formats.GFF3(input_gff)
  genes = [ e for e in gff.entries if (e.feature == params_geneFeature) ]
  with open(output_names, "w") as ofd:
    for e in [e for e in genes if (params_nameAttr in e.attr) ]:
      ofd.write("%s\t%d\t%d\t%s\n" % (e.seqid, e.start, e.end, e.attr[params_nameAttr]))
    #efor
  #ewith

  with open(output_regions, "w") as ofd:
    for e in genes:
      ofd.write("%s\t%d\t%d\t%s\n" % (e.seqid, e.start, e.end, e.attr["ID"]))
    #efor
  #ewith

  with open(output_coding_regions, "w") as ofd:
    for e in genes:
      geneID = e.attr["ID"]
      cds = gff.getChildren(geneID, feature=params_cdsFeature).entries
      for cr in cds:
        ofd.write("%s\t%d\t%d\t%s\n" % (cr.seqid, cr.start, cr.end, geneID))
      #efor
    #efor
  #ewith
