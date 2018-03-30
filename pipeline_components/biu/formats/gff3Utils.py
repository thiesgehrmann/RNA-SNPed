from .. import utils

import errno
import os
import csv
import gzip
from collections import namedtuple

import gzip

###############################################################################

idField     = lambda e: e.attr["ID"] if "ID" in e.attr else e.attr["Name"] if "Name" in e.attr else None
parentField = lambda e: e.attr["Parent"] if "Parent" in e.attr else None
nameField   = lambda e: e.attr["Name"] if "Name" in e.attr else None
seqIDField  = lambda e: e.seqid

###############################################################################

class GFF3Entry(object):

  seqid = None
  source = None
  feature = None
  start = None
  end = None
  score = None
  phase = None
  attr = None

  def __init__(self, row, idField=idField, parentField=parentField, nameField=nameField, **kwargs):
  
    self.__idField = idField
    self.__parentField = parentField
    self.__nameField = nameField

    def attrsplit(attr):
      spl = attr.split('=')
      if len(spl) == 1:
        return (attr, None)
      elif len(spl) > 2:
        return (spl[0], '='.join(spl[1:]))
      else:
        return (spl[0], spl[1])
      #fi
    #edef
  
    (self.seqid, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase, attr) = row
    if isinstance(attr, dict):
      self.attr = attr
    else:
      self.attr = dict([attrsplit(x.strip()) for x in attr.split(";") ])
    #fi
    self.start = int(self.start)
    self.end   = int(self.end)
  #edef

  @property
  def id(self):
    return self.__idField(self)
  #edef

  @property
  def parent(self):
    return self.__parentField(self)
  #edef

  @property
  def name(self):
    return self.__nameField(self)
  #edef

  def seq(self, fastaObject):
    if self.seqid not in fastaObject:
      utils.error("Sequence '%s' not found in provided fastaObject." % self.seqid)
      return None
    #fi
    if self.seqid not in fastaObject:
      utils.error("Sequence '%s' not found in provided fastaObject." % self.seqid)
      return None
    #fi

    substr = fastaObject[self.seqid][self.start-1:self.end]
    if self.strand == '-':
      substr = substr.revcomp()
    #fi

    return substr
  #edef

  def __attrString(self):
    return ';'.join([ '%s=%s' % (k,v) if (v is not None) else k for (k,v) in self.attr.items() ])

  def outputString(self):
    return "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s" % (self.seqid, self.source, self.feature, self.start, 
                                                   self.end, self.score, self.strand, self.phase,
                                                   self.__attrString() )
  #edef

  def copy(self):
    return GFF3Entry([self.seqid, self.source, self.feature, self.start, self.end, self.score, self.strand, self.phase, self.attr], idField=self.__idField, parentField=self.__parentField, nameField=self.__nameField)
  #edef

  def removeParent(self):
    self.attr = { k : v for (k,v) in self.attr.items() if v != self.parent }
  #edef

  def __str__(self):
    dstr = "GFF3Entry(seqid:%s, source:%s, feature:%s, start:%d, end:%d, score:%s, strand:%s, phase:%s, attr:%s)" % (self.seqid, self.source, self.feature, self.start,
                                                                                                                     self.end, self.score, self.strand, self.phase,
                                                                                                                     self.__attrString() )
    return dstr
  #edef

  def __repr__(self):
    return self.__str__()
#eclass

###############################################################################

class GFF3(object):
  entries = None
  seqids = None
  index = None
  topLevel = None

  __index = None
  __fileName = None

  def __init__(self, data, **kwargs):

    if isinstance(data, str):
      utils.dbm("GFF input source is file.")
      self.__fileName = data
      self.entries = GFF3.read(data, **kwargs)
    else:
      utils.dbm("GFF input source is list of GFF3Entries.")
      self.entries = data
    #fi
    self.seqids  = set([ e.seqid for e in self.entries])
    self.__index, self.topLevel = self._index()
  #edef

  def __iter__(self):
    return self.entries
  #edef

  def __str__(self):
    dstr  = "GFF3 object\n"
    dstr += " Where: %s\n" % (self.__fileName if self.__fileName is not None else hex(id(self)))
    dstr += " Entries: %d\n" % len(self.entries)
    dstr += " Top level statistics:\n"
    for featureType in self.topLevel:
      dstr += "  * %s : %d\n" % (featureType, len(self.topLevel[featureType]))
    #efor
    return dstr
  #edef
      

  def _index(self):
    internal_counter = 0

    idx = {}
    topLevel = {}
    for i, e in enumerate(self.entries):
      ID = e.id
      if ID is None:
        ID = "internal.%d" % internal_counter
        internal_counter += 1
      #fi

      parent = e.parent
      if ID not in idx:
        idx[ID] = [i, [] ]
      else:
        idx[ID] = [i, idx[ID][1] ]
      #fi
      if parent is None:
        if e.feature not in topLevel:
          topLevel[e.feature] = []
        #fi
        topLevel[e.feature].append(ID)
      else:
        if parent not in idx:
          idx[parent] = [ None, [] ]
        #fi
        idx[parent][1].append((i, ID))
      #fi
    #efor
    return idx, topLevel
  #edef

  def __len__(self):
    return len(self.entries)
  #edef

  def seq(self, ID, fastaObject):
    entries = self.getChildren(ID, feature="CDS").entries
    if len(entries) == 0:
      utils.error("Could not find ID '%s'" % ID)
      return None
    #fi
    entries = sorted(entries, key=lambda e: e.start)
    if entries[0].strand == '-':
      entries = entries[::-1]
    #fi

    return sum([ e.seq(fastaObject) for e in entries ])
  #edef

  def getIDEntry(self, ID):
    if ID in self.__index:
      return self.entries[self.__index[ID][0]]
    else:
      return None
    #fi
  #edef

  def getIDChildren(self, ID):
    if ID in self.__index:
      return self.__index[ID][1]
    else:
      utils.error("No children found for '%s'" % ID)
      return []
    #fi
  #edef
      
  def getChildren(self, ID, feature=None, depth=None, containParent=False):

    relEntries = [ self.getIDEntry(ID) ] if containParent else []
    ids = [ (0, c[0], c[1]) for c in self.getIDChildren(ID) ]

    while len(ids) > 0:
      cdepth, cidIndex, cid = ids.pop(0)
      relEntries.append(self.entries[cidIndex])
      if (cid in self.__index) and ( (depth is None) or (cdepth+1 < depth)):
        ids.extend([ (cdepth+1, c[0], c[1]) for c in self.getIDChildren(cid) ])
      #fi
    #ewhile

    if feature is not None:
      relEntries = [ r for r in relEntries if r.feature in feature ]
    #fi

    # Strip the parent tag if the parent == ID
    # "seqid, source, feature, start, end, score, strand, phase, attr"
    if not(containParent):
      E = []
      for e in relEntries:
        if e.parent == ID:
          e = e.copy()
          e.removeParent()
          E.append(e)
        else:
          E.append(e)
        #fi
      #efor
      relEntries = E
    #fi

    return GFF3(data = relEntries)
  #edef  

  def indexByInterval(self):
    try: 
      from intervaltree import Interval, IntervalTree
    
      t = { seqid: IntervalTree() for seqid in self.seqids }
      for e in [ e for e in self.entries if e.type.lower() == "mrna"]:
        t[e.seqid][e.start:e.end] = e
      #efor
      return t
    except ImportError:
      return {}
  #edef

  def getID(self, ID):
    if ID in self.__index:
      return self.entries[self.__index[ID][0]]
    else:
      return None
    #fi
  #edef

  def areTandem(self, id1, id2):
    f1 = self.getID(id1)
    f2 = self.getID(id2)
    if gene1.seqid != gene2.seqid:
      return False
    #fi

    inregion = set([ e[-1].attr["ID"] for e in self.interval[gene1.seqid][min(gene1.end,gene2.end):max(gene1.start,gene2.start)] ])
    if len(inregion - set([gene1.attr["ID"], gene2.attr["ID"]])) == 0:
      return True
    else:
      return False
    #fi
  #edef

  @property
  def dataFrame(self):
    import pandas as pd
    return pd.DataFrame([ (e.seqid, e.source, e.feature, e.start, e.end, e.score, e.phase, e.attr) for e in self.entries ],
                        columns=("seqid", "source", "feature", "start", "end", "score", "phase", "attr") )
  #edef

  def areSameStrand(self, id1, id2):
    f1 = self.getID(id1)
    f2 = self.getID(id2)
    return f1.strand == f2.strand
  #edef

  def write(self, fileName):
    with (gzip.open(fileName, "wb") if fileName[-2:] == "gz" else open(fileName, "w")) as gffFile:
      gffFile.write("#gff-version\t3\n")
      for e in self.entries:
        #gff3Entry = namedtuple("gff3Entry", "seqid, source, feature, start, end, score, strand, phase, attr")
        gffFile.write("%s\n" % e.outputString())
      #efor
    #ewith
  #edef

###############################################################################

  @staticmethod
  def read(filename, skipLines=0, maxLines=None, allowAdditionalColumns=False, **kwargs):
    G = []
    nLines = 0
    with (gzip.open(filename, "rt") if filename[-2:] == "gz" else open(filename, "r")) as gffFile:
      gffReader = csv.reader(gffFile, delimiter="\t", quotechar='"')
      for row in gffReader:
        nLines += 1
        if nLines < skipLines:
          continue
        elif (maxLines is not None) and nLines > maxLines:
          break
        #fi
  
        # Unless we allow it with allowAdditionalColumns, we require exactly 9 columns.
        ncolumns = len(row)
        if (ncolumns < 9) or ((ncolumns > 9) and not(allowAdditionalColumns)):
          continue
        #fi
        G.append(GFF3Entry(row[:9], **kwargs))
      #efor
    #ewith
    return G
  #edef

#eclass
###############################################################################
