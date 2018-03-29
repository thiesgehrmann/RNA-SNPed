#!/usr/bin/env python

import json
import sys
import inspect, os

###############################################################################

__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

if len(sys.argv) < 2:
  print("Error, incomplete arguments to checkInput.py")
  sys.exit(1)
#fi

configFile = sys.argv[1];

errors = []
warnings = []

config = {}
if not(os.path.isfile(sys.argv[1])):
  errors.append("ConfigFile '%s'doesn't exist!"% sys.argv[1])
else:
  try:
    with open(sys.argv[1],"r") as ifd:
      config = json.load(ifd)
    #ewith
  except json.decoder.JSONDecodeError as error:
    errors.append("JSON error: {0}".format(error))
    config = {}
  except:
    errors.append("Could not read JSON file.")
    config = {}
  #etry
#fi
dconfig = json.load(open("%s/defaults.json" % __INSTALL_DIR__,"r"))

###############################################################################

if ("samples" not in config) or (len(config["samples"]) == 0):
  errors.append("No samples defined")
else:
  fields = [ "ID", "group_id", "group_name", "sample_stage" ]
  for sample in config["samples"]:
    sampleData = config["samples"][sample]
    for field in fields:
      if field not in sampleData:
        errors.append("Field '%s' not define for sample '%s'." % ( field, sample))
      #fi
    #efor
    for readnumber in ["R1", "R2"]:
      if (readnumber not in sampleData) or (len(sampleData[readnumber]) == 0):
        errors.append("No '%s' read file defined for sample '%s'." % (readnumber, sample))
      #fi
      if not(os.path.isfile(sampleData[readnumber])):
        errors.append("FASTQ file '%s' for sample '%s' does not exist." % (fastq, sample))
      #fi
    #fi
  #efor
#fi

# Check also the tree format!

# Check that GFF file is present
if "gff_file" not in config:
  errors.append("GFF file not provided in 'gff_file'.")
else:
  if not(os.path.isfile(config["gff_file"])):
    errors.append("GFF file '%s' does not exist." % config["gff_file"])
  #fi
#fi

###############################################################################
    
for error in errors:
  print("ERROR: %s" % error)
#efor

for warning in warnings:
  print("WARNING: %s" % warning)
#efor

if len(errors) > 0:
  sys.exit(1)
#fi

sys.exit(0)

