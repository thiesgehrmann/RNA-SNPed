#!/bin/sh

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";

# If you have a special scala installation directory, please either add it to your PATH variable, or edit the following lines
#scala="/home/nfs/thiesgehrmann/scala-2.11.7/bin/scala"
#scalac="/home/nfs/thiesgehrmann/scala-2.11.7/bin/scalac"
scala="scala"
scalac="scalac"

scala_class_path="./:lib/htsjdk-2.5.0-SNAPSHOT-all.jar:lib/log4j-1.2.17.jar:lib/commons-math3-3.6.1.jar:lib/org.arabidopsis.interval.jar"
alias scala="$scala -classpath $scala_class_path";
alias scalac="$scalac -classpath $scala_class_path";

###############################################################################

echo "Building rnasnped..."
mkdir -p rnasnped
scalac -d ./ Fasta.scala GFF.scala SNP.scala cmpCalls.scala variantCaller.scala sampleTree.scala VCF.scala assignOrigin.scala filterVCF.scala phaseSNPs.scala slidingWindowSNPs2.scala mainClass.scala

###############################################################################

echo "Constructing jar..."
cat > rnasnped.mf << EOF
Main-Class: rnasnped.mainClass
Class-Path: $( echo $scala_class_path | cut --complement -d: -f1 | tr ':' '\n' | sed -e 's/^.*$/ & /')
 lib/scala-library.jar 
 lib/scala-parser-combinators_2.11-1.0.4.jar
EOF

jar -cmf rnasnped.mf rnasnped.jar $(find rnasnped | grep -e '.*class$')

