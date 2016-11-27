#!/bin/sh

infiles=`echo $1 | tr ',' ' '`
outfile=$2;
shift; shift;
tmpfile=$(mktemp /home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3_unique/vcf_merge_sort.XXXXXX)

cat $infiles > $tmpfile

sort $@  -t\t -k1,1 -k2,2n $tmpfile > $outfile
rm $tmpfile
