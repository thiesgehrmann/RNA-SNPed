#!/bin/sh

infiles=`echo $1 | tr ',' ' '`
outfile=$2;
shift; shift;
tmpfile=$(mktemp /tmp/vcf_merge_sort.XXXXXX)

cat $infiles > $tmpfile

sort $@  -t\t -k1,1 -k2,2n $tmpfile > $outfile
rm $tmpfile
