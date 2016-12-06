#!/bin/usr/sh

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
SCRIPTDIR="/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/RNA-Snep"
###############################################################################
#  CONFIGURATION

PICARD=/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/picard/picard-tools-2.5.0/picard.jar
var_caller_dir=$SCRIPTDIR/varCaller
scala_class_path="$var_caller_dir:$var_caller_dir/lib/htsjdk-2.5.0-SNAPSHOT-all.jar:$var_caller_dir/lib/log4j-1.2.17.jar:$var_caller_dir/lib_compile/scala-2.10.2/scala-compiler.jar:$var_caller_dir/lib/commons-math3-3.6.1.jar"
alias scala="/home/nfs/thiesgehrmann/scala-2.11.7/bin/scala -classpath $scala_class_path";
alias scalac="/home/nfs/thiesgehrmann/scala-2.11.7/bin/scalac -classpath $scala_class_path";

scala_cmd="/home/nfs/thiesgehrmann/scala-2.11.7/bin/scala -classpath $scala_class_path"


reference_file="${SCRIPTDIR}/example/genomes/refgenome.fasta";
gff_file="${SCRIPTDIR}/example/genomes/genome.gff"
prebuilt_sjdb="${SCRIPTDIR}/example/empty_splice_junction_db.tsv";
nonempty_sjdb="${SCRIPTDIR}/example/splice_junction_db.tsv";
data_file="$SCRIPTDIR/example/data_file.tsv"
tree_file="$SCRIPTDIR/example/tree.tsv"

#Gene annotation files
gene_coding_regions="${SCRIPTDIR}/example/gene_coding_regions.tsv"
gene_regions="${SCRIPTDIR}/example/gene_regions.tsv"
gene_names="${SCRIPTDIR}/example/gene_names.tsv"

OUTPUT_DIR='/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3/';
OUTPUT_DIR='/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/schco3_unique/'
OUTPUT_DIR="$SCRIPTDIR/example_output"

###############################################################################
###############################################################################
###############################################################################
  # Generate index using the SJDB we already have...

mkdir $OUTPUT_DIR/genome_index

STAR --runMode genomeGenerate \
     --runThreadN 2 \
     --sjdbOverhang 99 \
     --genomeDir $OUTPUT_DIR/genome_index \
     --genomeFastaFiles $reference_file \
     --sjdbFileChrStartEnd $prebuilt_sjdb;

###############################################################################
  # Generate reference index

samtools faidx $reference_file

###############################################################################
###############################################################################
###############################################################################
  # Align reads to the reference

while read -u10 line; do

  R1=`echo $line | cut -d\  -f1`;
  R2=`echo $line | cut -d\  -f2`;
  id=`echo $line | cut -d\  -f3`;

  star_cmd="STAR --genomeDir $OUTPUT_DIR/genome_index \
                 --runThreadN 10 \
                 --alignIntronMax 5000 \
                 --alignIntronMin 10 \
                 --outFileNamePrefix $OUTPUT_DIR/star.${id}. \
                 --outTmpDir /tmp/${star}.${id}.${RANDOM} \
                 --genomeLoad LoadAndRemove \
                 --readFilesIn $R1 $R2 \
                 --outSAMattributes All \
                 --outSAMstrandField intronMotif \
                 --outFilterMultimapNmax 1 \
                 --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                 --outReadsUnmapped Fastx \
                 --outSAMtype BAM SortedByCoordinate \
                 --limitBAMsortRAM 100000000000"

  # If you are using slurm, you can do this
  #echo "sauto short --mem 20000 -cmd \"source /opt/insy/env.el7/paths && $star_cmd\""
  #sauto short --job-name "star.$id" --mem 30000 -cmd "source /opt/insy/env.el7/paths && $star_cmd"
  $star_cmd

daone 10< <(cat $data_file | grep -v '^#') > $OUTPUT_DIR/alignment.stdout


###############################################################################
###############################################################################
###############################################################################
# Mark duplicates

PICARD=/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/picard/picard-tools-2.5.0/picard.jar

while read -u10 line; do
  id=`echo $line | cut -d\  -f3`;
 
  picard_cmd1="java -jar $PICARD AddOrReplaceReadGroups I=$OUTPUT_DIR/star.${id}.Aligned.sortedByCoord.out.bam O=$OUTPUT_DIR/picardRG.${id}.bam RGID=${id} RGLB=PA RGPL=illumina RGSM=${id} RGPU=HISEQ2500"
  picard_cmd2="java -jar $PICARD MarkDuplicates I=$OUTPUT_DIR/picardRG.${id}.bam O=$OUTPUT_DIR/picardMD.${id}.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$OUTPUT_DIR/picardMD.output_metrics.${id};"
  picard_cmd="$picard_cmd1;$picard_cmd2";
  # If you want to use slurm!
  #sauto short --job-name "picard_${id}" --mem 30000 --cpus-per-task 1 -cmd "$picard_cmd2"
  echo $id
  $picard_cmd1
  $picard_cmd2

done 10< <(cat $data_file | grep -v '^#')

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
  # Call variants with my thing!!!

# First make sure it compiles :P
cd $var_caller_dir && scalac Fasta.scala GFF.scala SNP.scala cmpCalls.scala variantCaller.scala
cd $SCRIPTDIR

cd $var_caller_dir && scalac sampleTree.scala GFF.scala VCF.scala assignOrigin.scala filterVCF.scala phaseSNPs.scala  slidingWindowSNPs2.scala
cd $SCRIPTDIR

while read -u10 line; do
  
  id=`echo $line | cut -d\  -f3`;
  output_file="$OUTPUT_DIR/varcalls.binom.${id}.vcf";
  vc_cmd="$scala_cmd -J-Xmx48G variantCaller $prebuilt_sjdb $reference_file $OUTPUT_DIR/picardMD.${id}.bam rnasnp.${id} $output_file";

  #sauto short --job-name "varCall_${id}" --mem 50000 -cmd "$vc_cmd"
  echo $vc_cmd;
  $vc_cmd

done 10< <(cat $data_file | grep -v '^#')

###############################################################################
###############################################################################
  # merge and sort

#First NX filter
while read -u10 line; do
  
  id=`echo $line | cut -d\  -f3`;
  output_file="$OUTPUT_DIR/varcalls.binom.NXfilt.${id}.vcf";
  nx_cmd="grep -v NX $OUTPUT_DIR/varcalls.binom.${id}.vcf > $output_file"
  
  #sauto short --job-name "varCall_${id}" --mem 50000 -cmd "$nx_cmd"
  #echo $vc_cmd;
  grep -v NX $OUTPUT_DIR/varcalls.binom.${id}.vcf > $output_file

done 10< <(cat $data_file | grep -v '^#')

#And THEN merge&sort
files=""
while read -u10 line; do
  id=`echo $line | cut -d\  -f3`;
  files="${files},${OUTPUT_DIR}/varcalls.binom.NXfilt.${id}.vcf"
done 10< <(cat $data_file | grep -v '^#')
outfile="$OUTPUT_DIR/varcalls.binom.sorted.merged.NXfilt.vcf"
ms_cmd="$SCRIPTDIR/merge_sort_vcf.sh $files $outfile --parallel 20 --buffer-size=390G"

#sauto bigmem --time 12:00:00 --job-name "ms_varcall" --mem 400000 -cmd "$ms_cmd"
$ms_cmd

###############################################################################
###############################################################################
###############################################################################
#  Assign origin to each SNP from RNA-SNP
vc_cmp_cmd="$scala_cmd -J-Xmx400G assignOrigin $OUTPUT_DIR/varcalls.binom.sorted.merged.NXfilt.vcf $tree_file  $OUTPUT_DIR/varcall.binom.origins.vcf"
#sauto bigmem --time 12:00:00 --job-name "vc_origin" --mem 400000 -cmd "$vc_cmp_cmd"
$vc_cmp_cmd


  # Count the number of bases with a depth greater than X per sample
for NR in 5 10 15 20 50; do
  ls $OUTPUT_DIR \
    | grep -e 'varcalls.binom.[0-9]\+[.]vcf$' \
    | while read x; do
      c_cmd="$scala_cmd filterVCF $OUTPUT_DIR/$x - infoMin DP ${NR} | grep -v '^#' | wc -l > $OUTPUT_DIR/$x.min${NR}.count"
      echo $c_cmd
      sauto short --job-name "countLOL" --mem 40000 -cmd "$c_cmd"
    done
done

for NR in 5 10 15 20 50; do
  ls $OUTPUT_DIR \
    | grep -e "varcalls.binom.[0-9]\+[.]vcf.min${NR}.count$" \
    | while read x; do
      n=`echo $x | cut -d. -f3`
      c=`cat $OUTPUT_DIR/$x`;
      echo $n,$c
    done > $OUTPUT_DIR/covered_regions.min${NR}.csv
done

###############################################################################
###############################################################################
###############################################################################

function count_info_feature(){
  tr '\t' '\n' | tr ';' '\n' | grep -e "^$1" | cut -d= -f2 | tr ',' '\n' | sort | uniq -c | sed -e 's/^[ ]\+//' | tr ' ' '\t' | sort -k1,1n
}

###############################################################################

  # Get only the SNPs that PASS
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.vcf $OUTPUT_DIR/varcall.binom.origins.pass.vcf pass

  # Annotate these SNPs with mating type loci, genes and named genes
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.pass.vcf - annotate $gene_names GN \
  | scala filterVCF - - annotate $gene_regions GID \
  | scala filterVCF - $OUTPUT_DIR/varcall.binom.origins.annotated.vcf annotate $gene_coding_regions GCR

  # Convert to GFF so we can view it in IGV
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf /dev/null toGFF3 $OUTPUT_DIR/varcall.binom.origins.annotated.gff

  # Produce the Tree for all the genes
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf /dev/null toDOT $tree_file Origin /dev/stdout | dot -Tpdf -o $OUTPUT_DIR/DOT.nonmatingtype.pdf

 # What are these SNPs like?
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf - homoSNP | count_info_feature NN
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf - hetSNP | count_info_feature NN

  # Add the read counts to the mutation counts
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf /dev/null toDOT $tree_file Origin - \
  | grep '^/[*].\+[*]/$' \
  | tr -d '*/' \
  | awk -v "OUTPUT_DIR=$OUTPUT_DIR" \
    '{ filename=OUTPUT_DIR "/star." $1 ".Log.final.out";
       if(system("[ -e " filename " ]") == 0) {
         cmd = "cat " filename " | grep \"Uniquely mapped reads number\" | rev | cut -f1 | rev"
         while (cmd | getline result) {
           print $0 "\t" result
         }
       }
     }' > $OUTPUT_DIR/sample_SNPs_reads.tsv

  # Annotate the deleteriousness of the SNPs
scala filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf $OUTPUT_DIR/varcall.binom.origins.GCR.vcf isolateRegions $gene_coding_regions GCR
python varCaller/filterVCFDeleteriousness.py $OUTPUT_DIR/varcall.binom.origins.GCR.vcf $OUTPUT_DIR/varcall.binom.origins.nomatingtype.GCR.deleteriousness.vcf $gff_file $reference_file

  # The SNPs in domains, how many are deleterious?
cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.GCR.deleteriousness.vcf | scala filterVCF - - infoStringEq DM "" not | count_info_feature DL

cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.GCR.deleteriousness.vcf \
  | scala phaseSNPs - $OUTPUT_DIR/varcall.binom.origins.nomatingtype.GCR.deleteriousness.phased.vcf GID

cat $OUTPUT_DIR/varcall.binom.origins.annotated.vcf \
  | scala filterVCF - - isolateRegions $gene_coding_regions GCR \
  | scala phaseSNPs - - GID

  # Get window for ALL hetero/homozygous SNPs in the wildtype
win_cmd="cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf \
  | grep 70s_wildtype \
  | $scala_cmd slidingWindowSNPs2 - $OUTPUT_DIR/windowSizes.wildtype.bedgraph 10000"
sauto short --job-name "wincmd_hh" --mem 50000 -cmd "$win_cmd"

  # Get window for ONLY heterozygous SNPs in the wildtype
win_cmd="cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf  \
  | grep 70s_wildtype \
  | $scala_cmd filterVCF - - hetSNP \
  | $scala_cmd slidingWindowSNPs2 - $OUTPUT_DIR/windowSizes.wildtype.het.bedgraph 10000"
sauto short --job-name "wincmd_h" --mem 50000 -cmd "$win_cmd"

  # Only heterozygous SNPs NOT in the wildtype
win_cmd="cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf \
  | grep -v  70s_wildtype \
  | $scala_cmd filterVCF - - hetSNP \
  | $scala_cmd slidingWindowSNPs2 - $OUTPUT_DIR/windowSizes.nowildtype.het.bedgraph 10000"
sauto short --job-name "wincmd_ha" --mem 50000 -cmd "$win_cmd"

  # Get the SNPs to estimate mutation rate
cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf \
  | scala filterVCF - /dev/null toDOT $tree_with_wildtype MutationRate - \
  | grep -e '^/[*]' \
  | tr -d '/*' \
  > $OUTPUT_DIR/mutation_rate_samples.tsv


  # Mutation rate only in the mating type loci!
cat $OUTPUT_DIR/varcall.binom.origins.matingtype.vcf \
  | scala filterVCF - /dev/null toDOT $tree_with_wildtype MTLSNPs - \
  | grep -e '^/[*]' \
  | tr -d '/*' \
  > $OUTPUT_DIR/mutation_rate_mtl.tsv

  # Look at mutations in highly mutated regions
cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf \
  | grep -v "70s" \
  | scala slidingWindowSNPs - - 10000 20 \
  | scala filterVCF $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf - isolateRegions - SW \
  | grep -v "70s" \
  | scala filterVCF - /dev/null toDOT $tree_with_wildtype SlidingWindowSNPs -  \
  | dot -Tpdf -o $OUTPUT_DIR/DOT.SlidingWindowSNPs.pdf

cat $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf \
  | grep -v "70s" \
  | scala slidingWindowSNPs - - 10000 20 \
  | scala filterVCF $OUTPUT_DIR/varcall.binom.origins.nomatingtype.vcf - isolateRegions - SW \
  | scala filterVCF - /dev/null toDOT $tree_with_wildtype MTLSNPs - \
  | grep -e '^/[*]' \
  | tr -d '/*' \
  > $OUTPUT_DIR/mutation_rate_slidingWindow.tsv

###############################################################################


