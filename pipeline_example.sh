#!/usr/bin/bash
set -e

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
###############################################################################
#  CONFIGURATION

PICARD=/home/nfs/thiesgehrmann/groups/w/phd/tasks/somatic_variation/picard/picard-tools-2.5.0/picard.jar
scala="/home/nfs/thiesgehrmann/scala-2.11.7/bin/scala";
scalac="/home/nfs/thiesgehrmann/scala-2.11.7/bin/scalac";
scalalib="/home/nfs/thiesgehrmann/scala-2.11.7/lib/scala-library.jar"

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

OUTPUT_DIR="$SCRIPTDIR/example_output"

###############################################################################
###############################################################################
###############################################################################

function msg {

  echo "###############################################################################"
  echo "# $@"
  echo "###############################################################################"

}

###############################################################################
###############################################################################
###############################################################################
  # Generate index using the SJDB we already have...

mkdir -p $OUTPUT_DIR/genome_index

msg "Generating STAR index"

STAR --runMode genomeGenerate \
     --runThreadN 2 \
     --sjdbOverhang 99 \
     --genomeDir $OUTPUT_DIR/genome_index \
     --genomeFastaFiles $reference_file \
     --sjdbFileChrStartEnd $prebuilt_sjdb;

###############################################################################
  # Generate reference index

msg "Generating Fasta index"
samtools faidx $reference_file

###############################################################################
###############################################################################
###############################################################################
  # Align reads to the reference

msg "Aligning reads"
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

done 10< <(cat $data_file | grep -v '^#') 
#> $OUTPUT_DIR/alignment.stdout


###############################################################################
###############################################################################
###############################################################################
# Mark duplicates

msg "Marking duplicates with PICARD"

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

msg "Call variants using RNASNPed"
while read -u10 line; do
  
  id=`echo $line | cut -d\  -f3`;
  output_file="$OUTPUT_DIR/varcalls.binom.${id}.vcf";
  vc_cmd="java -Xmx48G -jar rnasnped/rnasnped.jar variantCaller $prebuilt_sjdb $reference_file $OUTPUT_DIR/picardMD.${id}.bam rnasnp.${id} $output_file";

  #sauto short --job-name "varCall_${id}" --mem 50000 -cmd "$vc_cmd"
  echo $vc_cmd;
  $vc_cmd

done 10< <(cat $data_file | grep -v '^#')

###############################################################################
###############################################################################
  # merge and sort

msg "Merge and sort VCF files"

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

msg "Assign origin to each SNP"

#  Assign origin to each SNP from RNA-SNP
vc_cmp_cmd="java -Xmx400G -jar rnasnped/rnasnped.jar assignOrigin $OUTPUT_DIR/varcalls.binom.sorted.merged.NXfilt.vcf $tree_file  $OUTPUT_DIR/varcall.binom.origins.vcf"
#sauto bigmem --time 12:00:00 --job-name "vc_origin" --mem 400000 -cmd "$vc_cmp_cmd"
$vc_cmp_cmd


msg "Count base depths, you can remove this if you want..."

  # Count the number of bases with a depth greater than X per sample
for NR in 5 10 15 20 50; do
  ls $OUTPUT_DIR \
    | grep -e 'varcalls.binom.[0-9]\+[.]vcf$' \
    | while read x; do
      c_cmd="java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/$x - infoMin DP ${NR} | grep -v '^#' | wc -l > $OUTPUT_DIR/$x.min${NR}.count"
      $c_cmd
      #sauto short --job-name "countLOL" --mem 40000 -cmd "$c_cmd"
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

msg "Performing filtering operations on the VCF file"

  # Get only the SNPs that PASS
msg "Filter VCF by PASS"
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.vcf $OUTPUT_DIR/varcall.binom.origins.pass.vcf pass

  # Annotate these SNPs with mating type loci, genes and named genes
msg "Annotate SNPs"
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.pass.vcf - annotate $gene_names GN \
  | java -jar rnasnped/rnasnped.jar filterVCF - - annotate $gene_regions GID \
  | java -jar rnasnped/rnasnped.jar filterVCF - $OUTPUT_DIR/varcall.binom.origins.annotated.vcf annotate $gene_coding_regions GCR

  # Convert to GFF so we can view it in IGV
msg "Convert VCF -> GFF3 for IGV view"
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf /dev/null toGFF3 $OUTPUT_DIR/varcall.binom.origins.annotated.gff

  # Produce the Tree for all the genes
msg "Produce tree file containing SNP counts"
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf /dev/null toDOT $tree_file Origin - \
 | dot -Tpdf -o $OUTPUT_DIR/DOT.SNPs.pdf

 # What are these SNPs like?
msg "How many homozygous SNPs are there at each node in the sample tree?"
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf - homoSNP \
 | count_info_feature NN

msg "How many heterozygous SNPs are there at each node in the sample tree?"
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf - hetSNP \
 | count_info_feature NN

msg "Produce tsv file containing SNPs and read depth per node in sample tree"
  # Add the read counts to the mutation counts
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf /dev/null toDOT $tree_file Origin - \
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
msg "Determine deleteriousness of SNPs in coding regions"
java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf $OUTPUT_DIR/varcall.binom.origins.GCR.vcf isolateRegions $gene_coding_regions GCR
python rnasnped/filterVCFDeleteriousness.py $OUTPUT_DIR/varcall.binom.origins.GCR.vcf $OUTPUT_DIR/varcall.binom.origins.GCR.deleteriousness.vcf $gff_file $reference_file

msg "Phase deleterious SNPs"
cat $OUTPUT_DIR/varcall.binom.origins.GCR.deleteriousness.vcf \
  | java -jar rnasnped/rnasnped.jar phaseSNPs - $OUTPUT_DIR/varcall.binom.origins.GCR.deleteriousness.phased.vcf GID

msg "Phase all SNPs using VAFs"
cat $OUTPUT_DIR/varcall.binom.origins.annotated.vcf \
  | java -jar rnasnped/rnasnped.jar filterVCF - - isolateRegions $gene_coding_regions GCR \
  | java -jar rnasnped/rnasnped.jar phaseSNPs - $OUTPUT_DIR/varcall.binom.origins.annotated.phased.vcf  GID

  # Get window for ALL hetero/homozygous SNPs in the wildtype
msg "Calculate sliding window using all SNPs"
cat $OUTPUT_DIR/varcall.binom.origins.annotated.vcf \
  | java -jar rnasnped/rnasnped.jar slidingWindowSNPs2 - $OUTPUT_DIR/windowSizes.wildtype.bedgraph 10000

  # Get window for ONLY heterozygous SNPs in the wildtype
msg "Calculate sliding windows using only heterozygous SNPs"
cat $OUTPUT_DIR/varcall.binom.origins.annotated.vcf  \
  | java -jar rnasnped/rnasnped.jar filterVCF - - hetSNP \
  | java -jar rnasnped/rnasnped.jar slidingWindowSNPs2 - $OUTPUT_DIR/windowSizes.wildtype.het.bedgraph 10000

  # Only heterozygous SNPs NOT in the wildtype
msg "Calculate sliding window using only homozygous SNPs"
cat $OUTPUT_DIR/varcall.binom.origins.annotated.vcf \
  | java -jar rnasnped/rnasnped.jar filterVCF - - hetSNP \
  | java -jar rnasnped/rnasnped.jar slidingWindowSNPs2 - $OUTPUT_DIR/windowSizes.nowildtype.het.bedgraph 10000

  # Get the number of SNPs to estimate mutation rate
msg "Create tsv containing the number of SNPs in each node in the tree"
cat $OUTPUT_DIR/varcall.binom.origins.annotated.vcf \
  | java -jar rnasnped/rnasnped.jar filterVCF - /dev/null toDOT $tree_file MutationRate - \
  | grep -e '^/[*]' \
  | tr -d '/*' \
  > $OUTPUT_DIR/mutation_rate_samples.tsv

msg "Draw tree using only SNPs in highly mutated regions"
cat $OUTPUT_DIR/varcall.binom.origins.annotated.vcf \
  | java -jar rnasnped/rnasnped.jar slidingWindowSNPs2 - - 1000 \
  | grep -v -e "^#" -e "^track" \
  | awk '{ if( $4 > 16){print $0}}' \
  | java -jar rnasnped/rnasnped.jar filterVCF $OUTPUT_DIR/varcall.binom.origins.annotated.vcf - isolateRegions - SW \
  | java -jar rnasnped/rnasnped.jar filterVCF - /dev/null toDOT $tree_file SlidingWindowSNPs -  \
  | dot -Tpdf -o $OUTPUT_DIR/DOT.SlidingWindowSNPs.pdf


###############################################################################


