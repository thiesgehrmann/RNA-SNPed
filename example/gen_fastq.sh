#!/usr/bin/sh

# GENERATE FASTQ DATA. COMPLETELY NOT RELATED TO ACTUAL SEQUENCING DATA.
#PAIRED END SEQUENCING DATA

###############################################################################

function makepath() {

  output_r1="$1"
  output_r2="$2"
  declare -a path=("${!3}")
  
  #rm $output_r1
  #rm $output_r2

  for node in ${path[@]}; do
    echo $node
    t_r1=`mktemp`;
    t_r2=`mktemp`;

    # I disable all error data, just for the example
    wgsim -e 0 -d 100 -s 25 -N 100 -1 70 -2 70 -r 0 -R 0 -X 0 $node $t_r1 $t_r2
    cat $t_r1 >> $output_r1
    cat $t_r2 >> $output_r2
  done

}


###############################################################################

gc1=(genomes/parent.fasta genomes/child_1.fasta genomes/grand_child_1.1.fasta)
gc2=(genomes/parent.fasta genomes/child_1.fasta genomes/grand_child_1.2.fasta)
gc3=(genomes/parent.fasta genomes/child_2.fasta genomes/grand_child_2.1.fasta)
gc4=(genomes/parent.fasta genomes/child_2.fasta genomes/grand_child_2.2.fasta)

makepath data/sample_0_R1.fastq data/sample_0_R2.fastq gc1[@]
makepath data/sample_1_R1.fastq data/sample_1_R2.fastq gc2[@]
makepath data/sample_2_R1.fastq data/sample_2_R2.fastq gc3[@]
makepath data/sample_3_R1.fastq data/sample_3_R2.fastq gc4[@]
