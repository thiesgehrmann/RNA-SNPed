# Generated data

## Simulating accumulated variants

To simulate the accumulation of variants, I generated many different FASTA files hierarchically.

* A parental genome
 * A child genome with 4 variants
  1. A grand child genome with 3 variants
  2. A second grand child genome with 3 variants
 * A second child genome with 4 variants
  1. A third grand child genome with 4 variants
  2. A fourth grand child genome with 3 variants

These FASTA files can be found in the `genomes` directory. 

## Generating RNA-Seq data

The example contains generated RNA-Seq data, using  wgsim (for the example without any sequencing errors).
See `gen_fastq.sh` for information about the generation.
