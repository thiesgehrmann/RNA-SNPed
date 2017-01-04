# RNA-SNPed
RNA-SNPed is a method to detect SNPs in an Experimental Design from RNA-Seq data

## Dependencies
 The core tool, rnasnped, is standalone, but the full pipeline, including downstream analysis has several dependencies

 * Scala
 * Java
 * Picard
 * STAR RNA-Seq aligner
 * Samtools
 * Python
 * Ibidas
 * Pandas
 * Matplotlib
 * Numpy
 * Scipy
 * Seaborn

## Usage

To run the pipeline with the example dataset, simply execute `./pipeline_example.sh`.

## VarCaller tool

The core of the pipeline is the varCaller tool, written in scala.
It contains four different sub-tools:
 * variantCaller, which calls SNPs from aligned BAM files,
 * assignOrigin, which determins the origin of a SNP in a sample hierarchy,
 * filterVCF, which performs various filtering operations on VCF files, and
 * slidingWindowSNPs2, to count SNPs in a sliding window

### Building varCaller

The jar file is already provided in rnasnped/rnasnped.jar, but if you want to compile it yourself

To compile the tool, run the following command:
```shell
  cd src && ./build.sh && cd ..
```
