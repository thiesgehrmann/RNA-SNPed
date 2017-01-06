# RNA-SNPed
RNA-SNPed is a method to detect SNPs in an Experimental Design from RNA-Seq data.

## Dependencies
 The core tool, **rnasnped**, only requires Java, but the full pipeline, including downstream analysis has several dependencies
 * To run rnasnped:
   * Java
 * To compile rnasnped:
   * Scala
 * To run the pipeline:
   * Java
   * Picard
   * STAR RNA-Seq aligner
   * Samtools
   * To asses SNP deleteriousness
     * Python
     * Ibidas

## Usage

To run the pipeline with the example dataset, simply execute `./pipeline_example.sh`.

## **rnasnped** tool

The core of the pipeline is the varCaller tool, written in scala.
It contains four different sub-tools:
 * variantCaller, which calls SNPs from aligned BAM files,
 * assignOrigin, which determins the origin of a SNP in a sample hierarchy,
 * filterVCF, which performs various filtering operations on VCF files, and
 * slidingWindowSNPs2, to count SNPs in a sliding window
 * phaseSNPs, to phase SNPs within a gene based on the VAF

### Building rnasnped

The jar file is already provided in rnasnped/rnasnped.jar, but if you want to compile it yourself

To compile the tool, run the following command:
```shell
  cd rnasnped && ./build.sh && cd ..
```

###Usage

You can run the rnasnped tool using java
```bash
  java -jar rnasnped/rnasnped.jar <task>
```

where <task> takes the value of one of the five different tools listed above.
Additional help is available by adding help as the first parameter to the tool. e.g. `java -jar variantCaller help`

## Output VCF info fields

