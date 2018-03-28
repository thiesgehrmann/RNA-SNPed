import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__

###############################################################################

import json
dconfig = json.load(open("%s/defaults.json"% __PC_DIR__, "r"))
dconfig.update(config)

###############################################################################


__RUN_DIR__ = os.path.abspath(dconfig["outdir"])
__STAR_OUTDIR__ = "%s/star_align" % __RUN_DIR__
__PICARD_OUTDIR__ = "%s/picard" % __RUN_DIR__
__RNASNPED_OUTDIR__ = "%s/rnasnped" % __RUN_DIR__

###############################################################################
# Use STAR to align to reference genome

rule star_index:
  input:
    genomes = dconfig["reference_file"]
  output:
    index = "%s/star.idx" % __STAR_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  params:
    star_index_params = dconfig['star_index_params']
  shell: """
    mkdir {output.index}
    STAR --runMode genomeGenerate {params.star_index_params} --genomeDir {output.index} --genomeFastaFiles {input.genomes}
  """

rule star_pass_one_sample:
  input:
    r1 = lambda wildcards: dconfig["samples"][wildcards.sample]["R1"],
    r2 = lambda wildcards: dconfig["samples"][wildcards.sample]["R2"],
    index = rules.star_index.output.index
  output:
    sj = "%s/aln1.{sample}.SJ.out.tab" % __STAR_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  threads: 5
  params:
    rule_outdir = __STAR_OUTDIR__,
    star_params = dconfig["star_params"]
  shell: """
    STAR --runMode alignReads --runThreadN {threads} \
          {params.star_params} \
          --genomeDir {input.index} \
          --readFilesIn "{input.r1}" "{input.r2}" \
          --outFileNamePrefix "{params.rule_outdir}/aln1.{wildcards.sample}." \
          --outSAMmode None
  """
    

rule star_pass_one:
  input:
    sj = expand("%s/aln1.{sample}.SJ.out.tab" % __STAR_OUTDIR__, sample=dconfig["samples"].keys())
  output:
    sj = "%s/star_sj.tsv" % __STAR_OUTDIR__
  shell: """
    cut -f1,2,3 {input.sj} > "{output.sj}"
  """
 
rule star_pass_two_sample:
  input:
    index = rules.star_index.output.index,
    sj    = rules.star_pass_one.input.sj,
    r1 = lambda wildcards: dconfig["samples"][wildcards.sample]["R1"],
    r2 = lambda wildcards: dconfig["samples"][wildcards.sample]["R2"],
  output:
    stats = "%s/aln2.{sample}.Log.final.out" % __STAR_OUTDIR__,
    sam   = "%s/aln2.{sample}.Aligned.sortedByCoord.out.bam" % __STAR_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  threads: 5
  params:
    rule_outdir = __STAR_OUTDIR__,
    star_params = dconfig["star_params"],
    starIntronMotif = "" if dconfig["strandSpecific"] else "--outSAMstrandField intronMotif"
  shell: """
    STAR --runMode alignReads \
         --runThreadN {threads} \
         {params.star_params} {params.starIntronMotif} \
         --genomeDir {input.index} \
         --readFilesIn "{input.r1}" "{input.r2}" \
         --outFileNamePrefix "{params.rule_outdir}/aln2.{wildcards.sample}." \
         --outSAMtype BAM SortedByCoordinate \
         --sjdbFileChrStartEnd {input.sj}
  """

rule all_star:
  input:
    aln = expand("%s/aln2.{sample}.Aligned.sortedByCoord.out.bam" % __STAR_OUTDIR__, sample=dconfig["samples"].keys())

###############################################################################
# Add duplicate read group labels with PICARD

rule picard_readgroups:
  input:
    bam = lambda wildcards: "%s/aln2.%s.Aligned.sortedByCoord.out.bam" % (__STAR_OUTDIR__, wildcards.sample)
  output:
    bam = "%s/pixardRG.{sample}.bam" % __PICARD_OUTDIR__
  params:
    id = lambda wildcards : dconfig["samples"][wildcards.sample]["ID"]
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    picard AddOrReplaceReadGroups I="{input.bam}" O="{output.bam}" RGID={params.id} RGLB=PA RGPL=illumina RGSM={params.id} RGPU=HISEQ2500
  """

rule picard_duplicates:
  input:
    bam = lambda wildcards: "%s/pixardRG.%s.bam" % (__PICARD_OUTDIR__, wildcards.sample)
  output:
    bam = "%s/pixardMD.{sample}.bam" % __PICARD_OUTDIR__
  params:
    rule_dir = __PICARD_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    picard MarkDuplicates I="{input.bam}" O="{output.bam}" VALIDATION_STRINGENCY=SILENT M="{params.rule_dir}/picardMD.output_metrics.{wildcards.sample}"
  """

rule picard_index:
  input:
    bam = lambda wildcards: "%s/pixardMD.%s.bam" % (__PICARD_OUTDIR__, wildcards.sample)
  output:
    bai = "%s/pixardMD.{sample}.bam.bai" % __PICARD_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    samtools index "{input.bam}" "{output.bai}"
  """

rule all_picard:
  input:
    bam = expand("%s/pixardMD.{sample}.bam" % __PICARD_OUTDIR__, sample=dconfig["samples"].keys()),
    bai = expand("%s/pixardMD.{sample}.bam.bai" % __PICARD_OUTDIR__, sample=dconfig["samples"].keys())

###############################################################################
# Call variants

rule build_rnasnped:
  output:
    jar = "%s/rnasnped/rnasnped.jar" % __INSTALL_DIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    cd `dirname "{output.jar}"` && ./build.sh
  """

rule call_vars:
  input:
    ref = os.path.abspath(dconfig["reference_file"]),
    bam = lambda wildcards: "%s/pixardMD.%s.bam" % (__PICARD_OUTDIR__, wildcards.sample),
    bai = lambda wildcards: "%s/pixardMD.%s.bam.bai" % (__PICARD_OUTDIR__, wildcards.sample),
    sj = rules.star_pass_one.output.sj,
    jar = rules.build_rnasnped.output.jar
  output:
    vcf = "%s/varcalls.binom.{sample}.vcf" % __RNASNPED_OUTDIR__
  params:
    id = lambda wildcards: dconfig["samples"][wildcards.sample]["ID"]
  shell : """
    java -Xmx48G -jar "{input.jar}" variantCaller "{input.sj}" "{input.ref}" "{input.bam}" "rnasnp.{params.id}" "{output.vcf}"
  """

rule all_call_vars:
  input:
    vcf = expand("%s/varcalls.binom.{sample}.vcf" % __RNASNPED_OUTDIR__, sample=dconfig["samples"].keys())

rule filterNX:
  input:
    vcf = lambda wildcards: "%s/varcalls.binom.%s.vcf" % (__RNASNPED_OUTDIR__, wildcards.sample)
  output:
    vcf = "%s/varcalls.NXfilt.{sample}.vcf" % __RNASNPED_OUTDIR__
  shell: """
    grep -v NX "{input.vcf}" > "{output.vcf}"
  """

rule mergesort:
  input:
    vcf = expand("%s/varcalls.NXfilt.{sample}.vcf" % __RNASNPED_OUTDIR__, sample=dconfig["samples"].keys())
  output:
    vcf = "%s/varcalls.NXfilt.vcf" % __RNASNPED_OUTDIR__
  threads : 20
  params:
    input_vcf = ','.join(expand("%s/varcalls.NXfilt.{sample}.vcf" % __RNASNPED_OUTDIR__, sample=dconfig["samples"].keys())),
    mergesort = "%s/pipeline_components/merge_sort_vcf.sh" % __INSTALL_DIR__
  shell: """
    "{params.mergesort}" {params.input_vcf} "{output.vcf}" --parallel {threads} --buffer-size=390G
  """

rule assign_origin:
  input:
    jar = rules.build_rnasnped.output.jar,
    tree = dconfig["tree_file"],
    vcf = rules.mergesort.output.vcf
  output:
    vcf = "%s/varcall.binom.origins.vcf" % __RNASNPED_OUTDIR__
  shell: """
    java -Xmx400G -jar "{input.jar}" assignOrigin "{input.vcf}" "{input.tree}" "{output.vcf}"
  """

rule filter_pass:
  input:
    jar = rules.build_rnasnped.output.jar,
    vcf = rules.assign_origin.output.vcf
  output:
    vcf = "%s/varcall.binom.origins.pass.vcf" % __RNASNPED_OUTDIR__
  shell: """
    java -jar "{input.jar}" filterVCF "{input.vcf}" "{output.vcf}" pass
  """

###############################################################################
# Annotate the variants we found


###############################################################################
# Summarize the variants we found
