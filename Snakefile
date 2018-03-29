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
__ANNOTATED_OUTDIR__ = "%s/annotated" % __RUN_DIR__
__SUMMARY_OUTDIR__ = "%s/summary" % __RUN_DIR__

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
  conda: "%s/conda.yaml" % __PC_DIR__
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
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    java -Xmx400G -jar "{input.jar}" assignOrigin "{input.vcf}" "{input.tree}" "{output.vcf}"
  """

rule filter_pass:
  input:
    jar = rules.build_rnasnped.output.jar,
    vcf = rules.assign_origin.output.vcf
  output:
    vcf = "%s/varcall.binom.origins.pass.vcf" % __RNASNPED_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    java -jar "{input.jar}" filterVCF "{input.vcf}" "{output.vcf}" pass
  """

###############################################################################
# Annotate the variants we found

rule gene_annots:
  input:
    gff = dconfig["gff_file"],
  output:
    names = "%s/gene_names.tsv" % __ANNOTATED_OUTDIR__,
    regions = "%s/gene_regions.tsv" % __ANNOTATED_OUTDIR__,
    coding_regions = "%s/coding_regions.tsv" % __ANNOTATED_OUTDIR__,
  params:
    nameAttr = dconfig["GFF_name_attribute"],
    geneFeature = dconfig["GFF_gene_feature"],
    cdsFeature  = dconfig["GFF_cds_feature"]
  conda: "%s/conda.yaml" % __PC_DIR__
  run:
    from pipeline_components import biu

    gff = biu.formats.GFF3(input.gff)
    genes = [ e for e in gff.entries if (e.feature == params.geneFeature) ]
    with open(output.names, "w") as ofd:
      for e in [e for e in genes if (params.nameAttr in e.attr) ]:
        ofd.write("%s\t%d\t%d\t%s\n" % (e.seqid, e.start, e.end, e.attr[params.nameAttr]))
      #efor
    #ewith

    with open(output.regions, "w") as ofd:
      for e in genes:
        ofd.write("%s\t%d\t%d\t%s\n" % (e.seqid, e.start, e.end, e.attr["ID"]))
      #efor
    #ewith

    with open(output.coding_regions, "w") as ofd:
      for e in genes:
        geneID = e.attr["ID"]
        cds = gff.getChildren(geneID, feature=params.cdsFeature).entries
        for cr in cds:
          ofd.write("%s\t%d\t%d\t%s\n" % (cr.seqid, cr.start, cr.end, geneID))
        #efor
      #efor
    #ewith

rule annotate_genes:
  input:
    vcf = rules.filter_pass.output.vcf,
    gene_names = rules.gene_annots.output.names,
    gene_regions = rules.gene_annots.output.regions,
    coding_regions = rules.gene_annots.output.coding_regions,
    jar = rules.build_rnasnped.output.jar
  output:
    vcf = "%s/varcall.binom.origins.annotated.vcf" % __ANNOTATED_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    java -jar "{input.jar}" filterVCF "{input.vcf}" - annotate "{input.gene_names}" GN \
      | java -jar "{input.jar}" filterVCF - - annotate "{input.gene_regions}" GID \
      | java -jar "{input.jar}" filterVCF - "{output.vcf}" annotate "{input.coding_regions}" GCR
  """

rule annotate_deleterious_vars:
  input:
    jar = rules.build_rnasnped.output.jar,
    gff = dconfig["gff_file"],
    ref = dconfig["reference_file"],
    vcf = rules.annotate_genes.output.vcf
  output:
    vcf = "%s/varcall.binom.origins.deleteriousness.vcf" % __ANNOTATED_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  params:
    script = '%s/annotate_deleteriousness.py' % __PC_DIR__
  shell: """
    "{params.script}" "{input.vcf}" "{output.vcf}" "{input.gff}" "{input.ref}"
  """

###############################################################################
# Summarize the variants we found

rule dot_output:
  input:
    jar = rules.build_rnasnped.output.jar,
    vcf = rules.annotate_deleterious_vars.output.vcf,
    tree = dconfig["tree_file"]
  output:
    pdf = "%s/snps_tree.pdf" % __SUMMARY_OUTDIR__,
    png = "%s/snps_tree.png" % __SUMMARY_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    java -jar "{input.jar}" filterVCF "{input.vcf}" /dev/null toDOT "{input.tree}" Origin - \
      | dot -Tpdf -o "{output.pdf}"

    java -jar "{input.jar}" filterVCF "{input.vcf}" /dev/null toDOT "{input.tree}" Origin - \
      | dot -Tpng -o "{output.png}"
  """

rule variants_per_node:
  input:
    jar = rules.build_rnasnped.output.jar,
    vcf = rules.annotate_deleterious_vars.output.vcf,
    tree = dconfig["tree_file"]
  output:
    tsv = "%s/variants_per_node.tsv" % __SUMMARY_OUTDIR__
  conda: "%s/conda.yaml" % __PC_DIR__
  shell: """
    echo -e "ID\tName\tUnique\tUniqueYes\tUniqueMaybe\tTotal\tTotalYes\tTotalMaybe" > "{output.tsv}"

    java -jar "{input.jar}" filterVCF "{input.vcf}" /dev/null toDOT "{input.tree}" Counts - \
      | grep -e '^/[*]' \
      | tr -d '/*' \
      >> "{output.tsv}"
  """

rule effect_summary:
  input:
    vcf = rules.annotate_deleterious_vars.output.vcf,
  output:
    tsv = "%s/effect_counts.tsv" % __SUMMARY_OUTDIR__
  shell: """

    function count_info_feature(){{
      tr '\t' '\n' | tr ';' '\n' | grep -e "^$1" | cut -d= -f2 | tr ',' '\n' | sort | uniq -c | sed -e 's/^[ ]\+//' | tr ' ' '\t' | sort -k1,1n
    }}

    echo -e "#S = Synonymous\n#M = Missense\n#N = Nonsense\nEffect\tCount" > "{output.tsv}"

    cat "{input.vcf}" \
      | count_info_feature "DL" \
      | awk -F $'\t' 'BEGIN{{OFS=FS}} {{ print $2, $1 }}' \
      >> "{output.tsv}"

  """

###############################################################################
# Total output

rule all:
  input:
    vcf = rules.annotate_deleterious_vars.output.vcf,
    effect_summary = rules.effect_summary.output.tsv,
    vpn = rules.variants_per_node.output.tsv,
    pdf = rules.dot_output.output.pdf,
    png = rules.dot_output.output.png
    
