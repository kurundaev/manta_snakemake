# Run strelka on multicentric GBM
# 27 September 2016
# Julia X.M. Yin


configfile: "/home/julyin/analysis/template_scripts/config.yaml"

# Adjust these variables
SAMPLE_PRIM_LIST = ['DI', 'JG', 'WG']
SAMPLE_RECUR_LIST = ['EW', 'PM', 'RH']

NORMAL = "normal"
PRIMARY = "primary"
RECURRENT = "recurrent"

BAM_DIR = "/share/ScratchGeneral/julyin/cabaret/run/phase2/"
RAW_DIR = "/share/ClusterShare/thingamajigs/julyin/genome_gbm_variant_calling/manta/b37/raw-manta/"
FILT_DIR = "/share/ClusterShare/thingamajigs/julyin/genome_gbm_variant_calling/manta/b37/filtered-manta/"


ONCO_DIR = "/share/ClusterShare/thingamajigs/julyin/genome_gbm_variant_calling/strelka/b37/oncotator-strelka/"
VEP_DIR = "/share/ClusterShare/thingamajigs/julyin/genome_gbm_variant_calling/strelka/b37/vep-strelka/"


rule all:
    input:
        expand(FILT_DIR + "{sample}/{sample}." + NORMAL + ".vs." + PRIMARY + "/filtered_somaticSV.vcf", sample = SAMPLE_PRIM_LIST),
        expand(FILT_DIR + "{sample}/{sample}." + NORMAL + ".vs." + RECURRENT + "/filtered_somaticSV.vcf", sample = SAMPLE_RECUR_LIST)

rule manta_prim:
    input:
        normal_bam = BAM_DIR + "{sample}." + NORMAL + ".dedup.realn.bam",
        tumour_bam = BAM_DIR + "{sample}." + PRIMARY + ".dedup.realn.bam"
    output:
        RAW_DIR + "{sample}/{sample}." + NORMAL + ".vs." + PRIMARY + "/results/variants/somaticSV.vcf"
    params:
        outdir = VCF_DIR + "{sample}/{sample}." + NORMAL + ".vs." + PRIMARY
    message:
        """ [][][] Manta variant calling [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        bash /home/julyin/analysis/template_scripts/manta_snakemake/manta_sv_controlvnormal.sh \
            {input.normal_bam} {input.tumour_bam} {config[refs][reference_human_g1kv37_decoy]} {params.outdir}
        """

rule filter_prim:
    input:
        RAW_DIR + "{sample}/{sample}." + NORMAL + ".vs." + PRIMARY + "/results/variants/somaticSV.vcf"
    output:
        FILT_DIR + "{sample}/{sample}." + NORMAL + ".vs." + PRIMARY + "/filtered_somaticSV.vcf"
    message:
        """ [][][] Filter manta [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        python /home/julyin/analysis/template_scripts/manta_snakemake/pinese_manta_svfilter.py \
            --input {input} \
            --output {output} \
            --normal-n {sample} + NORMAL \
            --tumour-t {sample} + PRIMARY \
            --mintumourvaf 0.3
        """

rule manta_recur:
    input:
        normal_bam = BAM_DIR + "{sample}." + NORMAL + ".dedup.realn.bam",
        tumour_bam = BAM_DIR + "{sample}." + RECURRENT + ".dedup.realn.bam"
    output:
        RAW_DIR + "{sample}/{sample}." + NORMAL + ".vs." + RECURRENT + "/results/variants/somaticSV.vcf"
    params:
        outdir = VCF_DIR + "{sample}/{sample}." + NORMAL + ".vs." + RECURRENT
    message:
        """ [][][] Manta variant calling [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        bash /home/julyin/analysis/template_scripts/manta_snakemake/manta_sv_controlvnormal.sh \
            {input.normal_bam} {input.tumour_bam} {config[refs][reference_human_g1kv37_decoy]} {params.outdir}
        """

rule filter_recur:
    input:
        RAW_DIR + "{sample}/{sample}." + NORMAL + ".vs." + RECURRENT + "/results/variants/somaticSV.vcf"
    output:
        FILT_DIR + "{sample}/{sample}." + NORMAL + ".vs." + RECURRENT + "/filtered_somaticSV.vcf"
    message:
        """ [][][] Filter manta [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        python /home/julyin/analysis/template_scripts/manta_snakemake/pinese_manta_svfilter.py \
            --input {input} \
            --output {output} \
            --normal-n {sample} + NORMAL \
            --tumour-t {sample} + RECURRENT \
            --mintumourvaf 0.3
        """

#============================

rule cull_normal_oncotator:
    input:
        ONCO_DIR + "{sample}/{sample}." + NORMAL + ".vs." + TUMOUR + ".raw.{type}.oncotator.maf"
    output:
        ONCO_DIR + "{sample}/{sample}." + NORMAL + ".vs." + TUMOUR + ".{type}.oncotator.maf"
    message:
        """ [][][] Oncotator annotation [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        grep -v 'NORMAL' {input} > {output}
        """

rule vep_annotation:
    input:
        VCF_DIR + "{sample}/{sample}." + NORMAL + ".vs." + TUMOUR + "/myAnalysis/results/passed.somatic.{type}.vcf"
    output:
        vep = VEP_DIR + "{sample}/{sample}." + NORMAL + ".vs." + TUMOUR + ".{type}.vep.vcf",
        stats = VEP_DIR + "{sample}/{sample}." + NORMAL + ".vs." + TUMOUR + ".{type}.vep.stats.html"
    message:
        """ [][][] VEP annotation [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: 
        """
        module load marcow/perl/5.14.2
        module load gi/vep/76
        module load gi/samtools/0.1.19
        perl /share/ClusterShare/software/contrib/gi/vep/76/variant_effect_predictor.pl \
             --cache \
             --dir /share/ClusterShare/biodata/contrib/gi/vep \
             --port 3337 \
             --offline \
             --input_file {input} \
             --format vcf \
             --output_file {output.vep} \
             --vcf \
             --stats_file {output.stats} \
             --stats_text \
             --force_overwrite \
             --canonical \
             --fork 32 \
             --sift b \
             --polyphen b \
             --symbol \
             --numbers \
             --terms so \
             --biotype \
             --total_length \
             --plugin LoF,human_ancestor_fa:/share/ClusterShare/biodata/contrib/gi/LOFTEE/1.0/human_ancestor.fa.rz --fields Consequenc
        e,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,Feature_type,cDNA_position,CDS_
        position,Existing_variation,DISTANCE,STRAND,CLIN_SIG,LoF_flags,LoF_filter,LoF,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,
        CADD_raw,CADD_phred,Reliability_index
        """


rule vaf_extraction:
    input:
        expand(VEP_DIR + "strelka-vep-" + NORMAL + ".vs.{sample}/{sample}.snvs.vep.vcf", sample = SAMPLE_LIST)
        #g53 = VCF_DIR + "strelka-" + NORMAL + ".vs.G53.primary/myAnalysis/results/passed.somatic.snvs.vcf"
    output:
        "/home/julyin/genome_gbm_variant_calling/vaf/vaf.multicentric.G52.vs.G53.csv"
    message:
        """ [][][] VAF extraction [][][]
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        /opt/perl/bin/perl /home/julyin/analysis/multicentric_G52_G53/strelka/VAF_Generator_modified.pl \
            {input} {output}
            """


