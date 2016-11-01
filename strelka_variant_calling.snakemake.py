# Run oncotator on G89 patient and PDX
# 27 September 2016
# Julia X.M. Yin


configfile: "./config.yaml"

SAMPLE_LIST = ["G52.primary", "G53TissueDNA"]

# Adjust these variables
NORMAL = "G53controlDNA"
#G52 = "G52.primary"
#G53 = "G53.primary"
TUMOUR = "primary"
BAM_DIR = "/share/ScratchGeneral/julyin/cabaret/run/phase2/"
VCF_DIR = "/home/julyin/genome_gbm_variant_calling/mutect2/"
ANNOTATED_DIR = "/home/julyin/genome_gbm_variant_calling/mutect2_oncotator/"



DIR = "/home/julyin/genome_gbm_variant_calling/strelka-snvs-indels/"
SAMPLE_PDX = "strelka-G89controlDNA.dedup.realn.bam-vs-G89PDX.dedup.realn.bam"
SAMPLE_PATIENT = "strelka-G89controlDNA-vs-G89HFMTissueDNA"
OUT_DIR = "/home/julyin/genome_gbm_variant_calling/strelka_oncotator/" 

#DB_DIR=/share/ClusterShare/software/contrib/julyin/oncotator/oncotator_v1_ds_June112014
VAR_TYPE = ["snvs", "indels"]


SAM_PDX_DIR = DIR + SAMPLE_PDX + "/myAnalysis/results/"
SAM_PATIENT_DIR = DIR + SAMPLE_PATIENT + "/myAnalysis/results/"


rule all:
    input: 
        expand(OUT_DIR + "headerless.G89control.vs.G89PDX.passed.{type}.oncotator.maf", type=VAR_TYPE),
        expand(OUT_DIR + "headerless.G89control.vs.G89HFMTissue.passed.{type}.oncotator.maf", type=VAR_TYPE)


rule run_strelka:
    input:
        normal =
        g52 =
        g53 =
    output:
    message:
        """ [][][] Strelka variant calling of G52 and G53
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        mkdir $WORK_DIR
        cd $WORK_DIR
        cp $STRELKA_INSTALL_DIR/etc/strelka_config_bwa_default.ini strelka_config_bwa_default.ini
        $STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl \
            --normal=$NORMAL_PATH \
            --tumor=$TUMOUR_PATH \
            --ref=$REF_FILE \
            --config=strelka_config_bwa_default.ini \
            --output-dir=./myAnalysis
        cd ./myAnalysis
        make -j 8
        """






rule oncotator_PDX:
    input:
        SAM_PDX_DIR + "passed.somatic.{type}.vcf"
    output:
        OUT_DIR + "G89control.vs.G89PDX.passed.{type}.oncotator.maf"
    message:
        """ [][][] Oncotator annotation of G89control and PDX
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
    oncotator -i VCF --db-dir={config[db_path]} --canonical-tx-file={config[canonical_tx_file]} --output_format=TCGAMAF {input} {output} hg19
    """

rule oncotator_Patient:
    input:
        SAM_PATIENT_DIR + "passed.somatic.{type}.vcf"
    output:
        OUT_DIR + "G89control.vs.G89HFMTissue.passed.{type}.oncotator.maf"
    message:
        """ [][][] Oncotator annotation of G89control and patient
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
    oncotator -i VCF --db-dir={config[db_path]} --canonical-tx-file={config[canonical_tx_file]} --output_format=TCGAMAF {input} {output} hg19
    """


rule remove_comment_lines:
    input:
        OUT_DIR + "G89control.vs.{tumour}.passed.{type}.oncotator.maf"
    output:
        OUT_DIR + "headerless.G89control.vs.{tumour}.passed.{type}.oncotator.maf"
    message:
        """ [][][] Remove headers in maf files
        INPUTS:
        {input}
        OUTPUTS:
        {output}
        """
    shell: """
        sed '/^\#/d' {input} > {output} """
        











