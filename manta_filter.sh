#!/bin/bash


IN_PAT="/home/julyin/genome_gbm_variant_calling/manta/xenograft_G89/patient_tumour/results/variants/somaticSV.vcf"
OUT_PAT="/home/julyin/genome_gbm_variant_calling/manta/xenograft_G89/G89control.vs.G89HFMTissueDNA.manta.filtered.vcf"
NORMAL="G89controlDNA"
TUMOUR_PAT="G89HFMTissueDNA"

echo "[][] Running filter for $IN_PAT"
date

python pinese_manta_svfilter.py \
    --input $IN_PAT \
    --output $OUT_PAT \
    --normal $NORMAL \
    --tumour $TUMOUR_PAT \
    --mintumourvaf 0.3 \
    --mintumourdepth 15 \
    --mintumourdepth 15




IN_PDX="/home/julyin/genome_gbm_variant_calling/manta/xenograft_G89/filtered_mouse_PDX/results/variants/somaticSV.vcf"
TUMOUR_PDX="G89PDX_mousefilt"
OUT_PDX="/home/julyin/genome_gbm_variant_calling/manta/xenograft_G89/G89control.vs.G89PDX_mousefilt.manta.filtered.vcf"

echo "[][] Running filter for $IN_PDX"
date

python pinese_manta_svfilter.py \
    --input $IN_PDX \
    --output $OUT_PDX \
    --normal $NORMAL \
    --tumour $TUMOUR_PDX \
    --mintumourvaf 0.3 \
    --mintumourdepth 15 \
    --mintumourdepth 15


sed /^##/d $OUT_PAT > $OUT_PAT.headerless
sed /^##/d $OUT_PDX > $OUT_PDX.headerless



