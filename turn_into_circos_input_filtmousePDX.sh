#!/bin/bash

#VCF=/home/julyin/genome_gbm_variant_calling/manta/xenograft_G89/G89control.vs.G89HFMTissueDNA.manta.filtered.vcf
VCF=/home/julyin/genome_gbm_variant_calling/manta/xenograft_G89/G89control.vs.G89PDX_mousefilt.manta.filtered.vcf

DIR="/home/julyin/genome_gbm_variant_calling/manta/xenograft_G89/circos/control_vs_PDXmousefilt"

OUT_PARSER=$DIR/G89control.vs.G89PDX_mousefilt.manta.circos.raw.txt



mkdir -p $DIR

echo "Parse filtered manta file to show mutation types via XQ script"
echo "$VCF"
date

python quek_manta_parser.py $VCF > $OUT_PARSER



echo "Convert to circos format"

FILE="$DIR/chrM_rm.tmp"

grep -v hsM $OUT_PARSER > $FILE

echo "Number of insertions:" 
grep -c "INS" $FILE
grep -e INS $FILE | replace "INS" "fill_color=blue" > $DIR/ins.tmp


echo "Number of deletions:" 
grep -c "DEL" $FILE
grep -e DEL $FILE | replace "DEL" "fill_color=red" > $DIR/del.tmp

echo "Number of duplications:" 
grep -c "DUP" $FILE
grep -e DUP $FILE | replace "DUP" "fill_color=purple" > $DIR/dup.tmp

cat $DIR/ins.tmp $DIR/del.tmp $DIR/dup.tmp > $DIR/indel_dups.tmp
awk -F"\t" '{print $1 "\t" $2 "\t" $4 "\t" $5}' $DIR/indel_dups.tmp > $DIR/circos_indel_dups.txt


#=====
echo "Number of translocations:" 
grep -c "CTX" $FILE
grep -e CTX $FILE | replace "CTX" "color=green" > $DIR/ctx.tmp

echo "Number of inversions:" 
grep -c "INV" $FILE
grep -e INV $FILE | replace "INV" "color=orange" > $DIR/inv.tmp

cat $DIR/ctx.tmp $DIR/inv.tmp > $DIR/ctx_inv.tmp
awk -F"\t" '{print $1 "\t" $2  "\t" $2  "\t" $3  "\t" $4  "\t" $4 "\t" $5}' $DIR/ctx_inv.tmp > $DIR/circos_ctx_inv.txt

echo "Finished!"
date

#replace "DEL" "color=red" -- manta_raw.txt
#replace "INS" "color=green" -- manta_raw.txt
#replace "INV" "color=yellow" -- manta_raw.txt
#replace "DUP" "color=purple" -- manta_raw.txt
#replace "CTX" "color=blue" -- manta_raw.txt


