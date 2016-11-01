#!/bin/bash

#The configuration step creates a new workflow run script in the requested run directory.

. /etc/profile.d/modules.sh
module load kevyin/python/2.7.5

MANTA_INSTALL_PATH="/share/ClusterShare/software/contrib/julyin/manta-0.29.3/manta-0.29.3.centos5_x86_64"

#REF_FILE="/share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.5/hg19/ucsc.hg19.fasta"

NORMAL=$1
TUMOUR=$2
REF_FILE=$3
MANTA_ANALYSIS_PATH=$4

#MANTA_ANALYSIS_PATH="/share/ClusterScratch/julyin/manta_patienttumour"

${MANTA_INSTALL_PATH}/bin/configManta.py \
    --normalBam $NORMAL \
    --tumorBam $TUMOUR \
    --referenceFasta $REF_FILE \
    --runDir ${MANTA_ANALYSIS_PATH}


awk 'NR==1 {$0="#!/share/ClusterShare/software/contrib/kevyin/python/2.7.5/bin/python"}1' ${MANTA_ANALYSIS_PATH}/runWorkflow.py > ${MANTA_ANALYSIS_PATH}/runWorkflow_1.py
chmod 700 ${MANTA_ANALYSIS_PATH}/runWorkflow_1.py


echo "[][][] Start runWorkflow for $NORMAL vs $TUMOUR"
date

${MANTA_ANALYSIS_PATH}/runWorkflow_1.py -m local -j 8

echo "[][][] Finished runWorkflow for $NORMAL vs $TUMOUR"
date

