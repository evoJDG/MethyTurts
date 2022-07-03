#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -t 2-59


# Created by: Charley (adapted from James)
# Date: 03.07.2022
# Aim: To prepare BSBolt output for MethylKit by splitting methylation calls
# by C genomic context (CG, CHG, CHH, unmap)


### Prep environment ###

module load python # v.3.8.5
source /data/SBCS-EizaguirreLab/James/environments/my_bsbolt/bin/activate # v.1.5.0


### Extract sample IDs ###

Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/FullSampleList_no169.txt

SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)


### Set directories ###

METH_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/05_Methylation_Calling/meth_calls/$SAMPLE
OUT_DIR=$METH_DIR/out

mkdir -p $OUT_DIR
cd $METH_DIR


### Separate by C genomic context ###

INPUT_FILE=`ls $METH_DIR | grep *CGmap.gz`

head $INPUT_FILE | gzip -f -c -d $INPUT_FILE > $OUT_DIR/test.txt

gzip -f -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4=="CG")  { printf("%s\n", $0); } }' | gzip -c --best > $METH_DIR/${SAMPLE}.CG.map.gz)
gzip -f -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4=="CHG")  { printf("%s\n", $0); } }' | gzip -c --best > $METH_DIR/${SAMPLE}.CHG.map.gz)
gzip -f -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4=="CHH")  { printf("%s\n", $0); } }' | gzip -c --best > $METH_DIR/${SAMPLE}.CHH.map.gz)
gzip -f -c -d $INPUT_FILE | (awk -F $'\t' '{ if ($4!="CG" && $4!="CHG" && $4!="CHH")  { printf("%s\n", $0); } }' | gzip -c --best > $METH_DIR/${SAMPLE}.unmap.gz)

