#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=24:00:00
#$ -l h_vmem=1G
#$ -m bea

# Created by: Charley
# Date: 16.06.2022
# Aim: Copy trimmed fastq files from scratch to archive
# NB used *-* in order to select hatchling files only, adjust as needed

DATA_DIR=/data/scratch/btx902/Turtle_WGBS/fastqs
ARCH_DIR=/data/archive/archive-SBCS-EizaguirreLab/Turtle_WGBS/soapnuke/trimmed_cutadapt

rsync -alvrhz \
$DATA_DIR/*-* $ARCH_DIR

echo "All Done"

# After checking the copying was successful, delete trimmed files from scratch
# Also delete the merged fastq files (better to keep like original, can remake easily)

