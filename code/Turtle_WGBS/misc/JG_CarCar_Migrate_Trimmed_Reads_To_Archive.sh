#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:00:00
#$ -l h_vmem=8G
#$ -m bea
#$ -t 61-120

# Created by: James
# Date: 17.06.2022
# Aim: Copy trimmed fastq files from scratch to archive

# copying from the Raw_Data dir
DATA_DIR=/data/scratch/btx355/Chris_Grant_WGBS/Raw_Data
# to the archive dir
ARCH_DIR=/data/archive/archive-SBCS-EizaguirreLab/Turtle_WGBS/soapnuke/trimmed_cutadapt

# Samples listed in this file
Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/FullSampleList.txt

# sample to process
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)

# make the dir in the archive to receive trimmed files
mkdir -p ${ARCH_DIR}/${SAMPLE}

cd ${DATA_DIR}/${SAMPLE}

# remove original fq files
# rm ./*_1.fq.gz
# rm ./*_2.fq.gz

# remove merged fq files
# rm ./*trim_merged.fq.gz 

# copy all _trim.fq.gz files to arch dir
cp ${DATA_DIR}/${SAMPLE}/*_trim.fq.gz ${ARCH_DIR}/${SAMPLE}

# also copy bam file to 
# cp ${DATA_DIR}/${SAMPLE}/*.bam /data/scratch/btx355/Chris_Grant_WGBS/BSB_Aligned_bams

# rm  ${DATA_DIR}/${SAMPLE}/*.bam

# After checking the copying was successful, delete trimmed files from scratch
# Also delete the merged fastq files (better to keep like original, can remake easily)