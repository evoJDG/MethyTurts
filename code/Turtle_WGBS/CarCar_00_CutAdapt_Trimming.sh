#!/bin/bash
#$ -pe smp 5
#$ -l h_vmem=2G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -t 1-57

# Created by: Charley
# Date: 13.06.2022
# Aim: Trim fastq files for adapters and quality using cutadapt
# BGI have already trimmed, but repeating to be extra thorough


### Prep environment ### 

module load trimgalore # v.0.6.4 -> contains cutadapt v.2.10


### Extract sample IDs ###

Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/FullSampleList.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)


### Run cutadapt array ### 

# Adapter sequences provided by Yuxin Chen from BGI
# Quality filters same as used by Alice/BGI/literature
# Loop through each fastq pair within directory of each sample ID

for i in $SAMPLE/*_1.fq.gz

do
 echo 'Started: ' "${i:0:-8}"
 
 cutadapt \
 -e 0.1 -q 20 -m 20 -O 1 \
 -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
 -A AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
 -j ${NSLOTS} \
 -o "${i:0:-8}_1_trim.fq.gz" \
 -p "${i:0:-8}_2_trim.fq.gz" \
 $i "${i:0:-8}_2.fq.gz"
 
 echo 'Finished: ' "${i:0:-8}"
done
