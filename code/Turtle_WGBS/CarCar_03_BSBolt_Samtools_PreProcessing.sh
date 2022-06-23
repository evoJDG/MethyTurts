#!/bin/bash
#$ -l h_vmem=2G
#$ -pe smp 10 
#$ -l h_rt=240:00:00
#$ -o /data/SBCS-EizaguirreLab/Turtle_WGBS/04_Samtools_Pre_Processing/log_files/samtools_adults_pre_processing.stdout
#$ -e /data/SBCS-EizaguirreLab/Turtle_WGBS/04_Samtools_Pre_Processing/log_files/samtools_adults_pre_processing.stderr
#$ -t 60-119
#$ -tc 10

# Created by: Charley
# Date: 17.06.2022
# Aim: to perform de-duplication, sorting and indexing of bam files before methylation calling

########## note that the std out and std error will have to go to different files
# when we process the samples separately
# The -o and -e above need to be set accordingly


### Prep environment ###

module load samtools # v.1.9


### Set directories ###

# NB. this script expects all input bam files to be in the same dir 

BAM_DIR=/data/scratch/btx355/Chris_Grant_WGBS/bams # Alter path to input files as required
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Samtools_Pre_Processing/dedup_sorted_bams

cd $BAM_DIR


### Extract sample ID ###

# NB. using new sample list with hatchling 169-8 removed
# as it needs to be processed separately after BGI sort out reads issues
# -> Extract hatchlings with 1-59 and adults with 60-119

Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/FullSampleList_no169.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)

# Select correct line of list of files at each iteration
INPUT_FILE="${SAMPLE}.bam"


### Perform samtools processing steps ###

# Methylation calls for WGBS and targeted bisulfite sequencing 
# can be improved by removing PCR duplicate reads before calling methylation

# (1) Add mate score tag with samtools fixmate -m, in order to prep for duplicate removal, 
# with -p to disable proper pair check (as recommended on BSBolt docs)
# (2) Sort fixmates.bam file for methylation calling using samtools sort
# (3) Remove duplicate reads using samtools markdup
# Pipes used to suppress creation of intermediate files

samtools fixmate -@ ${NSLOTS} -p -m $INPUT_FILE - | \
samtools sort -@ ${NSLOTS} - | \
samtools markdup -@ ${NSLOTS} - - > $OUT_DIR/"${SAMPLE}.sorted.dedup.bam"

# (4) Index bam file using samtools index
samtools index -@ ${NSLOTS} $OUT_DIR/"${SAMPLE}.sorted.dedup.bam"

# (5) Generate alignment stats for bam files (optional)
STATS_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Samtools_Pre_Processing/stats

samtools stats -@ ${NSLOTS} \
$OUT_DIR/"${SAMPLE}.sorted.dedup.bam" > $STATS_DIR/"${SAMPLE}.sorted.dedup.bam.stats"

samtools flagstat -@ ${NSLOTS} \
$OUT_DIR/"${SAMPLE}.sorted.dedup.bam" > $STATS_DIR/"${SAMPLE}.sorted.dedup.bam.flagstat"

