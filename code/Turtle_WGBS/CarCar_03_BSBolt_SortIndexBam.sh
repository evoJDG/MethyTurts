#!/bin/bash
# Samtools sort index for BSBolt methylation call
#$ -N samtools
#$ -l h_vmem=32G
#$ -pe smp 4
#$ -l h_rt=240:0:0
#$ -o /data/SBCS-EizaguirreLab/Turtle_WGBS/04_Samtools_Pre_Processing/run_samtools_adults.stdout
#$ -e /data/SBCS-EizaguirreLab/Turtle_WGBS/04_Samtools_Pre_Processing/run_samtools_adults.stderr
#$ -V
#$ -t 61-120
#$ -tc 20


########## note that the std out and std error will have to go to different files
# when we process the samples separately
# The -o and -e above need to be set accordingly

module load samtools # samtools version 1.9

### extract sample ID
Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/FullSampleList.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)

########### this also needs to be amended for data stored in different dirs
BAMDIR=/data/scratch/btx355/Chris_Grant_WGBS/Raw_Data/$SAMPLE

cd $BAMDIR

# Select the correct line of list of files at each iteration
INPUT_FILE=${SAMPLE}.bam

# Methylation calls for WGBS and targeted bisulfite sequencing 
# can be improved by removing PCR duplicate reads before calling methylation.

# fixmates to prepare for duplicate removal, use -p to disable proper pair check
samtools fixmate -p -m $INPUT_FILE ${SAMPLE}.fixmates.bam

# sort then index bam file for methylation calling
samtools sort -@ 4 # set 4 threads\ 
-o ${SAMPLE}.sorted.bam\
 ${SAMPLE}.fixmates.bam

# remove duplicate reads
samtools markdup ${SAMPLE}.sorted.bam ${SAMPLE}.dup.bam

# index bam file for methylation calling
samtools index ${SAMPLE}.dup.bam
