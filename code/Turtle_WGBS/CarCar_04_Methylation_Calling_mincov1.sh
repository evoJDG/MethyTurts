#!/bin/bash
#$ -pe smp 6
#$ -l h_vmem=3G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -t 1-58


# Created by: Charley
# Date: 02.08.2022
# Aim: Perform BSBolt methylation calling on the sorted.dedup.bam files

# BSBolt methylation caller filters by 10X minimum coverage by default!
# This script filters by 1X minimum coverage

# This script also only outputs CpG sites using -CG option
# As it runs faster and don't need to perform meth splitting after
# Change if other C contexts are needed


### Prep environment ###

module load python # v.3.8.5
source /data/SBCS-EizaguirreLab/James/environments/my_bsbolt/bin/activate # v.1.5.0


### Extract sample IDs ###

Sample_File=/data/scratch/btx902/Turtle_WGBS/trimmed_fastqs/FullSampleList.txt

SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)


### Set directories ###

BAM_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Samtools_Pre_Processing/dedup_sorted_bams
DB=/data/SBCS-EizaguirreLab/Turtle_WGBS/02_BSBolt_Index/CarCar_v1_2021_12_bsbolt_DB
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/05_Methylation_Calling/meth_calls/mincov1/$SAMPLE
STATS_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/05_Methylation_Calling/stats/mincov1

cd $OUT_DIR


### Perform methylation calling ###

# Methylation calling is performed by counting the number of bisulfite 
# converted bases relative to the number of reads observed at each cytosine.
# Relative to the reference genome, methylation status at a cytosine and guanine can only
# be called using reads mapped to Watson and Crick strands respectively.

INPUT_FILE=$BAM_DIR/"${SAMPLE}.sorted.dedup.bam"

bsbolt CallMethylation \
-I $INPUT_FILE \
-O $OUT_DIR/"BSBolt_MethCall_${SAMPLE}" \
-DB $DB \
-CG -verbose -t ${NSLOTS} \
-min 1 > $STATS_DIR/"${SAMPLE}_BSBolt_methylation_stats.txt"

