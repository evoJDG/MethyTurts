#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=13G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -l highmem
#$ -t 1-60


# Created by: Charley
# Date: 10.06.2022
# Aim: Run BSBolt Align for all 120 samples (QC'ed and trimmed)
# Sample list location: /data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/FullSampleList.txt

### Prep environment ###
module unload anaconda3
module load python # Python v.3.8.5
source /data/SBCS-EizaguirreLab/James/environments/my_bsbolt/bin/activate # BSBolt v.1.5.0

### extract sample ID
Sample_File=/data/SBCS-EizaguirreLab/Turtle_WGBS/00_General/FullSampleList.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)

## set directories
DIR=/data/scratch/btx355/Chris_Grant_WGBS/Raw_Data/$SAMPLE # Note this has to be changed depending on user
DB=/data/SBCS-EizaguirreLab/Turtle_WGBS/02_BSBolt_Index/CarCar_v1_2021_12_bsbolt_DB/

### Run BSBolt Align for paired end WGBS reads ###

# for each sample, there are data from three or four lanes
# combine files from different lanes but keep paired ends separate

cat $DIR/*_1_trim.fq.gz> $DIR/${SAMPLE}_1_trim_merged.fq.gz

cat $DIR/*_2_trim.fq.gz> $DIR/${SAMPLE}_2_trim_merged.fq.gz

echo $SAMPLE 'merged' `date` >> /data/SBCS-EizaguirreLab/Turtle_WGBS/03_BSBolt_Alignments/SamplesMerged.txt

python3 -m bsbolt Align \
-F1 $DIR/${SAMPLE}_1_trim_merged.fq.gz \
-F2 $DIR/${SAMPLE}_2_trim_merged.fq.gz \
-DB $DB \
-O $DIR/$SAMPLE \
-t ${NSLOTS}

echo $SAMPLE 'aligned' `date` >> /data/SBCS-EizaguirreLab/Turtle_WGBS/03_BSBolt_Alignments/SamplesAligned.txt
