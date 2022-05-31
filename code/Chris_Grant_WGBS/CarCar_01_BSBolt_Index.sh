#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=128G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea
#$ -l highmem


# Created by: Charley
# Date: 30.05.2022
# Aim: Create BSBolt WGBS index for CarCar assembly


### Prep environment ###

module load python # Python v.3.8.5
source /data/SBCS-EizaguirreLab/Charley/environments/my_bsbolt/bin/activate # BSBolt v.1.5.0

OUTDIR=/data/scratch/btx902/Chris_Grant_WGBS

cd $OUTDIR


### WGBS index generation ###

python3 -m bsbolt Index \
-G CarCar_v1_2021_12.fasta \
-DB $OUTDIR/CarCar_v1_2021_12_bsbolt_DB

deactivate


