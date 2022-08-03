#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=128G
#$ -l h_rt=24:00:00
#$ -l highmem
#$ -m bea
#$ -o run_R01.1_prepBSBoltForMethylkit.R.stdout

module load R/4.1.1

SCRIPT_DIR=/data/SBCS-EizaguirreLab/Turtle_WGBS/06_MethylKit/Hatchery_CG_Meth_Calls/BSBolt_mincov1/scripts

Rscript $SCRIPT_DIR/R01.1_prepBSBoltForMethylkit.R
