#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition=prod_small
#SBATCH --time=0-02:00:00
##SBATCH --exclude bbpviz008
#SBATCH --job-name=Analyse
#SBATCH --account=proj55
#SBATCH --error=logs/analyse.log
#SBATCH --output=logs/analyse.log


#### SLURM Information
echo "SLURM information...
      partition: $MPIRUN_PARTITION
      number of nodes: $SLURM_JOB_NUM_NODES
      node id allocated: $SLURM_JOB_NODELIST
      number tasks per node: $SLURM_TASKS_PER_NODE
      directory from which sbatch was invoked: $SLURM_SUBMIT_DIR
      slurm memory binding: $SLURM_MEM_BIND
      slurm cpu binding: $SLURM_CPU_BIND
      job dependency: $SLURM_JOB_DEPENDENCY"

echo "You have been allocated: $SLURM_JOB_NODELIST"

set -e
set -x

module purge all
#module load bbplxviz/cscs-viz
#module load cmake
#module load slurm


export MODULEPATH=/gpfs/bbp.cscs.ch/project/proj55/software/bb5/local-20180425105504/share/modules:${MODULEPATH}; module load proj55

#export MODULEPATH=/gpfs/bbp.cscs.ch/project/proj55/software/cscsviz/local-20180425110422/share/modules:${MODULEPATH}; module load proj55
#module load /gpfs/bbp.cscs.ch/project/proj55/software/viz/local/share/modules/optframework-setup
#source /gpfs/bbp.cscs.ch/project/proj55/iavarone/software/workspace_BPO/venv/bin/activate
#module load /gpfs/bbp.cscs.ch/home/roessert/software/BluePyExtract/load_module


./analyse.sh 042bb46 1,2,3
#./analyse.sh 6730649 1,2,3
#./analyse.sh 809fecc 1,2,3
#./analyse.sh 4ff9be2 1,2,3
#./analyse.sh 760c8c4 1,2,3
#./analyse.sh 680bc29 1,2,3


wait
