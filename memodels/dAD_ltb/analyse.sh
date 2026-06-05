#!/bin/bash

set -e
set -x

module purge all
module load bbplxviz/cscs-viz
module load cmake
module load slurm
#source /gpfs/bbp.cscs.ch/project/proj55/iavarone/software/workspace_BPO/venv/bin/activate
#module load /gpfs/bbp.cscs.ch/project/proj38/roessert/software/viz/local/share/modules/optframework-setup
#module load /gpfs/bbp.cscs.ch/project/proj55/software/viz/local/share/modules/optframework-setup
#export MODULEPATH=/gpfs/bbp.cscs.ch/project/proj55/software/cscsviz/local-20180425110422/share/modules:${MODULEPATH}; module load proj55
export MODULEPATH=/gpfs/bbp.cscs.ch/project/proj55/software/bb5/local-20180425105504/share/modules:${MODULEPATH}; module load proj55
if [ "$#" -eq 2 ]; then
  RANK=0
else
  RANK=$3
fi

#python analyse.py --rundir="run/$1" --githash=$1 --seed=$2 --rank=${RANK} --usethreshold --stochdet=False --mmtest #--altmorph rp160229_A_idE
python analyse.py --rundir="run/$1" --githash=$1 --seed=$2 --rank=${RANK} --stochdet=False
# C290999C-I4 C070301B2 C080998A --live
# Bad Morph: sm110128a1-2_INT_idA
