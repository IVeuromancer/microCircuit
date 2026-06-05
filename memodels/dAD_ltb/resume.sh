#!/bin/bash


export GITHASH=71668e0
export RUNDIR="./run/${GITHASH}"

export RUNLOGPATH="`pwd`/runlog.json"
export MAX_NGEN=100
export MAX_SEED=3

JOBNAME="TCdAD_${GITHASH}"

#sbatch -n 64 -p cloud -t 2-00:00:00 -J ${JOBNAME} ipyparallel_resume.sbatch

#sbatch -n 244 -p prod -t 14:22:00 -C cpu -J ${JOBNAME} ipyparallel_resume.sbatch
sbatch -n 288 -p prod -t 1-00:00:00 --begin=2019-03-25T08:05:00 -C cpu -J ${JOBNAME} ipyparallel_resume.sbatch

#sbatch -n 32 -p prod -t 3-00:00:00 --reservation=singlecell-optimization -J ${JOBNAME} ipyparallel_resume.sbatch
