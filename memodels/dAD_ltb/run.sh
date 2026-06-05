#!/bin/bash

ETYPE=cAD_noscltb_legacy

export OFFSPRING_SIZE=100
export MAX_NGEN=100
export MAX_SEED=3
export ETYPE

git add -A && git commit --allow-empty -a -m "Running optimization ${ETYPE}"
export GITHASH=$(git rev-parse --short HEAD)
export RUNLOGPATH="`pwd`/runlog.json"
export RUNDIR="./run/${GITHASH}"
git archive --format=tar --prefix=${GITHASH}/ HEAD | (cd ./run/ && tar xf -)


sbatch -n 288 -p prod -t 1-00:00:00 -C cpu -J ${ETYPE} ipyparallel.sbatch
#sbatch -n 136 -p prod -t 3-00:00:00 -C cpu --qos=longjob -J ${ETYPE} ipyparallel.sbatch

