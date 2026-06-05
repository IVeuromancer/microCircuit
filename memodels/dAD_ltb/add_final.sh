#!/bin/bash

set -e
set -x

if [ "$#" -eq 2 ]; then
  BRANCH=`date +%Y%m%d%H%M%S`
  #git branch ${BRANCH} $1
  git branch ${BRANCH}
else
  BRANCH=$3
fi

kinit

python analyse.py --rundir="./run/$1" --githash=$1 --seed=$2 --branch=${BRANCH} --mmtest #--nosimulate

#git push --all -u
