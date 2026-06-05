#!/bin/bash

#SBATCH --account=proj55
#SBATCH --error=logs/ipyparallel_%j.log
#SBATCH --output=logs/ipyparallel_%j.log

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
module load bbplxviz/cscs-viz
module load cmake
module load slurm

export MODULEPATH=/gpfs/bbp.cscs.ch/project/proj55/software/bb5/local-20180425105504/share/modules:${MODULEPATH}; module load proj55

python -c 'import IPython; print IPython.__file__'

cd ${RUNDIR}

export USEIPYP=1
export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=ipyparallel.${SLURM_JOBID}
ipcontroller --init --ip='*' --sqlitedb --profile=${IPYTHON_PROFILE} &
sleep 30
srun ipengine --profile=${IPYTHON_PROFILE} &

pids=""
for seed in `seq 1 ${MAX_SEED}`
do
    python opt_model.py --resume --seed=${seed} --githash=${GITHASH} --max_ngen=${MAX_NGEN} --slurmid=${SLURM_JOBID} --runlog=${RUNLOGPATH} &
    pids="${pids} $!"
done

wait $pids
