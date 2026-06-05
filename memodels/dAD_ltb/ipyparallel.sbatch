#!/bin/bash

#SBATCH --account=proj55
#SBATCH --error=logs/ipyparallel_%j.log
#SBATCH --output=logs/ipyparallel_%j.log

##SBATCH --exclude=bbpviz001,bbpviz002,bbpviz003,bbpviz008,bbpviz009
##SBATCH --exclude=bbpcj027

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

module load slurm

if [[ ${HOSTNAME} =~ "r" ]]; then
  export MODULEPATH=/gpfs/bbp.cscs.ch/project/proj55/software/bb5/local-20180425105504/share/modules:${MODULEPATH}; module load proj55
else
    module purge all
    module load bbplxviz/cscs-viz
    module load cmake
    export MODULEPATH=/gpfs/bbp.cscs.ch/project/proj55/software/cscsviz/local-20180425110422/share/modules:${MODULEPATH}; module load proj55
fi


#module load /gpfs/bbp.cscs.ch/project/proj55/software/viz/local/share/modules/optframework-setup
#source /gpfs/bbp.cscs.ch/project/proj55/iavarone/software/workspace_BPO/venv/bin/activate
python -c 'import IPython; print IPython.__file__'

cd ${RUNDIR}

export USEIPYP=1
export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=ipyparallel.${SLURM_JOBID}
ipcontroller --init --ip='*' --sqlitedb --profile=${IPYTHON_PROFILE} &
sleep 30
srun ipengine --profile=${IPYTHON_PROFILE} &

CHECKPOINTS_DIR="./checkpoints/run.${GITHASH}"
mkdir -p ${CHECKPOINTS_DIR}

#rm -rf x86_64
#/gpfs/bbp.cscs.ch/project/proj55/iavarone/software/workspace_BPO/install/nrnpython/x86_64/bin/nrnivmodl mechanisms
nrnivmodl mechanisms

sleep 20

pids=""
for seed in `seq 1 ${MAX_SEED}`
do
    python opt_model.py --etype=${ETYPE} --start --offspring_size=${OFFSPRING_SIZE} --seed=${seed} --max_ngen=${MAX_NGEN} --githash=${GITHASH} --slurmid=${SLURM_JOBID} --runlog=${RUNLOGPATH} &
    pids="${pids} $!"
done

#./analyse.sh ${GITHASH} 1,2,3
wait $pids
