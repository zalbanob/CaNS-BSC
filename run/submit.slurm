#!/bin/bash
#SBATCH --job-name=cans
#SBATCH --account=ehpc1022
#SBATCH --qos=acc_ehpc
#SBATCH --partition=acc
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:4
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err

set -euo pipefail

module purge
module load nvidia-hpc-sdk/24.3
module load cuda/12.3

# cuDecomp runtime
export LD_LIBRARY_PATH=/gpfs/scratch/ehpc1022/CaNS/dependencies/cuDecomp/build/lib:${LD_LIBRARY_PATH:-}


export SLURM_CPU_BIND=none
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}


mpirun --bind-to none -np ${SLURM_NTASKS} ./cans
