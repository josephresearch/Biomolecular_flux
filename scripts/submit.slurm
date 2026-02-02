#!/bin/bash
#SBATCH --job-name=flux_run
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=05:10:00
#SBATCH --mail-user=yw9071@princeton.edu
#SBATCH --account=jerelle


module purge
module load intel/2024.2
module load intel-mkl/2024.2
module load intel-mpi/gcc/2021.13
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


srun $HOME/.local/bin/lmp_Aug2024 -in flux_run.lmp
