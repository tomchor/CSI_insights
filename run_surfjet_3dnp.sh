#!/bin/bash -l
#SBATCH --account=UMCP0012
#SBATCH --job-name=sj3d_SIjet4
#SBATCH --output=out/slurm.out
#SBATCH --error=out/%x.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:v100:1 # GPUs per Node
#SBATCH --mail-user=tchor@umd.edu
#SBATCH --mail-type=ALL

# Clear the environment from any previously loaded modules
module purge
module load gnu
module load cuda

/glade/u/apps/ch/opt/usr/bin/dumpenv # Dumps environment (for debugging with CISL support)

/glade/u/home/tomasc/repos/julia/julia --project -J ../Ocean_cuda_v0.48.0.so \
    surfjet_3dnp.jl --jet=SIjet4 --arch=GPU --factor=4 2>&1 | tee out/sj3d_SIjet4.out

