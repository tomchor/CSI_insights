#!/bin/bash -l
#SBATCH --account=UMCP0012
#SBATCH --job-name=i2CIjet01
#SBATCH --output=out/slurm.out
#SBATCH --error=out/%x.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --gres=gpu:v100:1 # GPUs per Node
#SBATCH --mail-user=tchor@umd.edu
#SBATCH --mail-type=ALL

# Clear the environment from any previously loaded modules
module purge
module load gnu
module load cuda

/glade/u/apps/ch/opt/usr/bin/dumpenv # Dumps environment (for debugging with CISL support)

/glade/u/home/tomasc/repos/julia/julia --project \
    intjet_2dnp.jl --jet=CIjet01 --arch=GPU --factor=1 2>&1 | tee out/ij2d_CIjet01.out

