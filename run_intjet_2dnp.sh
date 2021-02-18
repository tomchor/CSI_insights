#!/bin/bash -l
#SBATCH --account=UMCP0012
#SBATCH --job-name=jet_2dnp
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=8:00:00
##SBARCH --mem=0
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:v100:1 # GPUs per Node
##SBATCH --cpus-per-task=4
#SBATCH --mail-user=tchor@umd.edu
#SBATCH --mail-type=ALL

# Clear the environment from any previously loaded modules
module purge
module load gnu
module load cuda

/glade/u/apps/ch/opt/usr/bin/dumpenv # Dumps environment (for debugging with CISL support)

/glade/u/home/tomasc/repos/julia/julia jet_2dnp.jl --arch=GPU --factor=1 2>&1 | tee output.out

#module load julia
#julia inst_2dnp.jl --arch=GPU --factor=1 2>&1 | tee output.out
