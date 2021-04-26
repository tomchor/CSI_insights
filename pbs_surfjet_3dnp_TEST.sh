#!/bin/bash -l
#PBS -A UMCP0012
#PBS -N sj3dSIjet4
#PBS -k eod
#PBS -o logs/sj3d_SIjet4.out
#PBS -e logs/sj3d_SIjet4.err
#PBS -l walltime=24:00:00
#PBS -q casper
#PBS -l select=1:ncpus=1:ngpus=1
#PBS -l gpu_type=v100
#PBS -M tchor@umd.edu
#PBS -m abe

# Clear the environment from any previously loaded modules
module purge
module load gnu
module load cuda/11.0.3

#/glade/u/apps/ch/opt/usr/bin/dumpenv # Dumps environment (for debugging with CISL support)

/glade/u/home/tomasc/repos/julia_1.5.2/julia --project \
    test3d.jl --jet=SIjet4 --arch=GPU --factor=2 2>&1 | tee out/sj3d_SIjet4.out

