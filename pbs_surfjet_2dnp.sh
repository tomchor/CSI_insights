#!/bin/bash -l
#PBS -A UMCP0012
#PBS -N sj2dCIjet1
#PBS -k eod
#PBS -o logs/sj2d_CIjet1.out
#PBS -e logs/sj2d_CIjet1.err
#PBS -l walltime=16:00:00
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
    surfjet_2dnp.jl --jet=CIjet1 --arch=GPU --factor=1 2>&1 | tee out/sj2d_CIjet1.out

