#!/bin/bash -l
#PBS -A UMCP0012
#PBS -N sj2dSIjet4
#PBS -k eod
#PBS -o logs/sj2d_SIjet4.out
#PBS -e logs/sj2d_SIjet4.err
#PBS -l walltime=16:00:00
#PBS -q casper
#PBS -l select=1:ncpus=1:ngpus=1
#PBS -l gpu_type=v100
#PBS -M tchor@umd.edu
#PBS -m abe

# Clear the environment from any previously loaded modules
module purge
module load gnu
module load cuda

#/glade/u/apps/ch/opt/usr/bin/dumpenv # Dumps environment (for debugging with CISL support)

/glade/u/home/tomasc/repos/julia/julia --project \
    surfjet_2dnp.jl --jet=SIjet4 --arch=GPU --factor=2 2>&1 | tee out/sj2d_SIjet4.out

