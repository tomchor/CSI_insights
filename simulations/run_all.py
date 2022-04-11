from os import system
snames = [#"PNN_CIfront1",
          #"PNN_CIfront2",
          #"PNN_CIfront3",
          #"PNN_CIfront4",
          #"PNN_CIfront5",
          #"PNN_SIfront1",
          #"PNN_SIfront2",
          #"PNN_SIfront3",
          #"PNN_SIfront4",
          #"PNN_SIfront5",
          #"PNN_SIfront6",
          #"PNN_CIintjet01",
          #"FNN_CIfront1",
          #"FNN_CIfront3",
          #"FNN_SIfront4",
          #"PNN_CIfront1_f2",
          #"PNN_CIfront1_f4",
          #"PNN_CIfront1_f8",
          #"PNN_SIfront4_f2",
          #"PNN_SIfront4_f4",
          #"PNN_SIfront4_f8",
          #"PNN_CIfront1_AMD",
          "PNN_CIfront1_AMD_f2",
          "PNN_CIfront1_AMD_f4",
          "PNN_CIfront1_AMD_f8",
          #"PNN_SIfront4_AMD",
          "PNN_SIfront4_AMD_f2",
          "PNN_SIfront4_AMD_f4",
          "PNN_SIfront4_AMD_f8",
          ]

verbose = 1
aux_filename = "aux_pbs_ISI_jet.sh"
data_dir = "data/"
remove_checkpoints = False

pbs_script = \
r"""#!/bin/bash -l
#PBS -A UMCP0012
#PBS -N {0}
#PBS -k eod
#PBS -o logs/{0}.out
#PBS -e logs/{0}.err
#PBS -l walltime=24:00:00
#PBS -q casper
#PBS -l select=1:ncpus=1:ngpus=1
#PBS -l gpu_type=v100
#PBS -M tchor@umd.edu
#PBS -m abe

# Clear the environment from any previously loaded modules
module purge
module load ncarenv/1.3 gnu/9.1.0 ncarcompilers/0.5.0
module load netcdf/4.7.4 openmpi/4.1.0 julia/1.6.0
module load peak_memusage
module li

#/glade/u/apps/ch/opt/usr/bin/dumpenv # Dumps environment (for debugging with CISL support)

export JULIA_DEPOT_PATH="/glade/work/tomasc/.julia"

peak_memusage.exe julia --project \
    surfjet_np.jl --simname={0} --factor=1 2>&1 | tee out/{0}.out
"""

for sname in snames:
    if remove_checkpoints:
        cmd0 = f"rm data/chk.{sname}*.nc"
        if verbose>0: print(cmd0)
        system(cmd0)
    auxfile1 = pbs_script.format(sname)
    if verbose>1: print(auxfile1)

    with open(aux_filename, "w") as f:
        f.write(auxfile1)

    cmd1 = f"qsub {aux_filename}"
    if verbose>0: print(cmd1)
    system(cmd1)

    print()
