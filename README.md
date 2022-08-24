# Supporting code for "Insights into the mixing efficiency of submesoscale Centrifugal-Symmetric instabilities"

This repository contains the code necessary for generating the data in the paper "Insights into the
mixing efficiency of submesoscale Centrifugal-Symmetric instabilities". The code is fully written
in Julia, the main package is
[Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/), and simulations output
NetCDF files.

## Organization

- Files `Project.toml` and `Manifest.toml` ensure a reproducible Julia environment that should
  generate identical results to those in the paper, apart from the randomness in the initial noise and
machine imprecision.

- File `surfjet_np.jl` is the main script (it defines the grid, model, simulation, and finally runs
  the simulation)

- File `diagnostics.jl` is an auxiliary file that contains diagnostics that will be written to the
  output netcdf files

- File `jetinfo.jl` is an auxiliary file that contains the parameters for each simulation



## Running the simulations

Before running the simulations you'll need to download and install the necessary packages (which are
described in `Project.toml` and `Manifest.toml`). This is done automatically by Julia(instructions
[here](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)) and should not
be done manually!

After the packages are properly installed, in order to properly run the simulations you'll need to
specify a few options. The command below has an example:

```
julia --project surfjet_np.jl --simname=PNN_CIfront1 --factor=1
```

The flag `--project` is necessary to start Julia with the package versions defined in `Project.toml`
and `Manifest.toml`.

`surfjet_np.jl` is the main script.

`--simname=PNN_CIfront1` is the name of the simulation you want to run. The prefix `PNN` means it
should be run with a periodic, nonperiodic, nonperiodic. The other option is `FNN` which means flat,
nonperiodic, nonperiodic (which runs a 2D simulation). The rest of the name needs to be defined in
the file `jetinfo.jl`.

`--factor=1` is the "coarsening factor". Setting a value of 1 runs a simulation with full
resolution. Factor of 2 halves the grid spacing. Factor of 4 divides the spacings by 4 and so on.


This code will run either on a CPU or on a GPU, however, the full (i.e. factor=1) simulations are
large (~130 million points) so they are too large to run on a CPU and a GPU is needed.


