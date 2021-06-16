using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields
using CUDA: has_cuda

topology = (Periodic, Periodic, Periodic)
grid = RegularRectilinearGrid(size=(3,3,3),
                              x=(0, 1),
                              y=(0, 1),
                              z=(-1, 0), 
                              topology=topology)
println("\n", grid, "\n")
#-----




# Set up ICs and/or Background Fields
#++++
kick = 1e-3
u_ic(x, y, z) = 0.1 + kick*randn()
v_ic(x, y, z) = 0.05 + kick*randn()
u_g(x, y, z, t) = u_ic(x, y, z)
b_g(x, y, z, t) = v_ic(x, y, z)
#-----


# Define model!
#++++
model_kwargs = (architecture = has_cuda() ? GPU() : CPU(),
                grid = grid,
                buoyancy = BuoyancyTracer(),
                tracers = :b,
                )
model = IncompressibleModel(; model_kwargs..., closure=IsotropicDiffusivity())
println("\n", model, "\n")
#-----


# Adding the ICs
#++++
set!(model, u=u_ic, v=v_ic, b=u_ic)
#-----


#+++++ Diagnostics
using Oceananigans.Diagnostics: WindowedSpatialAverage
using Oceanostics.TurbulentKineticEnergyTerms: IsotropicViscousDissipationRate, IsotropicPseudoViscousDissipationRate
u, v, w = model.velocities

u_win = WindowedSpatialAverage(u; dims=(1, 2))
u_avg = AveragedField(u, dims=(1, 2))


ε = IsotropicPseudoViscousDissipationRate(model, u, v, w, model.closure.ν)
ε_win = WindowedSpatialAverage(ε, dims=(1,2))
ε_avg = AveragedField(ε, dims=(1, 2))

outputs = (; u_avg, u_win, ε_avg, ε_win)
#----


# Finally define Simulation!
#++++
start_time = 1e-9*time_ns()
using Oceanostics: SingleLineProgressMessenger
simulation = Simulation(model, Δt=1seconds, 
                        stop_time=17hours,
                        iteration_interval=5,
                        progress=SingleLineProgressMessenger(initial_wall_time_seconds=start_time))
println("\n", simulation, "\n")
#-----


# DIAGNOSTICS
#++++
include("diagnostics_test.jl")
as_background=false; f_0=1e-4; ρ₀=1023
global_attributes=NamedTuple()
checkpointer = construct_outputs(model, simulation, LES=false, simname="avgtest", frac=0)
#-----


# Run the simulation!
#+++++
println("\n", simulation,
        "\n",)

@printf("---> Starting run!\n")
run!(simulation)
#-----
