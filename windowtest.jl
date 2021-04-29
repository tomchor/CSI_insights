using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields

factor = 32
arch = CPU()

topology = (Periodic, Periodic, Periodic)
grid = RegularRectilinearGrid(size=(8,8,8),
                              x=(0, 1),
                              y=(0, 1),
                              z=(-1, 0), 
                              topology=topology)
println("\n", grid, "\n")
#-----




# Set up ICs and/or Background Fields
#++++
kick = 1e-6
println("\nSetting geostrophic jet as an INITIAL CONDITION\n")
u_ic(x, y, z) = + kick*randn()
v_ic(x, y, z) = + kick*randn()

bg_fields = NamedTuple()
#-----


# Define model!
#++++
model_kwargs = (architecture = arch,
                grid = grid,
                )
model = IncompressibleModel(; model_kwargs..., closure=nothing)
println("\n", model, "\n")
#-----


# Adding the ICs
#++++
set!(model, u=u_ic, v=v_ic)
#-----


#+++++ Diagnostics
using Oceananigans.Diagnostics: WindowedSpatialAverage
u, v, w, = model.velocities

slicer = FieldSlicer(j=3:6, with_halos=false)
hor_window_average(F) = WindowedSpatialAverage(F; dims=(1, 2), field_slicer=slicer)
hor_mixed_average(F) = WindowedSpatialAverage(AveragedField(F; dims=1); dims=(1, 2), field_slicer=slicer)

mixed_u = hor_mixed_average(u)
window_u = hor_window_average(u)
outputs = (; mixed_u, window_u)
#----

# Define time-stepping
#++++
u_scale = maximum(u)
Δt = 1/5 * min(grid.Δx, grid.Δy) / u_scale
wizard = TimeStepWizard(cfl=0.8,
                        diffusive_cfl=0.5,
                        Δt=Δt, max_change=1.02, min_change=0.2, max_Δt=Inf, min_Δt=0.1seconds)
#----


# Finally define Simulation!
#++++
start_time = 1e-9*time_ns()
using Oceanostics: SingleLineProgressMessenger
simulation = Simulation(model, Δt=10seconds, 
                        stop_time=100seconds,
                        wall_time_limit=23.5hours,
                        iteration_interval=5,
                        progress=SingleLineProgressMessenger(initial_wall_time_seconds=start_time),
                        stop_iteration=Inf,)
println("\n", simulation, "\n")
#-----


# DIAGNOSTICS
#++++
using Oceananigans.OutputWriters: NetCDFOutputWriter, TimeInterval
simulation.output_writers[:snap] = NetCDFOutputWriter(model, outputs,
                                                        schedule = TimeInterval(20second),
                                                        filepath = "snap.window.nc",
                                                        mode = "c")
#-----


# Run the simulation!
#+++++
println("\n", simulation,
        "\n",)

@printf("---> Starting run!\n")
run!(simulation)
#-----
