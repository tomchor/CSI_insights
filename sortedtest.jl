using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields

factor = 32
arch = CPU()

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
                buoyancy = BuoyancyTracer(),
                tracers = :b,
                )
model = IncompressibleModel(; model_kwargs..., closure=nothing)
println("\n", model, "\n")
#-----


# Adding the ICs
#++++
set!(model, u=u_ic, v=v_ic, b=u_ic)
#-----

#++++ SortedField by operation
struct SortingOperand{F, D}
          field :: F
           dims :: D
end
SortedComputedField = ComputedField{X, Y, Z, S, <:RandomOperand} where {X, Y, Z, S}

function compute!(ψ::SortedComputedField)
    arch = architecture(ψ)
    ψ .= 0
    ψ .= ψ.operand.amplitude * CUDA.rand(size(ψ)...)
    fill_halo_regions!(ψ, arch)
    return nothing
end

struct RandomOperand
    amplitude :: Float64
end

const RandomComputedField = ComputedField{X, Y, Z, S, <:RandomOperand} where {X, Y, Z, S}

function compute!(ψ::RandomComputedField)
    arch = architecture(ψ)
    ψ .= 0
    ψ .= ψ.operand.amplitude * CUDA.rand(size(ψ)...)
    fill_halo_regions!(ψ, arch)
    return nothing
end

# Build random streamfunction
arch = GPU()
dψdt = ComputedField(Face, Face, Center, RandomOperand(forcing_amplitude), arch, grid)
#----



#++++ SortedField by custom output
function flattenedsort(A, dim_order::Union{Tuple, AbstractVector})
    return reshape(sort(permutedims(A, dim_order)[:]), (grid.Nx, grid.Ny, grid.Nz))
end

function sorted_b_output(model; average=false)
    b = model.tracers.b
    sorted_B = flattenedsort(interior(b), [3,2,1])
    if !average
        return sorted_B
    else
        return dropdims(mean(sorted_B, dims=(1,2)), dims=(1,2))
    end
end
mean_sorted_b_output = (mod)->sorted_b_output(mod; average=true)
#----


#+++++ Diagnostics
using Oceananigans.Diagnostics: WindowedSpatialAverage
u, v, w, = model.velocities
b = model.tracers.b

outputs = (b=b,
           mean_b=AveragedField(b, dims=(1,2)), 
           sorted_b=sorted_b_output,
           mean_sorted_b=mean_sorted_b_output,
          )
dims = Dict("sorted_b" => ("xC", "yC", "zC"), 
            "profile" => ("zC",),
            "mean_sorted_b" => ("zC",),
            );
#----


# Finally define Simulation!
#++++
start_time = 1e-9*time_ns()
using Oceanostics: SingleLineProgressMessenger
simulation = Simulation(model, Δt=10seconds, 
                        stop_time=100seconds,
                        iteration_interval=5,
                        progress=SingleLineProgressMessenger(initial_wall_time_seconds=start_time))
println("\n", simulation, "\n")
#-----


# DIAGNOSTICS
#++++
using Oceananigans.OutputWriters: NetCDFOutputWriter, TimeInterval
simulation.output_writers[:snap] = NetCDFOutputWriter(model, outputs,
                                                        schedule = TimeInterval(20second),
                                                        filepath = "snap.window.nc",
                                                        dimensions = dims,
                                                        mode = "c")
#-----


# Run the simulation!
#+++++
println("\n", simulation,
        "\n",)

@printf("---> Starting run!\n")
run!(simulation)
#-----
