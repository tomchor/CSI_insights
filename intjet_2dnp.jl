using ArgParse
using Printf
#using Pkg
#Pkg.add(path="/glade/u/home/tomasc/repos/Oceananigans.jl/", rev="tc/kernelcomputedfield")
using Oceananigans; oc = Oceananigans
using Oceananigans.Utils
using Oceananigans.Advection: WENO5, UpwindBiasedThirdOrder
using Oceananigans.OutputWriters, Oceananigans.Fields
using SpecialFunctions: erf


# Parse initial arguments
#++++
"Returns a dictionary of command line arguments."
function parse_command_line_arguments()
    settings = ArgParseSettings()
    @add_arg_table! settings begin

        "--factor"
            help = "Factor to divide Nh and Nz for"
            default = 16
            arg_type = Int

        "--arch"
            help = "CPU or GPU"
            default = :CPU
            arg_type = Symbol
    end
    return parse_args(settings)
end
args = parse_command_line_arguments()
factor = args["factor"]
if args["arch"] == :CPU
    arch = CPU()
elseif args["arch"] == :GPU
    arch = GPU()
else
    error("--arch keyword should be CPU or GPU")
end

@printf("Starting Oc with a dividing factor of %d and a %s architecture\n", 
        factor,
        arch,)
#-----


# Get simulation parameters
#++++
LES = false
include("jetinfo.jl")

simulation_nml = Simulation_list(f0=1e-4, Ny=4096, Nz=1024, Ly=8000, Lz=500).CIjet01
@unpack name, f0, u₀, N2_inf, Ny, Nz, Ly, Lz, σy, σz, y₀, z₀ = simulation_nml

simname = @sprintf("instFNN_%s", name)
#-----


# Calculate secondary parameters
#++++
b₀ = u₀ * f0
T_inertial = 2*π/f0

global_attributes = merge(simulation_nml,
                          (LES=Int(LES),
                          u0=u₀, y0=y₀, z0=z₀,
                          b0 = b₀, T_inertial = T_inertial,),
                         )
println("\n", global_attributes, "\n")
#-----


# Set GRID
#++++
Ny = Ny÷factor; Nz = Nz÷factor
Nx = 1
Lx = 6 * (Ly / Ny) * Nx
topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(size=(Nx, Ny, Nz),
                            x=(0, Lx),
                            y=(0, Ly),
                            z=(-Lz/2, +Lz/2), 
                            topology=topology)
println("\n", grid, "\n")
#-----

# Setting BCs
#++++
ubc = UVelocityBoundaryConditions(grid, 
                                  top = BoundaryCondition(Flux, 0),
                                  bottom = BoundaryCondition(Flux, 0),
                                 )
vbc = VVelocityBoundaryConditions(grid, 
                                  top = BoundaryCondition(Flux, 0),
                                  bottom = BoundaryCondition(Flux, 0),
                                 )
wbc = WVelocityBoundaryConditions(grid, 
                                 )
bbc = TracerBoundaryConditions(grid, 
                               bottom = BoundaryCondition(Gradient, N2_inf),
                               top = BoundaryCondition(Gradient, N2_inf),
                               )
#-----


# Set-up sponge layer
#++++
heaviside(X) = ifelse(X < 0, zero(X), one(X))
mask2nd(X) = heaviside(X) * X^2
mask3rd(X) = heaviside(X) * (-2*X^3 + 3*X^2)
const Hz = grid.Lz
const Hy = grid.Ly
const frac = 8

function bottom_mask(x, y, z)
    z₁ = -Hz/2; z₀ = z₁ + Hz/frac
    return mask3rd((z - z₀)/(z₁ - z₀))
end
function top_mask(x, y, z)
    z₁ = +Hz/2; z₀ = z₁ - Hz/frac
    return mask3rd((z - z₀)/(z₁ - z₀))
end
function north_mask(x, y, z)
    y₁ = Hy; y₀ = y₁ - Hy/frac
    return mask3rd((y - y₀)/(y₁ - y₀))
end
function south_mask(x, y, z)
    y₁ = 0; y₀ = y₁ + Hy/frac
    return mask3rd((y - y₀)/(y₁ - y₀))
end

full_mask(x, y, z) = bottom_mask(x, y, z) + top_mask(x, y, z) #+ north_mask(x, y, z) + south_mask(x, y, z)
full_sponge = Relaxation(rate=1/30minute, mask=full_mask, target=0)
#-----



# Set up ICs
#++++
using SpecialFunctions: erfc
const n2 = N2_inf

const Hz = grid.Lz
const Hy = grid.Ly
kick = 1e-4

gaussian(ψ) = exp(-ψ^2)
intgaussian(ψ) = √π/2 * (erf(ψ) + 1)
umask(Y, Z) = gaussian(Y) * gaussian(Z)
bmask(Y, Z) = (σy * intgaussian(Y)) * (-2Z * gaussian(Z) / σz)

u_ic(x, y, z) = u₀ * umask((y-y₀)/σy, (z-z₀)/σz) + kick*randn()
v_ic(x, y, z) = + kick*randn()
b_ic(x, y, z) = -b₀ * bmask((y-y₀)/σy, (z-z₀)/σz) + n2 * (z+Hz/2) + 1e-8*randn()
#-----


# Define our model!
#++++
import Oceananigans.TurbulenceClosures: AnisotropicDiffusivity, SmagorinskyLilly, IsotropicDiffusivity
if LES
    closure = SmagorinskyLilly()
else
#    closure = IsotropicDiffusivity(ν=1e-5, κ=1e-5)
    closure = AnisotropicDiffusivity(νh=5e-2, κh=5e-2, νz=5e-4, κz=5e-4)
end
model = IncompressibleModel(architecture = arch,
                            grid = grid,
                            advection = UpwindBiasedThirdOrder(),
                            #advection = WENO5(),
                            timestepper = :RungeKutta3,
                            closure = closure,
                            coriolis = FPlane(f=f0),
                            tracers = (:b,),
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (b=bbc, u=ubc, v=vbc, w=wbc),
                            forcing = (v=full_sponge, w=full_sponge,),
                            )
println("\n", model, "\n")
#-----



# Adding the IC
#++++
set!(model, u=u_ic, v=v_ic, b=b_ic)

v̄ = sum(model.velocities.v.data.parent) / (grid.Nx * grid.Ny * grid.Nz)
model.velocities.v.data.parent .-= v̄


if false
    using Plots; pyplot()
    x = oc.Grids.xnodes(Center, grid)
    y = oc.Grids.ynodes(Center, grid)
    z = oc.Grids.znodes(Center, grid)

    contourf(y, z, interior(model.tracers.b)[1,:,:]', levels=30)
end
#-----


# Define time-stepping and printing
#++++
u_scale = abs(u₀)
Δt = 0.1 * max(grid.Δx, grid.Δy) / u_scale
wizard = TimeStepWizard(cfl=0.8,
                        diffusive_cfl=0.8,
                        Δt=Δt, max_change=1.1, min_change=0.1, max_Δt=Inf, min_Δt=0.05seconds)

advCFL = oc.Diagnostics.AdvectiveCFL(wizard)
difCFL = oc.Diagnostics.DiffusiveCFL(wizard)
start_time = time_ns()
function progress(sim)
    msg = @printf("i: % 6d,    sim time: %10s,    wall time: %10s,    Δt: %10s,    diff CFL: %.2e,    adv CFL: %.2e\n",
                  sim.model.clock.iteration,
                  prettytime(sim.model.clock.time),
                  prettytime(1e-9 * (time_ns() - start_time)),
                  prettytime(sim.Δt.Δt),
                  difCFL(sim.model),
                  advCFL(sim.model),
                  )
    return msg
end
#-----

# Finally define Simulation!
#++++
simulation = Simulation(model, Δt=wizard, 
                        #stop_time=5*T_inertial,
                        stop_time=15*T_inertial,
                        iteration_interval=5, progress=progress,
                        stop_iteration=Inf,)
#-----




# START DIAGNOSTICS
#++++
import Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face

u, v, w = model.velocities
b = model.tracers.b
p = sum(model.pressures)

# Start calculation of snapshot variables
#++++
if LES
    νx, νy, νz = model.diffusivities.νₑ
else
    if :ν in fieldnames(typeof(model.closure))
        νx = νy = νz = model.closure.ν
    else
        νx = νy = model.closure.νx
        νz = model.closure.νz
    end
end

dbdz = ∂z(b)
ω_x = ∂y(w) - ∂z(v)
ω_y = ∂z(u) - ∂x(w)
ω_z = ∂x(v) - ∂y(u)

wb_res = @at (Center, Center, Center) w*b
wb_sgs = @at (Center, Center, Center) νz * dbdz

include("diagnostics.jl")

tke = KernelComputedField(Center, Center, Center, compute_turbulent_kinetic_energy!, model;
                          field_dependencies=(u, v, w, 0, 0, 0))

ε_ani = KernelComputedField(Center, Center, Center, compute_ani_viscous_dissipation!, model;
                            field_dependencies=(νx, νy, νz, u, v, w))

Ri = KernelComputedField(Center, Center, Face, compute_richardson_number!, model;
                         field_dependencies=(u, v, b, 0, 0, 0))

Ro = KernelComputedField(Face, Face, Center, compute_rossby_number!, model;
                         field_dependencies=(u, v, 0, 0, model.coriolis.f))

PV = KernelComputedField(Center, Center, Face, compute_pv_from_Ro_Ri!, model;
                         field_dependencies=(Ri, Ro, dbdz, model.coriolis.f))

dwpdz = KernelComputedField(Center, Center, Center, compute_vertical_pressure_term!, model;
                            field_dependencies=(w, p), parameters=1027)

outputs_netcdf = (u=u, 
                  v=v,
                  w=w,
                  b=b,
                  p=ComputedField(p),
                  wb_res=ComputedField(wb_res),
                  wb_sgs=ComputedField(wb_sgs),
                  dwpdz=ComputedField(dwpdz),
                  dbdz=ComputedField(dbdz),
                  ω_x=ComputedField(ω_x),
                  ω_z=ComputedField(ω_z),
                  tke=tke,
                  ε=ε_ani,
                  Ro=Ro,
                  Ri=Ri,
                  PV=PV,
                  )
#-----


# Video (low def) SNAPSHOTS
#++++
simulation.output_writers[:out_writer] =
    NetCDFOutputWriter(model, outputs_netcdf,
                       filepath = @sprintf("out.%s.nc", simname),
                       schedule = TimeInterval(6hours),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float64},
                       field_slicer = FieldSlicer(i=1),
                      )
#-----

# Analysis (high def) SNAPSHOTS
#++++
simulation.output_writers[:vid_writer] =
    NetCDFOutputWriter(model, outputs_netcdf,
                       filepath = @sprintf("vid.%s.nc", simname),
                       schedule = TimeInterval(60minutes),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float32},
                       field_slicer = FieldSlicer(i=1, j=1:Ny, k = 1+Nz÷frac : Nz-Nz÷frac),
                      )
#-----

# AVG outputs
#++++
avg_b = AveragedField(b, dims=(1, 2))
avg_tke = AveragedField(tke, dims=(1, 2))
avg_ε = AveragedField(ε_ani, dims=(1, 2))
avg_wb_res = AveragedField(wb_res, dims=(1, 2))
avg_wb_sgs = AveragedField(wb_sgs, dims=(1, 2))
avg_dwpdz = AveragedField(dwpdz, dims=(1, 2))

outputs = (b=avg_b,
           tke=avg_tke, 
           ε=avg_ε,
           wb_res=avg_wb_res,
           wb_sgs=avg_wb_sgs,
           dwpdz=avg_dwpdz,
          )
simulation.output_writers[:avg_writer] =
    NetCDFOutputWriter(model, outputs,
                       filepath = @sprintf("avg.%s.nc", simname),
                       schedule = AveragedTimeInterval(2minutes; window=1.9minutes, stride=1),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float64},
                      )
#-----
#-----


# Run the simulation!
#+++++
println("\n", simulation,
        "\n",)

@printf("---> Starting run!\n")
run!(simulation)
#-----

