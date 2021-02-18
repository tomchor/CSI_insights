using ArgParse
using Printf
using Oceananigans; oc = Oceananigans
using Oceananigans.Utils
using Oceananigans.Advection: UpwindBiasedThirdOrder
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
            default = 32
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
as_background=true
include("jetinfo.jl")

simulation_nml = SurfaceJetSimulations().SIjet4
@unpack name, f0, u₀, N2_inf, N2_pyc, Ny, Nz, Ly, Lz, σy, σz, y₀, z₀ = simulation_nml

simname = @sprintf("instPNN_%s", name)
#-----


# Calculate secondary parameters
#++++
b₀ = u₀ * f0
T_inertial = 2*π/f0

global_attributes = merge(simulation_nml,
                          (
                          u_0=u₀, y_0=y₀, z_0=z₀,
                          b0 = b₀, T_inertial = T_inertial,),
                         )
println("\n", global_attributes, "\n")
#-----


# Set GRID
#++++  GRID
Nx = Ny÷32
Lx = (Ly / Ny) * Nx
topology = (Periodic, Bounded, Bounded)
grid = RegularCartesianGrid(size=(Nx÷factor, Ny÷factor, Nz÷factor),
                            x=(0, Lx),
                            y=(0, Ly),
                            z=(-Lz, 0), 
                            topology=topology)
println("\n", grid, "\n")
#-----


# Set up Geostrophic flow
#++++++
heaviside(X) = ifelse(X < 0, zero(X), one(X))
const n2_inf = N2_inf
const n2_pyc = N2_pyc
const Hz = grid.Lz
const Hy = grid.Ly
const sig_z = σz
const sig_y = σy
const u_0 = u₀
const y_0 = y₀
const z_0 = z₀
const z_c = -40
const z_m = z_c - n2_pyc/n2_inf*(z_c+Hz)
const f_0 = f0
fy(ψ) = exp(-ψ^2)
intgaussian(ψ) = √π/2 * (erf(ψ) + 1)
umask(Y, Z) = Z * fy(Y)
bmask(Y, Z) = (1/sig_z) * (sig_y * intgaussian(Y))

u_g(x, y, z, t) = +u_0 * umask((y-y_0)/sig_y, ((z-z_0)/sig_z +1))
background_strat(z) = ifelse(z < z_c, 
                             n2_pyc * (z+Hz),
                             n2_inf * (z-z_m))
b_g(x, y, z, t) = -f_0 * u_0 * bmask((y-y_0)/sig_y, ((z-z_0)/sig_z +1)) + background_strat(z)
dudz_g(x, y, z, t) = +u_0 * (1/sig_z) * fy((y-y_0)/sig_y)
#-----

# Setting BCs
#++++
if as_background
    @inline surface_grad(x, y, t) = -dudz_g(x, y, 0, t)
    @inline bottom_grad(x, y, t) = -dudz_g(x, y, -Hz, t)
    U_top_bc = GradientBoundaryCondition(surface_grad)
    U_bot_bc = GradientBoundaryCondition(bottom_grad)
    B_bc = GradientBoundaryCondition(0)
else
    U_top_bc = FluxBoundaryCondition(0)
    U_bot_bc = FluxBoundaryCondition(0)
    B_bc = GradientBoundaryCondition(N2_inf)
end

ubc = UVelocityBoundaryConditions(grid, 
                                  top = U_top_bc,
                                  bottom = BoundaryCondition(Flux, 0),
                                  )
vbc = VVelocityBoundaryConditions(grid, 
                                  top = BoundaryCondition(Flux, 0),
                                  bottom = BoundaryCondition(Flux, 0),
                                 )
wbc = WVelocityBoundaryConditions(grid, 
                                 )
bbc = TracerBoundaryConditions(grid, 
                               bottom = B_bc,
                               top = B_bc,
                               )
#-----


# Set-up sponge layer
#++++
mask2nd(X) = heaviside(X) * X^2
mask3rd(X) = heaviside(X) * (-2*X^3 + 3*X^2)
const Hy = grid.Ly
const frac = 8

function bottom_mask(x, y, z)
    z₁ = -Hz; z₀ = z₁ + Hz/frac
    return mask2nd((z - z₀)/(z₁ - z₀))
end
function top_mask(x, y, z)
    z₁ = +Hz; z₀ = z₁ - Hz/frac
    return mask2nd((z - z₀)/(z₁ - z₀))
end
function north_mask(x, y, z)
    y₁ = Hy; y₀ = y₁ - Hy/frac
    return mask2nd((y - y₀)/(y₁ - y₀))
end
function south_mask(x, y, z)
    y₁ = 0; y₀ = y₁ + Hy/frac
    return mask2nd((y - y₀)/(y₁ - y₀))
end

full_mask(x, y, z) = north_mask(x, y, z) + south_mask(x, y, z)# + bottom_mask(x, y, z)
full_sponge = Relaxation(rate=1/10minute, mask=full_mask, target=0)
#-----



# Set up ICs and/or Background Fields
#++++
kick = 1e-6
if as_background
    println("\nSetting geostrophic jet as BACKGROUND\n")
    u_ic(x, y, z) = 0 #+ kick*randn()
    v_ic(x, y, z) = 0 #+ kick*randn()
    b_ic(x, y, z) = + 1e-8*randn()

    bg_fields = (u=u_g, b=b_g,)
else
    println("\nSetting geostrophic jet as an INITIAL CONDITION\n")
    u_ic(x, y, z) = u_g(x, y, z, 0) + kick*randn()
    v_ic(x, y, z) = + kick*randn()
    b_ic(x, y, z) = b_g(x, y, z, 0) + 1e-8*randn()

    bg_fields = NamedTuple()
end
#-----


# Define our model!
#++++
import Oceananigans.TurbulenceClosures: SmagorinskyLilly, AnisotropicMinimumDissipation
import Oceananigans.TurbulenceClosures: AnisotropicDiffusivity, IsotropicDiffusivity
closure_LES = SmagorinskyLilly(C=0.23)
#closure = AnisotropicMinimumDissipation()
closure_cns = AnisotropicDiffusivity(νh=8e-3, κh=8e-3, νz=1e-3, κz=1e-3)
model_kwargs = (architecture = arch,
                            grid = grid,
                            advection = UpwindBiasedThirdOrder(),
                            timestepper = :RungeKutta3,
                            coriolis = FPlane(f=f0),
                            tracers = (:b,),
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (b=bbc, u=ubc, v=vbc, w=wbc),
                            forcing = (u=full_sponge, v=full_sponge, w=full_sponge, b=full_sponge),
                            background_fields = bg_fields,
                            )
model = IncompressibleModel(; model_kwargs..., closure=closure_cns)
println("\n", model, "\n")
#-----

# Adding the ICs
#++++
set!(model, u=u_ic, v=v_ic, b=b_ic)

v̄ = sum(model.velocities.v.data.parent) / (grid.Nx * grid.Ny * grid.Nz)
model.velocities.v.data.parent .-= v̄


if false
    using Plots; pyplot()
    xC = oc.Grids.xnodes(Center, grid)
    yC = oc.Grids.ynodes(Center, grid)
    zC = oc.Grids.znodes(Center, grid)
    zF = oc.Grids.znodes(Face, grid)

    contourf(y, z, interior(model.tracers.b)[1,:,:]', levels=30)
end
#-----


# Define time-stepping and printing
#++++
u_scale = abs(u₀)
Δt = min(grid.Δy, grid.Δz) / u_scale
wizard = TimeStepWizard(cfl=0.15,
                        diffusive_cfl=0.5,
                        Δt=Δt, max_change=1.1, min_change=0.2, max_Δt=Inf, min_Δt=0.1seconds)

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
                        stop_time=10*T_inertial,
                        iteration_interval=5, progress=progress,
                        stop_iteration=400,)
#-----


# START DIAGNOSTICS
#++++

# Preamble
#+++++ Preamble
import Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face

u, v, w = model.velocities
b = model.tracers.b
p = sum(model.pressures)

U = model.background_fields.velocities.u
B = model.background_fields.tracers.b

if as_background
    u_tot = u + U
    b_tot = b + B
else
    u_tot = u
    b_tot = b
end
#----


# Start calculation of snapshot variables
#++++
if false
    νₑ = νz = model.diffusivities.νₑ
else
    if :ν in fieldnames(typeof(model.closure))
        νx = νy = νz = model.closure.ν
    else
        νx = νy = model.closure.νx
        νz = model.closure.νz
    end
end

dbdz = ∂z(b_tot)
ω_x = ∂y(w) - ∂z(v)
ω_y = ∂z(u_tot) - ∂x(w)
#ω_z = ∂x(v) - ∂y(u_tot)

wb_res = @at (Center, Center, Center) w*b
#wb_sgs = @at (Center, Center, Center) νz * dbdz

include("diagnostics.jl")

tke = KernelComputedField(Center, Center, Center, kinetic_energy_ccc!, model;
                          field_dependencies=(u, v, w))

if false
    ε = KernelComputedField(Center, Center, Center, isotropic_viscous_dissipation_ccc!, model;
                            field_dependencies=(νₑ, u, v, w))
else
    ε = KernelComputedField(Center, Center, Center, anisotropic_viscous_dissipation_ccc!, model;
                            field_dependencies=(νx, νy, νz, u, v, w))
end

Ri = KernelComputedField(Center, Center, Face, richardson_number_ccf!, model;
                         field_dependencies=(u_tot, v, b_tot, 0, 0, 0))

Ro = KernelComputedField(Face, Face, Center, rossby_number_ffc!, model;
                         field_dependencies=(u_tot, v, 0, 0, f_0))

PV = KernelComputedField(Face, Face, Face, ertel_potential_vorticity_fff!, model;
                         field_dependencies=(u_tot, v, w, b_tot, f_0))

dwpdz = KernelComputedField(Center, Center, Center, vertical_pressure_distribution_ccc!, model;
                            field_dependencies=(w, p), parameters=1027)

SP_y = KernelComputedField(Center, Center, Center, shear_production_y_ccc!, model;
                           field_dependencies=(u, v, w, U))

SP_z = KernelComputedField(Center, Center, Center, shear_production_z_ccc!, model;
                           field_dependencies=(u, v, w, U))

outputs_snap = (u=u,
                v=v,
                w=w,
                b=b,
                p=ComputedField(p),
                wb_res=ComputedField(wb_res),
                #wb_sgs=ComputedField(wb_sgs),
                dwpdz=ComputedField(dwpdz),
                dbdz=ComputedField(dbdz),
                ω_x=ComputedField(ω_x),
                #ω_z=ComputedField(ω_z),
                tke=tke,
                ε=ε,
                Ro=Ro,
                Ri=Ri,
                PV=PV,
                SP_y=SP_y,
                SP_z=SP_z,
                )

if false
    outputs_snap = merge(outputs_snap, (ν_e=νₑ,))
end
if as_background
    outputs_snap = merge(outputs_snap, (u_tot=ComputedField(u_tot),
                                            b_tot=ComputedField(b_tot),))
end
#-----


# Analysis (high def) SNAPSHOTS
#++++
simulation.output_writers[:out_writer] =
    NetCDFOutputWriter(model, outputs_snap,
                       filepath = @sprintf("out.%s.nc", simname),
                       schedule = TimeInterval(6hours),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float64},
                       field_slicer = FieldSlicer(i=1),
                      )
#-----


# Video (low def) SNAPSHOTS
#++++
simulation.output_writers[:vid_writer] =
    NetCDFOutputWriter(model, outputs_snap,
                       filepath = @sprintf("vid.%s.nc", simname),
                       schedule = TimeInterval(90minutes),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float32},
                       field_slicer = FieldSlicer(i=1, with_halos=false),
                      )
#-----

# AVG outputs
#++++
x_average(F) = AveragedField(F, dims=(1,))
outputs_avg = map(x_average, outputs_snap)

simulation.output_writers[:avg_writer] =
    NetCDFOutputWriter(model, outputs_avg,
                       filepath = @sprintf("avg.%s.nc", simname),
                       schedule = AveragedTimeInterval(2minutes; window=1.9minutes, stride=1),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float64},
                       field_slicer = FieldSlicer(j=(grid.Ny÷frac):Int(grid.Ny*(1-1/frac)), with_halos=false),
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








# Define LES model
#+++++
function swap_model(old_model::IncompressibleModel; kwargs...)
    new_model_specification = Dict(property_name => getproperty(old_model, property_name) 
                                   for property_name in propertynames(old_model))

    for (property_name, property) in kwargs
        new_model_specification[property_name] = property
    end
    
    return IncompressibleModel(; new_model_specification...)
end

model = swap_model(model; closure=closure_LES,
                   background_fields=bg_fields)
println("\n", model, "\n")
#-----


# Define time-stepping and printing
#++++
u_scale = abs(u₀)
Δt = simulation.Δt.Δt
wizard = TimeStepWizard(cfl=0.1,
                        diffusive_cfl=0.1,
                        Δt=Δt, max_change=1.05, min_change=0.2, max_Δt=Inf, min_Δt=0.1seconds)

advCFL = oc.Diagnostics.AdvectiveCFL(wizard)
difCFL = oc.Diagnostics.DiffusiveCFL(wizard)
cfl(t) = t<10hours ? 0.1 : 0.2
function progress(sim)
    wizard.cfl = cfl(sim.model.clock.time)
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


# REDFINE Simulation!
#++++
simulation = Simulation(model, Δt=wizard, 
                        stop_time=10*T_inertial,
                        iteration_interval=5, progress=progress,
                        stop_iteration=Inf,)
#-----


# RESTART DIAGNOSTICS
#++++

# Preamble
#+++++ Preamble
import Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face

u, v, w = model.velocities
b = model.tracers.b
p = sum(model.pressures)

U = model.background_fields.velocities.u
B = model.background_fields.tracers.b

if as_background
    u_tot = u + U
    b_tot = b + B
else
    u_tot = u
    b_tot = b
end
#----


# Start calculation of snapshot variables
#++++
if true
    νₑ = νz = model.diffusivities.νₑ
else
    if :ν in fieldnames(typeof(model.closure))
        νx = νy = νz = model.closure.ν
    else
        νx = νy = model.closure.νx
        νz = model.closure.νz
    end
end

dbdz = ∂z(b_tot)
ω_x = ∂y(w) - ∂z(v)
ω_y = ∂z(u_tot) - ∂x(w)
ω_z = ∂x(v) - ∂y(u_tot)

wb_res = @at (Center, Center, Center) w*b
#wb_sgs = @at (Center, Center, Center) νz * dbdz

include("diagnostics.jl")

tke = KernelComputedField(Center, Center, Center, kinetic_energy_ccc!, model;
                          field_dependencies=(u, v, w))

if true
    ε = KernelComputedField(Center, Center, Center, isotropic_viscous_dissipation_ccc!, model;
                            field_dependencies=(νₑ, u, v, w))
else
    ε = KernelComputedField(Center, Center, Center, anisotropic_viscous_dissipation_ccc!, model;
                            field_dependencies=(νx, νy, νz, u, v, w))
end

Ri = KernelComputedField(Center, Center, Face, richardson_number_ccf!, model;
                         field_dependencies=(u_tot, v, b_tot, 0, 0, 0))

Ro = KernelComputedField(Face, Face, Center, rossby_number_ffc!, model;
                         field_dependencies=(u_tot, v, 0, 0, f_0))

PV = KernelComputedField(Face, Face, Face, ertel_potential_vorticity_fff!, model;
                         field_dependencies=(u_tot, v, w, b_tot, f_0))

dwpdz = KernelComputedField(Center, Center, Center, vertical_pressure_distribution_ccc!, model;
                            field_dependencies=(w, p), parameters=1027)

SP_y = KernelComputedField(Center, Center, Center, shear_production_y_ccc!, model;
                           field_dependencies=(u, v, w, U))

SP_z = KernelComputedField(Center, Center, Center, shear_production_z_ccc!, model;
                           field_dependencies=(u, v, w, U))

outputs_snap = (u=u,
                  v=v,
                  w=w,
                  b=b,
                  p=ComputedField(p),
                  wb_res=ComputedField(wb_res),
                  #wb_sgs=ComputedField(wb_sgs),
                  dwpdz=ComputedField(dwpdz),
                  dbdz=ComputedField(dbdz),
                  ω_x=ComputedField(ω_x),
                  #ω_z=ComputedField(ω_z),
                  tke=tke,
                  ε=ε,
                  Ro=Ro,
                  Ri=Ri,
                  PV=PV,
                  SP_y=SP_y,
                  SP_z=SP_z,
                  )

if true
    outputs_snap = merge(outputs_snap, (ν_e=νₑ,))
end
if as_background
    outputs_snap = merge(outputs_snap, (u_tot=ComputedField(u_tot),
                                            b_tot=ComputedField(b_tot),))
end
#-----


# Analysis (high def) SNAPSHOTS
#++++
simulation.output_writers[:out_writer] =
    NetCDFOutputWriter(model, outputs_snap,
                       filepath = @sprintf("out.%s.nc", simname),
                       schedule = TimeInterval(6hours),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float64},
                       field_slicer = FieldSlicer(i=1),
                      )
#-----


# Video (low def) SNAPSHOTS
#++++
simulation.output_writers[:vid_writer] =
    NetCDFOutputWriter(model, outputs_snap,
                       filepath = @sprintf("vid.%s.nc", simname),
                       schedule = TimeInterval(90minutes),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float32},
                       field_slicer = FieldSlicer(i=1, with_halos=false),
                      )
#-----

# AVG outputs
#++++
outputs_avg = map(x_average, outputs_snap)
simulation.output_writers[:avg_writer] =
    NetCDFOutputWriter(model, outputs_avg,
                       filepath = @sprintf("avg.%s.nc", simname),
                       schedule = AveragedTimeInterval(5minutes; window=4.9minutes, stride=1),
                       mode = "c",
                       global_attributes = global_attributes,
                       array_type = Array{Float64},
                       field_slicer = FieldSlicer(j=(grid.Ny÷frac):Int(grid.Ny*(1-1/frac)), with_halos=false),
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
