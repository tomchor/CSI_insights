using ArgParse
using Printf
using Oceananigans; oc = Oceananigans
using Oceananigans.Utils
using Oceananigans.Units
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
            default = 16
            arg_type = Int

        "--arch"
            help = "CPU or GPU"
            default = "CPU"
            arg_type = String

        "--jet"
            help = "Name of jet in jetinfo.jl"
            default = :CIjet01
            arg_type = Symbol
    end
    return parse_args(settings)
end
args = parse_command_line_arguments()
factor = args["factor"]
arch = eval(Meta.parse(args["arch"]*"()"))
jet = args["jet"]


@printf("Starting Oc with jet %s, a dividing factor of %d and a %s architecture\n", 
        args["jet"],
        factor,
        arch,)
#-----


# Get simulation parameters
#++++
LES = false
as_background=true
include("jetinfo.jl")

simulation_nml = getproperty(InteriorJetSimulations(), jet)
@unpack name, f0, u₀, N2_inf, N2_pyc, Ny, Nz, Ly, Lz, σy, σz, y₀, z₀, νh, νz = simulation_nml

simname = @sprintf("FNN_%s", name)
#-----


# Calculate secondary parameters
#++++
b₀ = u₀ * f0
ρ₀ = 1027
T_inertial = 2*π/f0

global_attributes = merge(simulation_nml,
                          (LES=Int(LES),
                          u_0=u₀, y_0=y₀, z_0=z₀,
                          b0 = b₀, T_inertial = T_inertial,),
                         )
println("\n", global_attributes, "\n")
#-----

# Set GRID
#++++  GRID
Nx = 1
Lx = 6 * (Ly / Ny) * Nx
topology = (Periodic, Bounded, Bounded)
grid = RegularRectilinearGrid(size=(Nx, Ny÷factor, Nz÷factor),
                            x=(0, Lx),
                            y=(0, Ly),
                            z=(-Lz, 0), 
                            topology=topology)
println("\n", grid, "\n")
#-----


# Set up Geostrophic flow
#++++++
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

gaussian(ψ) = exp(-ψ^2)
intgaussian(ψ) = √π/2 * (erf(ψ) + 1)
umask(Y, Z) = gaussian(Y) * gaussian(Z)
bmask(Y, Z) = (sig_y * intgaussian(Y)) * (-2Z * gaussian(Z) / sig_z)

u_g(x, y, z, t) = u_0 * umask((y-y_0)/sig_y, (z-z_0)/sig_z)
b_g(x, y, z, t) = -u_0*f_0 * bmask((y-y_0)/sig_y, (z-z_0)/sig_z) + n2_inf * (z+Hz/2)
dudz_g(x, y, z, t) = + u_0 * gaussian((y-y_0)/sig_y) * (-2(z-z_0)/sig_z^2) * gaussian((z-z_0)/sig_z)
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
                                  bottom = U_bot_bc,
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
@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
@inline mask2nd(X) = heaviside(X) * X^2
const Hy = grid.Ly
const frac = 8

@inline function bottom_mask(x, y, z)
    z₁ = -Hz; z₀ = z₁ + Hz/frac
    return mask2nd((z - z₀)/(z₁ - z₀))
end
@inline function top_mask(x, y, z)
    z₁ = +Hz; z₀ = z₁ - Hz/frac
    return mask2nd((z - z₀)/(z₁ - z₀))
end
@inline function north_mask(x, y, z)
    y₁ = Hy; y₀ = y₁ - Hy/frac
    return mask2nd((y - y₀)/(y₁ - y₀))
end
@inline function south_mask(x, y, z)
    y₁ = 0; y₀ = y₁ + Hy/frac
    return mask2nd((y - y₀)/(y₁ - y₀))
end

full_mask(x, y, z) = north_mask(x, y, z) + south_mask(x, y, z)# + bottom_mask(x, y, z)
rate = 1/5hours
if as_background
    full_sponge_0 = Relaxation(rate=rate, mask=full_mask, target=0)
    forcing = (u=full_sponge_0, v=full_sponge_0, w=full_sponge_0, b=full_sponge_0)
else
    full_sponge_0 = Relaxation(rate=rate, mask=full_mask, target=0)
    full_sponge_u = Relaxation(rate=rate, mask=full_mask, target=u_g)
    full_sponge_b = Relaxation(rate=rate, mask=full_mask, target=b_g)
    forcing = (u=full_sponge_u, v=full_sponge_0, w=full_sponge_0, b=full_sponge_b)
end
#-----



# Set up ICs and/or Background Fields
#++++
kick = 1e-6
if as_background
    println("\nSetting geostrophic jet as BACKGROUND\n")
    u_ic(x, y, z) = + kick*randn()
    v_ic(x, y, z) = + kick*randn()
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
if LES
    import Oceananigans.TurbulenceClosures: SmagorinskyLilly, AnisotropicMinimumDissipation
    closure = SmagorinskyLilly(C=0.23)
else
    import Oceananigans.TurbulenceClosures: AnisotropicDiffusivity, IsotropicDiffusivity
    closure = AnisotropicDiffusivity(νh=νh, κh=νh, νz=νz, κz=νz)
end
model = IncompressibleModel(architecture = arch,
                            grid = grid,
                            advection = UpwindBiasedThirdOrder(),
                            timestepper = :RungeKutta3,
                            closure = closure,
                            coriolis = FPlane(f=f0),
                            tracers = (:b,),
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (b=bbc, u=ubc, v=vbc, w=wbc),
                            forcing = forcing,
                            background_fields = bg_fields,
                            )
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
Δt = 0.1 * min(grid.Δx, grid.Δy) / u_scale
wizard = TimeStepWizard(cfl=0.15,
                        diffusive_cfl=0.6,
                        Δt=Δt, max_change=1.01, min_change=0.02, max_Δt=Inf, min_Δt=0.2seconds)
#-----

# Finally define Simulation!
#++++
start_time = 1e-9*time_ns()
using Oceanostics: ProgressMessenger
simulation = Simulation(model, Δt=wizard, 
                        stop_time=20*T_inertial,
                        iteration_interval=5, progress=ProgressMessenger(LES=LES,
                                                                         initial_wall_time_seconds=start_time),
                        stop_iteration=Inf,)
#-----

# START DIAGNOSTICS
#++++
import Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Diagnostics: WindowedSpatialAverage
using Oceanostics.FlowDiagnostics: richardson_number_ccf!, rossby_number_ffc!, ertel_potential_vorticity_fff!
using Oceanostics.TurbulentKineticEnergyTerms: KineticEnergy, 
                                               IsotropicViscousDissipation, AnisotropicViscousDissipation,
                                               PressureRedistribution_y, PressureRedistribution_z

const ρ0 = ρ₀
include("diagnostics.jl")
construct_outputs(model, simulation, LES=LES, simname=simname)
#-----


# Run the simulation!
#+++++
println("\n", simulation,
        "\n",)

@printf("---> Starting run!\n")
run!(simulation)
#-----

