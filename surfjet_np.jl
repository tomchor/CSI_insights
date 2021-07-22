using Pkg
Pkg.instantiate()
using ArgParse
using Printf
using Oceananigans
@info "Loaded Oceananigans"
using Oceananigans.Units
using Oceananigans.Advection: WENO5
using Oceananigans.OutputWriters, Oceananigans.Fields
using SpecialFunctions: erf
using CUDA: has_cuda


# Read and parse initial arguments
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
            default = "CPU"
            arg_type = String

        "--fullname"
            help = "Setup and name of jet in jetinfo.jl"
            default = "S2d_SIjet1"
            arg_type = String
    end
    return parse_args(settings)
end
args = parse_command_line_arguments()
factor = args["factor"]
arch = eval(Meta.parse(args["arch"]*"()"))
fullname = args["fullname"]


try
    global setup, jet, modifier = split(fullname, "_")
    global AMD = modifier=="AMD" ? true : false
    global noflux = modifier=="NF" ? true : false
catch e
    global setup, jet = split(fullname, "_")
    global AMD = false
    global noflux = false
end
ndims = parse(Int, strip(setup, ['S', 'd']))
jet = Symbol(jet)

if ndims ∉ [2,3] # Validade input
    throw(AssertionError("Dimensions must be 2 or 3"))
end

@printf("Starting %sd jet %s with a dividing factor of %d and a %s architecture\n", 
        ndims,
        jet,
        factor,
        arch,)
#-----


# Get simulation parameters
#++++
as_background=false
include("jetinfo.jl")

if ndims==3 # 3D LES simulation
    simulation_nml = getproperty(SurfaceJetSimulations(Ny=400*2^4, Nz=2^7, ThreeD=true), jet)
    prefix = "PNN"
    LES = true
else # 2D DNS simulation
    simulation_nml = getproperty(SurfaceJetSimulations(), jet)
    prefix = "FNN"
    LES = false
end
@unpack name, f0, u₀, N2_inf, N2_pyc, Ny, Nz, Ly, Lz, σy, σz, y₀, z₀, νz, sponge_frac = simulation_nml

if @isdefined modifier
    simname = "$(prefix)_$(name)_$(modifier)"
else
    simname = "$(prefix)_$(name)"
end
pickup = any(startswith("chk.$simname"), readdir("data"))
#-----

# Set GRID
#++++  GRID
if ndims==3
    Nx = Ny÷64
    Lx = (Ly / Ny) * Nx
else
    Nx = factor
    Lx = 6 * (Ly / Ny) * Nx
end
topology = (Periodic, Bounded, Bounded)

grid = RegularRectilinearGrid(size=(Nx÷factor, Ny÷factor, Nz÷factor),
                              x=(0, Lx),
                              y=(0, Ly),
                              z=(-Lz, 0), 
                              topology=topology)
println("\n", grid, "\n")
#-----


# Calculate secondary parameters
#++++
b₀ = u₀ * f0
ρ₀ = 1027
T_inertial = 2*π/f0
y_r = y₀ + √2/4 * σy
z_r = 0
Ro_r = - √2 * u₀ * (z₀/σz-1) * exp(-1/8) / (2*f0*σy)
Ri_r = N2_inf * σz^2 * exp(1/4) / u₀^2
νh = νz * (grid.Δy / grid.Δz)^(4/3)

secondary_params = merge((LES=Int(LES), u_0=u₀, y_0=y₀, z_0=z₀, b0=b₀), 
                         (;y_r, z_r, Ro_r, Ri_r, T_inertial, νh))

global_attributes = merge(simulation_nml, secondary_params)
println("\n", global_attributes, "\n")
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
@inline fy(ψ) = exp(-ψ^2)
@inline intgaussian(ψ) = √π/2 * (erf(ψ) + 1)
@inline umask(Y, Z) = Z * fy(Y)
@inline bmask(Y, Z) = (1/sig_z) * (sig_y * intgaussian(Y))

u_g(x, y, z, t) = +u_0 * umask((y-y_0)/sig_y, ((z-z_0)/sig_z +1))
@inline background_strat(z) = ifelse(z < z_c, 
                             n2_pyc * (z+Hz),
                             n2_inf * (z-z_m))
b_g(x, y, z, t) = -f_0 * u_0 * bmask((y-y_0)/sig_y, ((z-z_0)/sig_z +1)) + background_strat(z)
@inline dudz_g(x, y, z, t) = +u_0 * (1/sig_z) * fy((y-y_0)/sig_y)
#-----

# Setting BCs
#++++
if as_background
    throw(ArgumentError("background isn't used anymore!"))
else
    U_top_bc = FluxBoundaryCondition(0)
    U_bot_bc = FluxBoundaryCondition(0)
    if noflux
        B_bc = FluxBoundaryCondition(0)
    else
        B_bc = GradientBoundaryCondition(N2_inf)
    end
end

ubc = FieldBoundaryConditions(top = U_top_bc,
                              bottom = U_bot_bc,
                              )
vbc = FieldBoundaryConditions(top = FluxBoundaryCondition(0),
                              bottom = FluxBoundaryCondition(0),
                               )
wbc = FieldBoundaryConditions()
bbc = FieldBoundaryConditions(bottom = B_bc,
                              top = B_bc,
                              )
#-----


# Set-up sponge layer
#++++
@inline heaviside(X) = ifelse(X < 0, zero(X), one(X))
@inline mask2nd(X) = heaviside(X) * X^2
const frac = sponge_frac

@inline function north_mask(x, y, z)
    y₁ = Hy; y₀ = y₁ - Hy*frac
    return mask2nd((y - y₀)/(y₁ - y₀))
end
@inline function south_mask(x, y, z)
    y₁ = 0; y₀ = y₁ + Hy*frac
    return mask2nd((y - y₀)/(y₁ - y₀))
end
@inline full_mask(x, y, z) = north_mask(x, y, z) + south_mask(x, y, z)

const rate = 1/10minutes
full_sponge_0 = Relaxation(rate=rate, mask=full_mask, target=0)
if as_background
    throw(ArgumentError("background isn't used anymore!"))
else
    full_sponge_u = Relaxation(rate=rate, mask=full_mask, target=u_g)
    full_sponge_b = Relaxation(rate=rate, mask=full_mask, target=b_g)
    forcing = (u=full_sponge_u, v=full_sponge_0, w=full_sponge_0)
end
#-----



# Set up ICs and/or Background Fields
#++++
const amplitude = 1e-6
if as_background
    throw(ArgumentError("background isn't used anymore!"))
else
    println("\nSetting geostrophic jet as an INITIAL CONDITION\n")
    u_ic(x, y, z) = u_g(x, y, z, 0) + amplitude*randn()
    v_ic(x, y, z) = + amplitude*randn()
    w_ic(x, y, z) = + amplitude*randn()
    b_ic(x, y, z) = b_g(x, y, z, 0) #+ 1e-8*randn()

    bg_fields = NamedTuple()
end
#-----


# Define model!
#++++
if LES
    import Oceananigans.TurbulenceClosures: SmagorinskyLilly, AnisotropicMinimumDissipation
    νₘ, κₘ = 1.0e-6, 1.5e-7
    if AMD
        closure = AnisotropicMinimumDissipation(ν=νₘ, κ=κₘ)
    else
        closure = SmagorinskyLilly(C=0.16, ν=νₘ, κ=κₘ)
    end
else
    import Oceananigans.TurbulenceClosures: AnisotropicDiffusivity, IsotropicDiffusivity
    closure = AnisotropicDiffusivity(νh=νh, κh=νh, νz=νz, κz=νz)
end
model_kwargs = (architecture = arch,
                grid = grid,
                advection = WENO5(),
                timestepper = :RungeKutta3,
                coriolis = occursin("slosh", lowercase(simname)) ? FPlane(f=0) : FPlane(f=f0), # do we want sloshing?
                tracers = (:b,),
                buoyancy = BuoyancyTracer(),
                boundary_conditions = (b=bbc, u=ubc, v=vbc, w=wbc),
                forcing = forcing,
                background_fields = bg_fields,
                )
model = NonhydrostaticModel(; model_kwargs..., closure=closure)
println("\n", model, "\n")
#-----


# Adding the ICs
#++++
set!(model, u=u_ic, v=v_ic, w=w_ic, b=b_ic)

v̄ = sum(model.velocities.v.data.parent) / (grid.Nx * grid.Ny * grid.Nz)
model.velocities.v.data.parent .-= v̄
#-----


# Define time-stepping
#++++
u_scale = abs(u₀)
Δt = 1/2 * min(grid.Δz, grid.Δy) / u_scale
wizard = TimeStepWizard(cfl=0.9,
                        diffusive_cfl=0.9,
                        Δt=Δt, max_change=1.02, min_change=0.2, max_Δt=Inf, min_Δt=0.1seconds)
#----


# Finally define Simulation!
#++++
if ndims==3 # 3D LES simulation
    stop_time = min(10*T_inertial, 20days)
else # 2D DNS simulation
    stop_time = min(3*T_inertial, 20days)
end
include("diagnostics.jl")
start_time = 1e-9*time_ns()
using Oceanostics: SingleLineProgressMessenger
simulation = Simulation(model, Δt=wizard, 
                        stop_time=stop_time,
                        wall_time_limit=23.5hours,
                        iteration_interval=5,
                        progress=SingleLineProgressMessenger(LES=LES, initial_wall_time_seconds=start_time),
                        stop_iteration=Inf,)
println("\n", simulation, "\n")
#-----


# DIAGNOSTICS
#++++
const ρ0 = ρ₀
checkpointer = construct_outputs(model, simulation, LES=LES, simname=simname, frac=frac)
#-----


# Run the simulation!
#+++++
println("\n", simulation,
        "\n",)
if has_cuda() run(`nvidia-smi`) end


@printf("---> Starting run!\n")
run!(simulation, pickup=true)

using Oceananigans.OutputWriters: write_output!
write_output!(checkpointer, model)
#-----
