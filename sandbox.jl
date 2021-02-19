using Oceananigans; using Oceananigans.Utils

# Get simulation parameters
#++++
const f0 = 1e-4
const dUdy = 2e-2
const dUdz = -1e-2
const dbdy = 1e-6
const dbdz = 3e-6
#-----



# Set GRID
#+++++
topology = (Periodic, Periodic, Periodic)
grid = RegularCartesianGrid(size=(16, 16, 4), extent=(1000, 1000, 200), 
                            topology=topology)
#-----


u_g(x, y, z, t) = dUdy * y + dUdz * z
b_g(x, y, z, t) = dbdy * y + dbdz * z


model = IncompressibleModel(architecture = CPU(), grid = grid,
                            coriolis = FPlane(f=f0),
                            tracers = (:b,),
                            buoyancy = BuoyancyTracer(),
                            background_fields = (u=u_g, b=b_g,)
                            )
println("\n", model, "\n")


u_ic(x, y, z) = 2
v_ic(x, y, z) = 3
w_ic(x, y, z) = 5
b_ic(x, y, z) = 7

set!(model, u=u_ic, b=b_ic, v=v_ic, w=w_ic)
#-----


# START DIAGNOSTICS
#++++
import Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face
using Oceananigans.Fields: compute!

u, v, w = model.velocities
b = model.tracers.b
p = sum(model.pressures)

U = model.background_fields.velocities.u
B = model.background_fields.tracers.b

u_tot = u + U
b_tot = b + B


#include("/glade/u/home/tomasc/repos/KCF_diagnostics/tke_budget_terms.jl")
include("diagnostics.jl")



pause
SPy_true = -u_ic(0, 0, 0)*v_ic(0, 0, 0) * dUdy
SPy_kern = KernelComputedField(Center, Center, Center, shear_production_y_ccc!, model;
                           field_dependencies=(u, v, w, U))
compute!(SPy_true)
compute!(SPy_kern)


SPz_true = -u_ic(0, 0, 0)*w_ic(0, 0, 0) * dUdz
SPz_kern = KernelComputedField(Center, Center, Center, shear_production_z_ccc!, model;
                           field_dependencies=(u, v, w, U))
compute!(SPz_true)
compute!(SPz_kern)

tke_true = ComputedField(1/2*(u^2 + v^2 + w^2))
tke_kern = KernelComputedField(Center, Center, Center, kinetic_energy_ccc!, model;
                           field_dependencies=(u, v, w))
compute!(tke_true)
compute!(tke_kern)
