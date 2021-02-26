using Oceananigans
using Oceananigans.Utils
using Oceananigans.Fields

Lx = 150; Ly = 6000; Lz = 80
topology = (Periodic, Bounded, Bounded)
grid = RegularRectilinearGrid(size=(1, 512, 8), x=(0, Lx), y=(0, Ly), z=(-Lz, 0), 
                              topology=(Periodic, Bounded, Bounded))


model = IncompressibleModel(architecture = GPU(),
                            grid = grid,
                            )

w_ic(x, y, z) = 0.01*y
v_ic(x, y, z) = 0.01*x
set!(model, w=w_ic, v=v_ic)


import Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Grids: Center, Face


u, v, w = model.velocities

function naive_calc()
    p = sum(model.pressures)
    wp = @at (Center, Center, Face) w*p
    dwpdz = (1/1024) * ∂z(wp)
    println(dwpdz)
    return ComputedField(dwpdz^2 + dwpdz)
end

function nested_calc()
    p = ComputedField(sum(model.pressures))
    wp = ComputedField(@at (Center, Center, Face) w*p)
    dwpdz = (1/1024) * ∂z(wp)
    println(dwpdz)
    return ComputedField(dwpdz)
end

function crazy_calc()
    p = ComputedField(sum(model.pressures))
    wp = ComputedField(@at (Center, Center, Face) w*p)
    dwpdz = (1/1024) * ∂z(wp)
    println(dwpdz)
    dwpdz = ComputedField(dwpdz)
    dwpdz2 = ComputedField(dwpdz^2)
    return ComputedField(dwpdz2+dwpdz)
end

function naive_tke()
    return ComputedField((u^2+v^2+w^2)/2)
end

function nested_tke()
    u2 = ComputedField(u^2)
    v2 = ComputedField(v^2)
    w2 = ComputedField(w^2)
    return ComputedField((u2+v2+w2)/2)
end
#dwpdz_naive = naive_calc()
#dwpdz_nested = nested_calc()


pause
using KernelAbstractions: @index, @kernel
using Oceananigans.Operators
include("Operators.jl")
#using Operators: ∂zᶜ!

p = ComputedField(sum(model.pressures))
vw = @at (Center, Center, Face) w*v
wp = ComputedField(@at (Center, Center, Face) w*p)


dwpdz = KernelComputedField(Center, Center, Center, ϕ∂zᶜ!, model;
                            field_dependencies=(wp,), parameters=1/1024)
compute!(dwpdz)


dwpdz_2 = ComputedField(∂z(w*p)/1024)
compute!(dwpdz_2)

dwpdz_3 = ComputedField(dwpdz_2^3)
compute!(dwpdz_3)



function nested_comp_calc()
    p = ComputedField(sum(model.pressures))
    compute!(p)
    wp = ComputedField(@at (Center, Center, Face) w*p)
    compute!(wp)
    dwpdz = (1/1024) * ∂z(wp)
    println(dwpdz)
    return ComputedField(dwpdz)
end

dwpdz_nestedcomp = nested_comp_calc()

pause

dvwdz = KernelComputedField(Center, Center, Center, ∂zᶜ!, model;
                            field_dependencies=(ww,))
compute!(dvwdz)


dvwdz_2 = ComputedField(∂z(vw))
compute!(dvwdz_2)

