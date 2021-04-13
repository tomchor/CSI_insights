using Printf
using KernelAbstractions: @index, @kernel

using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.Diagnostics: WindowedSpatialAverage
using Oceananigans.Grids: Center, Face

using Oceanostics.FlowDiagnostics: IsotropicBuoyancyMixingRate, AnisotropicBuoyancyMixingRate
using Oceanostics.TurbulentKineticEnergyTerms: KineticEnergy, 
                                               IsotropicViscousDissipationRate, IsotropicPseudoViscousDissipationRate,
                                               AnisotropicViscousDissipationRate, AnisotropicPseudoViscousDissipationRate,
                                               PressureRedistribution_y, PressureRedistribution_z


#++++ KERNEL COMPUTED FIELDS

#++++ Shear production terms
@inline ψ²(i, j, k, grid, ψ) = @inbounds ψ[i, j, k]^2
@inline ψ′²(i, j, k, grid, ψ, Ψ) = @inbounds (ψ[i, j, k] - Ψ[i, j, k])^2

@kernel function shear_production_y_ccc!(shear_production, grid, u, v, w, U)
    i, j, k = @index(Global, NTuple)
    v_int = ℑyᵃᶜᵃ(i, j, k, grid, v) # C, F, C  → C, C, C

    ∂yU = ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᵃᶠᵃ, U) # F, C, C  → F, F, C  → C, C, C
    uv = ℑxᶜᵃᵃ(i, j, k, grid, u) * v_int
    uv∂yU = uv * ∂yU

    @inbounds shear_production[i, j, k] = -uv∂yU
end



@kernel function shear_production_z_ccc!(shear_production, grid, u, v, w, U)
    i, j, k = @index(Global, NTuple)
    w_int = ℑzᵃᵃᶜ(i, j, k, grid, w) # C, C, F  → C, C, C

    ∂zU = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᵃᵃᶠ, U) # F, C, C  → F, C, F  → C, C, C
    uw = ℑxᶜᵃᵃ(i, j, k, grid, u) * w_int
    uw∂zU = uw * ∂zU

    @inbounds shear_production[i, j, k] = -uw∂zU
end
#-----



#++++ PV components
@kernel function ertel_potential_vorticity_vertical_fff!(PV, grid, u, v, b, f₀)
    i, j, k = @index(Global, NTuple)

    dVdx =  ℑzᵃᵃᶠ(i, j, k, grid, ∂xᶠᵃᵃ, v) # C, F, C  → F, F, C → F, F, F
    dUdy =  ℑzᵃᵃᶠ(i, j, k, grid, ∂yᵃᶠᵃ, u) # F, C, C  → F, F, C → F, F, F
    dbdz = ℑxyᶠᶠᵃ(i, j, k, grid, ∂zᵃᵃᶠ, b) # C, C, C  → C, C, F → F, F, F
    pv_z = (f₀ + dVdx - dUdy) * dbdz

    @inbounds PV[i, j, k] = pv_z
end



@kernel function ertel_potential_vorticity_horizontal_fff!(PV, grid, u, v, w, b, f₀)
    i, j, k = @index(Global, NTuple)

    dWdy =  ℑxᶠᵃᵃ(i, j, k, grid, ∂yᵃᶠᵃ, w) # C, C, F  → C, F, F  → F, F, F
    dVdz =  ℑxᶠᵃᵃ(i, j, k, grid, ∂zᵃᵃᶠ, v) # C, F, C  → C, F, F  → F, F, F
    dbdx = ℑyzᵃᶠᶠ(i, j, k, grid, ∂xᶠᵃᵃ, b) # C, C, C  → F, C, C  → F, F, F
    pv_x = (dWdy - dVdz) * dbdx # F, F, F

    dUdz =  ℑyᵃᶠᵃ(i, j, k, grid, ∂zᵃᵃᶠ, u) # F, C, C  → F, C, F → F, F, F
    dWdx =  ℑyᵃᶠᵃ(i, j, k, grid, ∂xᶠᵃᵃ, w) # C, C, F  → F, C, F → F, F, F
    dbdy = ℑxzᶠᵃᶠ(i, j, k, grid, ∂yᵃᶠᵃ, b) # C, C, C  → C, F, C → F, F, F
    pv_y = (dUdz - dWdx) * dbdy # F, F, F

    @inbounds PV[i, j, k] = pv_x + pv_y
end
#----

#-----


#++++ CREATE SNAPSHOT OUTPUTS
function get_outputs_tuple(model; LES=false)
    
    #++++ Preamble
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
    
    
    if LES
        νₑ = νz = model.diffusivities.νₑ
        κₑ = κz = model.diffusivities.νₑ # Assumes Pr=1
    else
        if model.closure isa IsotropicDiffusivity
            νx = νy = νz = model.closure.ν
            κx = κy = κz = model.closure.κ[:b]
        elseif model.closure isa AnisotropicDiffusivity
            νx = model.closure.νx; νy = model.closure.νy; νz = model.closure.νz
            κx = model.closure.κx[:b]; κy = model.closure.κy[:b]; κz = model.closure.κz[:b]
        end
    end
    #----
    
    #++++ Create scratch space
    ccc_scratch = Field(Center, Center, Center, model.architecture, model.grid)
    fcc_scratch = Field(Face, Center, Center, model.architecture, model.grid)
    ccf_scratch = Field(Center, Center, Face, model.architecture, model.grid)
    cff_scratch = Field(Center, Face, Face, model.architecture, model.grid)
    fff_scratch = Field(Face, Face, Face, model.architecture, model.grid)
    #----

    # Start calculation of snapshot variables
    #++++
    dbdz = @at (Center, Center, Face) ∂z(b_tot)
    ω_x = ∂y(w) - ∂z(v)
    
    wb_res = @at (Center, Center, Center) w*b
    tke = KineticEnergy(model, u, v, w, data=ccc_scratch.data)
    
    if LES
        ε = IsotropicViscousDissipationRate(model, u, v, w, νₑ, data=ccc_scratch.data)
        ε2 = IsotropicPseudoViscousDissipationRate(model, u, v, w, νₑ, data=ccc_scratch.data)
        χ = IsotropicBuoyancyMixingRate(model, b, κₑ, n2_inf, data=ccc_scratch.data)
    else
        ε = AnisotropicViscousDissipationRate(model, u, v, w, νx, νy, νz, data=ccc_scratch.data)
        ε2 = AnisotropicPseudoViscousDissipationRate(model, u, v, w, νx, νy, νz, data=ccc_scratch.data)
        χ = AnisotropicBuoyancyMixingRate(model, b, κx, κy, κz, n2_inf, data=ccc_scratch.data)
    end

    PV_ver = KernelComputedField(Face, Face, Face, ertel_potential_vorticity_vertical_fff!, model;
                                 computed_dependencies=(u_tot, v, b_tot), 
                                 parameters=f_0, data=fff_scratch.data)
    
    PV_hor = KernelComputedField(Face, Face, Face, ertel_potential_vorticity_horizontal_fff!, model;
                                 computed_dependencies=(u_tot, v, w, b_tot), 
                                 parameters=f_0, data=fff_scratch.data)
    
    dvpdy_ρ = PressureRedistribution_y(model, v, p, ρ₀, data=ccc_scratch.data)
    dwpdz_ρ = PressureRedistribution_z(model, w, p, ρ₀, data=ccc_scratch.data)

    
    SP_y = KernelComputedField(Center, Center, Center, shear_production_y_ccc!, model;
                               computed_dependencies=(u, v, w, U),
                               data=ccc_scratch.data,
                              )
    
    SP_z = KernelComputedField(Center, Center, Center, shear_production_z_ccc!, model;
                               computed_dependencies=(u, v, w, U),
                               data=ccc_scratch.data,
                              )
    #-----
    
    
    # Assemble the outputs tuple
    #++++
    outputs = (u=u,
               v=v,
               w=w,
               b=b,
               p=ComputedField(p, data=ccc_scratch.data),
               wb_res=ComputedField(wb_res, data=ccc_scratch.data),
               dwpdz_ρ=dwpdz_ρ,
               dvpdy_ρ=dvpdy_ρ,
               dbdz=ComputedField(dbdz, data=ccf_scratch.data),
               ω_x=ComputedField(ω_x, data=cff_scratch.data),
               tke=tke,
               ε=ε,
               ε2=ε2,
               χ=χ,
               PV_ver=PV_ver,
               PV_hor=PV_hor,
               SP_y=SP_y,
               SP_z=SP_z,
               )
    
    if LES
        outputs = merge(outputs, (ν_e=νₑ,))
    end
    if as_background
        outputs = merge(outputs, (u_tot=ComputedField(u_tot, data=fcc_scratch.data),
                                  b_tot=ComputedField(b_tot, data=ccc_scratch.data),))
    end
    #----

    return outputs
end
#----


#++++ Construct outputs into simulation
function construct_outputs(model, simulation; 
                           LES=false, simname="TEST", frac=1/8,
                           )

    if any(startswith("chk.$simname"), readdir("data"))
        @info "Checkpoint for $simname found. Assuming this is a pick-up simulation! Setting mode to append."
        mode = "a"
    else
        @info "No checkpoint for $simname found. Setting mode to clobber."
        mode = "c"
    end
    
    # Output (high def) SNAPSHOTS
    #++++
    @info "Setting up out writer"
    outputs_snap = get_outputs_tuple(model, LES=LES)
    simulation.output_writers[:out_writer] =
        NetCDFOutputWriter(model, outputs_snap,
                           filepath = @sprintf("data/out.%s.nc", simname),
                           schedule = TimeInterval(6hours),
                           mode = mode,
                           global_attributes = global_attributes,
                           array_type = Array{Float64},
                          )
    #-----
    
    
    # Video (low def) SNAPSHOTS
    #++++
    @info "Setting up vid writer"
    delete(nt::NamedTuple{names}, keys) where names = NamedTuple{filter(x -> x ∉ keys, names)}(nt)
    
    outputs_vid = delete(outputs_snap, (:SP_y, :SP_z, :dwpdz_ρ, :dvpdy_ρ, :p))
    
    simulation.output_writers[:vid_writer] =
        NetCDFOutputWriter(model, outputs_vid,
                           filepath = @sprintf("data/vid.%s.nc", simname),
                           schedule = TimeInterval(90minutes),
                           mode = mode,
                           global_attributes = global_attributes,
                           array_type = Array{Float32},
                           field_slicer = FieldSlicer(i=1, with_halos=false),
                           )
    #-----
    
    
    # AVG outputs
    #++++
    @info "Setting up avg writer"
    x_average(F) = AveragedField(F, dims=(1,))
    xz_average(F) = AveragedField(F, dims=(1,3))
    
    slicer = FieldSlicer(j=Int(grid.Ny*frac):Int(grid.Ny*(1-1*frac)), with_halos=false)
    hor_window_average(F) = WindowedSpatialAverage(F; dims=(1, 2), field_slicer=slicer)
    
    outputs_avg = map(hor_window_average, outputs_snap)
    
    simulation.output_writers[:avg_writer] =
        NetCDFOutputWriter(model, outputs_avg,
                           filepath = @sprintf("data/avg.%s.nc", simname),
                           schedule = TimeInterval(2minutes),
                           mode = mode,
                           global_attributes = global_attributes,
                           array_type = Array{Float64},
                          )
    #-----


    # Checkpointer
    #+++++
    @info "Setting up chk writer"
    simulation.output_writers[:chk_writer] = checkpointer = 
                                             Checkpointer(model;
                                             dir="data/",
                                             prefix = @sprintf("chk.%s", simname),
                                             schedule = TimeInterval(18hours),
                                             force = true,
                                             cleanup = true,
                                             )
    #-----

    return checkpointer
end
#-----




