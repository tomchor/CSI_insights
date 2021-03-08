using KernelAbstractions: @index, @kernel
using Oceananigans.Operators

#++++ KERNEL COMPUTED FIELDS
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




#++++ CREATE SNAPSHOT OUTPUTS
function get_outputs_tuple(; LES=false, model=model)
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
    
    
    # Start calculation of snapshot variables
    #++++
    if LES
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
    
    wb_res = @at (Center, Center, Center) w*b
    tke = KineticEnergy(model, u, v, w)
    
    if LES
        ε = IsotropicViscousDissipation(model, νₑ, u, v, w)
    else
        ε = AnisotropicViscousDissipation(model, νx, νy, νz, u, v, w)
    end

    PV_ver = KernelComputedField(Face, Face, Face, ertel_potential_vorticity_vertical_fff!, model;
                                 computed_dependencies=(u_tot, v, b_tot), 
                                 parameters=f_0)
    
    PV_hor = KernelComputedField(Face, Face, Face, ertel_potential_vorticity_horizontal_fff!, model;
                                 computed_dependencies=(u_tot, v, w, b_tot), 
                                 parameters=f_0)
    
    dvpdy_ρ = PressureRedistribution_y(model, v, p, ρ₀)
    dwpdz_ρ = PressureRedistribution_z(model, w, p, ρ₀)

    
    SP_y = KernelComputedField(Center, Center, Center, shear_production_y_ccc!, model;
                               computed_dependencies=(u, v, w, U))
    
    SP_z = KernelComputedField(Center, Center, Center, shear_production_z_ccc!, model;
                               computed_dependencies=(u, v, w, U))
    #-----
    
    
    # Assemble the outputs tuple
    #++++
    outputs = (u=u,
               v=v,
               w=w,
               b=b,
               p=ComputedField(p),
               wb_res=ComputedField(wb_res),
               dwpdz_ρ=dwpdz_ρ,
               dvpdy_ρ=dvpdy_ρ,
               dbdz=ComputedField(dbdz),
               ω_x=ComputedField(ω_x),
               tke=tke,
               ε=ε,
               PV_ver=PV_ver,
               PV_hor=PV_hor,
               SP_y=SP_y,
               SP_z=SP_z,
               )
    
    if LES
        outputs = merge(outputs, (ν_e=νₑ,))
    end
    if as_background
        outputs = merge(outputs, (u_tot=ComputedField(u_tot),
                                  b_tot=ComputedField(b_tot),))
    end
    #----

    return outputs
end
#----


#++++ Construct outputs into simulation
function construct_outputs(; LES=false, model=model)
    
    # Output (high def) SNAPSHOTS
    #++++
    outputs_snap = get_outputs_tuple(LES=LES, model=model)
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
    delete(nt::NamedTuple{names}, keys) where names = NamedTuple{filter(x -> x ∉ keys, names)}(nt)
    
    outputs_vid = delete(outputs_snap, (:SP_y, :SP_z, :dwpdz_ρ, :dvpdy_ρ, :p))
    
    simulation.output_writers[:vid_writer] =
        NetCDFOutputWriter(model, outputs_vid,
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
    xz_average(F) = AveragedField(F, dims=(1,3))
    
    slicer = FieldSlicer(j=(grid.Ny÷frac):Int(grid.Ny*(1-1/frac)), with_halos=false)
    hor_window_average(F) = WindowedSpatialAverage(F; dims=(1, 2), field_slicer=slicer)
    
    outputs_avg = map(hor_window_average, outputs_snap)
    
    simulation.output_writers[:avg_writer] =
        NetCDFOutputWriter(model, outputs_avg,
                           filepath = @sprintf("avg.%s.nc", simname),
                           schedule = AveragedTimeInterval(10minutes; window=9.9minutes, stride=1),
                           mode = "c",
                           global_attributes = global_attributes,
                           array_type = Array{Float64},
                          )
    #-----
end
#-----
