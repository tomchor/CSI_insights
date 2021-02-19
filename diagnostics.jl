using KernelAbstractions: @index, @kernel
using Oceananigans.Operators

@inline ψ²(i, j, k, grid, ψ) = @inbounds ψ[i, j, k]^2
@inline ψ′²(i, j, k, grid, ψ, Ψ) = @inbounds (ψ[i, j, k] - Ψ[i, j, k])^2

#@kernel function kinetic_energy_ccc!(tke, grid, u, v, w)
#    i, j, k = @index(Global, NTuple)
#
#    @inbounds tke[i, j, k] = (
#                              ℑxᶜᵃᵃ(i, j, k, grid, ψ², u) +
#                              ℑyᵃᶜᵃ(i, j, k, grid, ψ², v) +
#                              ℑzᵃᵃᶜ(i, j, k, grid, ψ², w)
#                             ) / 2
#end
#
#
#
#@kernel function isotropic_viscous_dissipation_ccc!(ϵ, grid, ν, u, v, w)
#    i, j, k = @index(Global, NTuple)
#
#    Σˣˣ = ∂xᶜᵃᵃ(i, j, k, grid, u)
#    Σʸʸ = ∂yᵃᶜᵃ(i, j, k, grid, v)
#    Σᶻᶻ = ∂zᵃᵃᶜ(i, j, k, grid, w)
#
#    Σˣʸ = (ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᵃᶠᵃ, u) + ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᵃᵃ, v)) / 2
#    Σˣᶻ = (ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᵃᵃᶠ, u) + ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᵃᵃ, w)) / 2
#    Σʸᶻ = (ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᵃᵃᶠ, v) + ℑyzᵃᶜᶜ(i, j, k, grid, ∂yᵃᶠᵃ, w)) / 2
#
#    @inbounds ϵ[i, j, k] = ν[i, j, k] * 2 * (Σˣˣ^2 + Σʸʸ^2 + Σᶻᶻ^2 + 2 * (Σˣʸ^2 + Σˣᶻ^2 + Σʸᶻ^2))
#end
#
#
#
#@kernel function anisotropic_viscous_dissipation_ccc!(ϵ, grid, νx, νy, νz, u, v, w)
#    i, j, k = @index(Global, NTuple)
#
#    ddx² = ∂xᶜᵃᵃ(i, j, k, grid, u)^2 + ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᵃᵃ, v)^2 + ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᵃᵃ, w)^2
#
#    ddy² = ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᵃᶠᵃ, u)^2 + ∂yᵃᶜᵃ(i, j, k, grid, v)^2 + ℑyzᵃᶜᶜ(i, j, k, grid, ∂yᵃᶠᵃ, w)^2
#
#    ddz² = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᵃᵃᶠ, u)^2 + ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᵃᵃᶠ, v)^2 + ∂zᵃᵃᶜ(i, j, k, grid, w)^2
#
#    @inbounds ϵ[i, j, k] = νx[i,j,k]*ddx² + νy[i,j,k]*ddy² + νz[i,j,k]*ddz²
#end
#
#
#
#@kernel function richardson_number_ccf!(Ri, grid, u, v, b, dUdz_bg, dVdz_bg, N2)
#    i, j, k = @index(Global, NTuple)
#
#    dBdz = ∂zᵃᵃᶠ(i, j, k, grid, b) + N2 # dbdz(c, c, f)
#    dUdz_tot = ℑxᶜᵃᵃ(i, j, k, grid, ∂zᵃᵃᶠ, u) + dUdz_bg # dudz(f, c, f) → dudz(c, c, f)
#    dVdz_tot = ℑyᵃᶜᵃ(i, j, k, grid, ∂zᵃᵃᶠ, v) + dVdz_bg # dvdz(c, f, f) → dvdz(c, c, f)
#
#    @inbounds Ri[i, j, k] = dBdz / (dUdz_tot^2 + dVdz_tot^2)
#end
#
#
#
#
#
#@kernel function rossby_number_ffc!(Ro, grid, u, v, dUdy_bg, dVdx_bg, f₀)
#    i, j, k = @index(Global, NTuple)
#
#    dUdy_tot = ∂yᵃᶠᵃ(i, j, k, grid, u) + dUdy_bg
#    dVdx_tot = ∂xᶠᵃᵃ(i, j, k, grid, v) + dVdx_bg
#
#    @inbounds Ro[i, j, k] = (dVdx_tot - dUdy_tot) / f₀
#end
#
#
#
#@kernel function compute_pv_from_Ro_Ri!(PV, grid, Ri, Ro, N², f₀)
#    i, j, k = @index(Global, NTuple)
#
#    Ro_int = ℑxyzᶜᶜᶠ(i, j, k, grid, Ro)
#
#    @inbounds PV[i, j, k] = N²[i, j, k]*f₀ * (1 + Ro[i, j, k] - 1/Ri[i, j, k])
#    #@inbounds PV[i, j, k] = N²[i, j, k]*f₀ * (1 + Ro_int - 1/Ri[i, j, k])
#end
#
#
#@kernel function potential_vorticity_in_thermal_wind_fff!(PV, grid, u, v, b, f₀)
#    i, j, k = @index(Global, NTuple)
#
#    dVdx =  ℑzᵃᵃᶠ(i, j, k, grid, ∂xᶠᵃᵃ, v) # F, F, C → F, F, F
#    dUdy =  ℑzᵃᵃᶠ(i, j, k, grid, ∂yᵃᶠᵃ, u) # F, F, C → F, F, F
#    dbdz = ℑxyᶠᶠᵃ(i, j, k, grid, ∂zᵃᵃᶠ, b) # C, C, F → F, F, F
#
#    pv_barot = (f₀ + dVdx - dUdy) * dbdz
#
#    dUdz = ℑyᵃᶠᵃ(i, j, k, grid, ∂zᵃᵃᶠ, u) # F, C, F → F, F, F
#    dVdz = ℑxᶠᵃᵃ(i, j, k, grid, ∂zᵃᵃᶠ, v) # C, F, F → F, F, F
#
#    pv_baroc = -f₀*(dUdz^2 + dVdz^2)
#
#    @inbounds PV[i, j, k] = pv_barot[i, j, k] + pv_baroc[i, j, k]
#end
#
#
#
#@kernel function ertel_potential_vorticity_fff!(PV, grid, u, v, w, b, f₀)
#    i, j, k = @index(Global, NTuple)
#
#    dWdy =  ℑxᶠᵃᵃ(i, j, k, grid, ∂yᵃᶠᵃ, w) # C, C, F  → C, F, F  → F, F, F
#    dVdz =  ℑxᶠᵃᵃ(i, j, k, grid, ∂zᵃᵃᶠ, v) # C, F, C  → C, F, F  → F, F, F
#    dbdx = ℑyzᵃᶠᶠ(i, j, k, grid, ∂xᶠᵃᵃ, b) # C, C, C  → F, C, C  → F, F, F
#    pv_x = (dWdy - dVdz) * dbdx # F, F, F
#
#    dUdz =  ℑyᵃᶠᵃ(i, j, k, grid, ∂zᵃᵃᶠ, u) # F, C, C  → F, C, F → F, F, F
#    dWdx =  ℑyᵃᶠᵃ(i, j, k, grid, ∂xᶠᵃᵃ, w) # C, C, F  → F, C, F → F, F, F
#    dbdy = ℑxzᶠᵃᶠ(i, j, k, grid, ∂yᵃᶠᵃ, b) # C, C, C  → C, F, C → F, F, F
#    pv_y = (dUdz - dWdx) * dbdy # F, F, F
#
#    dVdx =  ℑzᵃᵃᶠ(i, j, k, grid, ∂xᶠᵃᵃ, v) # C, F, C  → F, F, C → F, F, F
#    dUdy =  ℑzᵃᵃᶠ(i, j, k, grid, ∂yᵃᶠᵃ, u) # F, C, C  → F, F, C → F, F, F
#    dbdz = ℑxyᶠᶠᵃ(i, j, k, grid, ∂zᵃᵃᶠ, b) # C, C, C  → C, C, F → F, F, F
#    pv_z = (f₀ + dVdx - dUdy) * dbdz
#
#    @inbounds PV[i, j, k] = pv_x + pv_y + pv_z
#end
#
#
#
#
#
#@kernel function wp_ccc!(wp, grid, w, p)
#    i, j, k = @index(Global, NTuple)
#
#    @inbounds wp[i, j, k] = ℑzᵃᵃᶜ(i, j, k, grid, w) * p[i, j, k]
#end
#
#
#@kernel function vertical_pressure_distribution_ccc!(dwpdz, grid, w, p, ρ₀)
#    i, j, k = @index(Global, NTuple)
#
#    wp = ℑzᵃᵃᶠ(i, j, k, grid, p) * w[i, j, k] # C, C, F
#
#    @inbounds dwpdz[i, j, k] = (1/ρ₀) * ∂zᵃᵃᶜ(i, j, k, grid, wp) # C, C, F  → C, C, C
#end
#

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


function custom_average(model)
    U = AveragedField(model.velocities.u, dims=1)
    compute!(U)    
    U_slice_averaged = mean(U.data[1, 20:50, 1:grid.Nz], dims=2)
    return U_slice_averaged
end
