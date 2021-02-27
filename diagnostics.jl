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


#++++ Windowed Average
import Oceananigans.OutputWriters: slice_parent
using Statistics: mean

struct WindowedSpatialAverage{F, S, D}
          field :: F
   field_slicer :: S
           dims :: D
end


WindowedSpatialAverage(field; dims, field_slicer=FieldSlicer()) = WindowedSpatialAverage(field, field_slicer, dims)

function (wsa::WindowedSpatialAverage)(model)
    compute!(wsa.field)
    window = slice_parent(wsa.field_slicer, wsa.field)
    return dropdims(mean(window, dims=wsa.dims), dims=wsa.dims)
end


#++++ Make dims detection 
using NCDatasets: defVar
using Oceananigans.Fields: reduced_location
import Oceananigans.OutputWriters: xdim, ydim, zdim, define_output_variable!
function define_output_variable!(dataset,
                                 wtsa::Union{WindowedSpatialAverage, WindowedTimeAverage{<:WindowedSpatialAverage}},
                                 name, array_type, compression, attributes, dimensions)
    wsa = wtsa isa WindowedTimeAverage ? wtsa.operand : wtsa
    LX, LY, LZ = reduced_location(location(wsa.field), dims=wsa.dims)

    output_dims = tuple(xdim(LX)..., ydim(LY)..., zdim(LZ)...)
    defVar(dataset, name, eltype(array_type), (output_dims..., "time"),
           compression=compression, attrib=attributes)
    return nothing
end
#----

#----



#++++ PV components

@kernel function ertel_potential_vorticity_barotropic_fff!(PV, grid, u, v, b, f₀)
    i, j, k = @index(Global, NTuple)

    dVdx =  ℑzᵃᵃᶠ(i, j, k, grid, ∂xᶠᵃᵃ, v) # C, F, C  → F, F, C → F, F, F
    dUdy =  ℑzᵃᵃᶠ(i, j, k, grid, ∂yᵃᶠᵃ, u) # F, C, C  → F, F, C → F, F, F
    dbdz = ℑxyᶠᶠᵃ(i, j, k, grid, ∂zᵃᵃᶠ, b) # C, C, C  → C, C, F → F, F, F
    pv_z = (f₀ + dVdx - dUdy) * dbdz

    @inbounds PV[i, j, k] = pv_z
end



@kernel function ertel_potential_vorticity_baroclinic_fff!(PV, grid, u, v, w, b, f₀)
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


#----
