using Printf
using Statistics: mean
import NCDatasets as NCD

using Oceananigans.AbstractOperations: @at, ∂x, ∂y, ∂z
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.Fields: ComputedField, KernelComputedField
using Oceananigans.Diagnostics: WindowedSpatialAverage
using Oceananigans.Grids: Center, Face

using Oceanostics: KineticEnergy, 
                   IsotropicViscousDissipationRate,
                   AnisotropicPseudoViscousDissipationRate,
                   YPressureRedistribution, ZPressureRedistribution,
                   YShearProduction, ZShearProduction
using Oceanostics.FlowDiagnostics: ErtelPotentialVorticityᶠᶠᶠ



#++++ Unpack model variables
ccc_scratch = Field(Center, Center, Center, model.architecture, model.grid)
fcc_scratch = Field(Face, Center, Center, model.architecture, model.grid)
ccf_scratch = Field(Center, Center, Face, model.architecture, model.grid)
cff_scratch = Field(Center, Face, Face, model.architecture, model.grid)
fff_scratch = Field(Face, Face, Face, model.architecture, model.grid)

u, v, w = model.velocities
b = model.tracers.b
p = sum(model.pressures)

U, V, W, = Oceananigans.Fields.BackgroundVelocityFields((u=u_g,), model.grid, model.clock)
B, = Oceananigans.Fields.BackgroundTracerFields((b=b_g,), (:b,), model.grid, model.clock)

if LES
    νₑ = νz = model.diffusivity_fields.νₑ
    if AMD
        κₑ = κz = model.diffusivity_fields.κₑ.b
    else
        κₑ = κz = ComputedField(model.diffusivity_fields.κₑ.b)
    end
else
    if model.closure isa IsotropicDiffusivity
        νx = νy = νz = model.closure.ν
        κx = κy = κz = model.closure.κ[:b]
    elseif model.closure isa AnisotropicDiffusivity
        νx = model.closure.νx; νy = model.closure.νy; νz = model.closure.νz
        κx = model.closure.κx[:b]; κy = model.closure.κy[:b]; κz = model.closure.κz[:b]
    end
end
#---

#++++ CREATE SNAPSHOT OUTPUTS
function get_outputs_tuple(model; LES=false)
    
    #++++ Create mask fields
    mask_u = Oceananigans.Fields.FunctionField{Face, Center, Center}(full_mask, model.grid)
    mask_v = Oceananigans.Fields.FunctionField{Center, Face, Center}(full_mask, model.grid)
    mask_w = Oceananigans.Fields.FunctionField{Center, Center, Face}(full_mask, model.grid)
    #----


    # Start calculation of snapshot variables
    #++++
    dbdz = @at (Center, Center, Face) ∂z(b)
    ω_x = ∂y(w) - ∂z(v)
    
    wb_res = @at (Center, Center, Center) w*b
    tke = ComputedField(KineticEnergy(model), data=ccc_scratch.data)
    
    if LES
        ε = ComputedField(IsotropicViscousDissipationRate(model, u, v, w, νₑ), data=ccc_scratch.data)
    else
        ε = ComputedField(AnisotropicPseudoViscousDissipationRate(model, u, v, w, νx, νy, νz), data=ccc_scratch.data)
    end

    u_dissip = ComputedField((@at (Center, Center, Center) u * rate * mask_u * (u - U)), data=ccc_scratch.data)
    v_dissip = ComputedField((@at (Center, Center, Center) v * rate * mask_v * v), data=ccc_scratch.data)
    w_dissip = ComputedField((@at (Center, Center, Center) w * rate * mask_w * w), data=ccc_scratch.data)
    sponge_dissip = @at (Center, Center, Center) (u_dissip + v_dissip + w_dissip)

    PV = ComputedField(ErtelPotentialVorticityᶠᶠᶠ(model), data=fff_scratch.data)

    dvpdy = ComputedField(YPressureRedistribution(model), data=ccc_scratch.data)
    dwpdz = ComputedField(ZPressureRedistribution(model), data=ccc_scratch.data)
    
    shearprod_y = ComputedField(YShearProduction(model, u-U, v, w, U, 0, 0), data=ccc_scratch.data)
    shearprod_z = ComputedField(ZShearProduction(model, u-U, v, w, U, 0, 0), data=ccc_scratch.data)
    #-----
    
    
    # Assemble the outputs tuple
    #++++
    outputs = (u=u,
               v=v,
               w=w,
               b=b,
               p=ComputedField(p, data=ccc_scratch.data),
               wb_res=ComputedField(wb_res, data=ccc_scratch.data),
               PV=PV,
               dwpdz=dwpdz,
               dvpdy=dvpdy,
               dbdz=ComputedField(dbdz, data=ccf_scratch.data),
               ω_x=ComputedField(ω_x, data=cff_scratch.data),
               tke=tke,
               ε=ε,
               shearprod_y=shearprod_y,
               shearprod_z=shearprod_z,
               sponge_dissip=ComputedField(sponge_dissip, data=ccc_scratch.data),
               )
    
    if LES
        outputs = merge(outputs, (ν_e=νₑ, κ_e=κₑ))
    end
    #----

    return outputs
end
#----

#++++ Write background variables to dataset
function write_to_ds(dsname, varname, data; coords=("xC", "yC", "zC"), dtype=Float64)
    ds = NCD.NCDataset(dsname, "a")
    newvar = NCD.defVar(ds, varname, dtype, coords)
    newvar[:,:,:] = Array(data)
    NCD.close(ds)
end

function save_UB(dsname)
    Ucf = ComputedField(U+0*u, data=fcc_scratch.data); compute!(Ucf)
    Bcf = ComputedField(B+0*b, data=ccc_scratch.data); compute!(Ucf)

    write_to_ds(dsname, "U", interior(Ucf), coords=("xF", "yC", "zC"))
    write_to_ds(dsname, "B", interior(Bcf), coords=("xC", "yC", "zC"))
end
#----


#++++ Construct outputs into simulation
function construct_outputs(model, simulation; 
                           LES=false, simname="TEST", frac=1/16,
                           )

    #++++ Check for checkpoints
    if any(startswith("chk.$(simname)_iteration"), readdir("data"))
        @info "Checkpoint for $simname found. Assuming this is a pick-up simulation! Setting mode to append."
        mode = "a"
    else
        @info "No checkpoint for $simname found. Setting mode to clobber."
        mode = "c"
    end
    #----

    #++++ Define sorting b function
    function flattenedsort(A, dim_order::Union{Tuple, AbstractVector})
        return reshape(sort(Array(permutedims(A, dim_order)[:])), (grid.Nx, grid.Ny, grid.Nz))
    end

    function sort_b(model; average=false)
        b = model.tracers.b
        sorted_B = flattenedsort(interior(b), [3,2,1])
        if !average
            return sorted_B
        else
            return dropdims(mean(sorted_B, dims=(1,2)), dims=(1,2))
        end
    end
    mean_sort_b = (mod)->sort_b(mod; average=true)
    #-----

    
    # Output (high def) SNAPSHOTS
    #++++
    @info "Setting up out writer"
    outputs_snap = get_outputs_tuple(model, LES=LES)
    simulation.output_writers[:out_writer] = out_writer =
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
    
    outputs_vid = delete(outputs_snap, (:dwpdz, :dvpdy, :p))
    if ndims==2 outputs_vid = merge(outputs_vid, (b_sorted=sort_b,)) end
    dims = Dict("b_sorted" => ("xC", "yC", "zC"),)
    
    simulation.output_writers[:vid_writer] =
        NetCDFOutputWriter(model, outputs_vid,
                           filepath = @sprintf("data/vid.%s.nc", simname),
                           schedule = TimeInterval(60minutes),
                           mode = mode,
                           global_attributes = global_attributes,
                           dimensions = dims,
                           array_type = Array{Float32},
                           field_slicer = FieldSlicer(i=1, with_halos=false),
                           )
    #-----
    
    
    # AVG outputs
    #++++
    @info "Setting up avg writer"
    x_average(F) = AveragedField(F, dims=(1,))
    xy_average(F) = AveragedField(F, dims=(1,2))
    xz_average(F) = AveragedField(F, dims=(1,3))
    
    #slicer = FieldSlicer(j=Int(grid.Ny*frac):Int(grid.Ny*(1-1*frac)), with_halos=false)
    hor_window_average(F) = WindowedSpatialAverage(F; dims=(1, 2))
    
    outputs_avg = map(hor_window_average, outputs_snap)
    #outputs_avg = map(xy_average, outputs_snap)
    outputs_avg = merge(outputs_avg, (b_sorted=mean_sort_b,))
    dims = Dict("b_sorted" => ("zC",),)
    
    simulation.output_writers[:avg_writer] =
        NetCDFOutputWriter(model, outputs_avg,
                           filepath = @sprintf("data/avg.%s.nc", simname),
                           schedule = TimeInterval(2minutes),
                           mode = mode,
                           global_attributes = global_attributes,
                           dimensions = dims,
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


    #++++ Write Background fields
    if mode=="c"
        save_UB("data/out.$simname.nc")
        save_UB("data/avg.$simname.nc")
    end
    #----

    return checkpointer
end
#-----


