import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIsurfjet1",
          "PNN_CIsurfjet2",
          "PNN_CIsurfjet3",
          "PNN_CIsurfjet4",
          "PNN_CIsurfjet5",
          "PNN_SIsurfjet1",
          "PNN_SIsurfjet2",
          "PNN_SIsurfjet3",
          "PNN_SIsurfjet4",
          "PNN_SIsurfjet5",
          "PNN_SIsurfjet6",
          "PNN_CIintjet01",
          ]
#----


#++++ Options
plot = True
#----

effslist = []
for sname in snames:
    #++++ Open datasets avg and vid
    print(f"Opening {sname}")
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #----

    #++++ Define integration, grid and metrics using pynanigans + xgcm
    indef_integrate = pn.utils.regular_indef_integrate
    #-----


    # Get important variables
    #++++
    if False:
        variables_avg = ["wb_res", "ε", "dvpdy_ρ", "dwpdz_ρ", "b_tot", "tke"]
        variables_vid = ["wb_res", "ε", "b_tot", "tke"]

        avg = avg[variables_avg]
        vid = vid[variables_vid]
    #---

    # BPE calculation
    #+++++++
    avg["BPE"] = -(avg.b_sorted * avg.zC)
    avg["dBPEdt"] = avg.BPE.load().differentiate("time")/avg.T_inertial
    #-------

 
    #++++ Get important times and perform average
    times_vid = [0, 0.5, 2, 3, 5, 10]
    times_avg = slice(None, None, None)
    
    avg = avg.sel(time=times_avg)
    vid = vid.sel(time=times_vid, method="nearest")
    avg_0d = avg.mean(("zC"))
    #-----
    
    
    # Calculations of efficiency
    #+++++
    ε_p = avg_0d.dBPEdt
    ε_k = avg_0d.ε

    denom = avg_0d.dBPEdt + ε_k
    denom2 = avg_0d.dBPEdt + ε_k + avg_0d.sponge_dissip
    
    ε_p_int = indef_integrate(ε_p, dim="time")
    ε_k_int = indef_integrate(ε_k, dim="time")
    denom_int = indef_integrate(denom, dim="time")
    denom_int2 = indef_integrate(denom2, dim="time")
    
    γ = ε_p / denom
    γ2 = ε_p / denom2

    Γ = ε_p_int / denom_int
    Γ2 = ε_p_int / denom_int2

    γ_coeff = ε_p / ε_k
    Γ_coeff = ε_p_int / ε_k_int
    #-----
    
    #++++
    ds_eff = xr.Dataset(dict(γ=γ, Γ=Γ, 
                             γ2=γ2, Γ2=Γ2,
                             γ_coeff=γ_coeff, Γ_coeff=Γ_coeff,
                             BPE=avg_0d.BPE, 
                             dBPEdt=avg_0d.dBPEdt,
                             ε=avg_0d.ε,
                             sponge_dissip=avg_0d.sponge_dissip,
                             ε_p_int=ε_p_int,
                             denom_int=denom_int, denom_int2=denom_int2,
                             ))
    #----

    #++++ Include last values to make my life easier
    ds_eff["γ_last"] = ds_eff.γ.isel(time=-1)
    ds_eff["Γ_last"] = ds_eff.Γ.isel(time=-1)
    ds_eff["γ2_last"] = ds_eff.γ2.isel(time=-1)
    ds_eff["Γ2_last"] = ds_eff.Γ2.isel(time=-1)
    ds_eff["γ_coeff_last"] = ds_eff.γ_coeff.isel(time=-1)
    #----

    #++++ Append reference quantities to list
    ds_eff.load()
    ds_eff = ds_eff.expand_dims("simulation").assign_coords(simulation=[sname])
    effslist.append(ds_eff)
    #----

    #++++ Save to disk
    ds_eff.to_netcdf(f"data/efficiencies_{sname}.nc")
    #----
    
#++++ Create and save dataset
alleffs = xr.concat(effslist, dim="simulation", combine_attrs="drop")
alleffs.to_netcdf("data/alleffs.nc")
#----

