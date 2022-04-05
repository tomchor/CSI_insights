import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
π = np.pi

#++++ Define directory and simulation name
path = f"simulations/data/"
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
          #"PNN_CIintjet01",
          "PNN_CIsurfjet1_f2",
          "PNN_CIsurfjet1_f4",
          "PNN_CIsurfjet1_f8",
          "PNN_SIsurfjet4_f2",
          "PNN_SIsurfjet4_f4",
          "PNN_SIsurfjet4_f8",
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
    times_avg = slice(None, None, None)
    
    avg = avg.sel(time=times_avg)
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

    #++++ Calculate the Ozmidov scale
    t_max = float(avg_0d.where(avg_0d.ε == avg_0d.ε.max(), drop=True).time)
    def_times = [0, 0.5, 2, 3, 4, 5, 8, 10]
    vid_reltimes = vid.sel(time=slice(t_max, None, 5))
    #times_vid = np.sort(np.append(t_max,  [el for el in def_times if el > t_max ]))
    #vid_reltimes = vid.sel(time=times_vid, method="nearest")


    #++++++ Max(<ε>) calculation
    if False:
        print("Calculating ε_max by percentile")
        ε_max = np.percentile(vid_reltimes.ε.sel(time=t_max, method="nearest"), 99.)

    elif False:
        print("Calculating ε_max by coarsening")
        y_window = int(100//vid_reltimes.ΔyF.max()) # 100 m resolution in horizontal
        z_window = int(10//vid_reltimes.ΔzF.max()) # 10 m resolution in vertical

        ε_max = vid_reltimes.ε.coarsen(yC=y_window, zC=z_window, boundary="trim").mean().max()
    else:
        print("Don't calculate ε_max")
        ε_max = 1e-10
    #------

    #++++ Filter out places with low dissipation
    mask_ccc = vid_reltimes.ε >= ε_max
    mask_ccf = grid_vid.interp(vid_reltimes.ε, 'z', boundary="fill") >= ε_max

    ε_filt = vid_reltimes.ε.where(mask_ccc, drop=True)
    dbdz_filt = vid_reltimes.dbdz.where(mask_ccf, drop=True)
    dbdz_filt2 = grid_vid.interp(vid_reltimes.dbdz, 'z').where(mask_ccc, drop=True)

    N_filt = np.sqrt(dbdz_filt)
    N_filt2 = np.sqrt(dbdz_filt2)
    #----

    #+++++ Ozmidov scale calculation
    N_inf = np.sqrt(vid_reltimes.N2_inf)
    N_avg = np.sqrt(dbdz_filt.pnmean())
    ε_avg = ε_filt.pnmean()

    Lo_inf = 2*π*np.sqrt(ε_avg / N_inf**3) # Using √(<ε>/N_inf³)
    Lo_avg = 2*π*np.sqrt(ε_avg / N_avg**3) # Using √(<ε>/<N>³)
    Lo = 2*π * (np.sqrt(ε_filt / N_filt2**3)).pnmean() # Using < √(ε/N³) >
    #-----
    #----

    
    #++++ Create dataset
    ds_eff = xr.Dataset(dict(γ=γ, Γ=Γ, 
                             γ2=γ2, Γ2=Γ2,
                             γ_coeff=γ_coeff, Γ_coeff=Γ_coeff,
                             BPE=avg_0d.BPE, 
                             dBPEdt=avg_0d.dBPEdt,
                             ε=avg_0d.ε,
                             sponge_dissip=avg_0d.sponge_dissip,
                             ε_p_int=ε_p_int,
                             denom_int=denom_int, denom_int2=denom_int2,
                             Lo_inf=Lo_inf,
                             Lo_avg=Lo_avg,
                             Lo=Lo,
                             Δz = float(vid_reltimes.ΔzC.median()),
                             Δy = float(vid_reltimes.ΔyC.median()),
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
    ds_eff.to_netcdf(f"data_post/efficiencies_{sname}.nc")
    #----
    
#++++ Create and save dataset
alleffs = xr.concat(effslist, dim="simulation", combine_attrs="drop")
alleffs.to_netcdf("data_post/alleffs.nc")
#----

