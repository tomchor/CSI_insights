import numpy as np
import xarray as xr
from aux00_utils import open_simulation
π = np.pi


#++++ Define directory and simulation name
path = f"simulations/data/"
snames = [#"FNN_CIfront1",
          #"FNN_SIfront4",
          "PNN_CIfront1",
          "PNN_CIfront2",
          "PNN_CIfront3",
          "PNN_CIfront4",
          "PNN_CIfront5",
          "PNN_SIfront1",
          "PNN_SIfront2",
          "PNN_SIfront3",
          "PNN_SIfront4",
          "PNN_SIfront5",
          "PNN_SIfront6",
          "PNN_CIfront1_f2",
          "PNN_CIfront1_f4",
          "PNN_CIfront1_f8",
          "PNN_SIfront4_f2",
          "PNN_SIfront4_f4",
          "PNN_SIfront4_f8",
          "PNN_CIfront1_AMD",
          "PNN_CIfront1_AMD_f2",
          "PNN_CIfront1_AMD_f4",
          "PNN_CIfront1_AMD_f8",
          "PNN_SIfront4_AMD",
          "PNN_SIfront4_AMD_f2",
          "PNN_SIfront4_AMD_f4",
          "PNN_SIfront4_AMD_f8",
          "PNN_CIfront1_AMD2",
          "PNN_CIfront1_AMD2_f2",
          "PNN_CIfront1_AMD2_f4",
          "PNN_CIfront1_AMD2_f8",
          "PNN_SIfront4_AMD2",
          "PNN_SIfront4_AMD2_f2",
          "PNN_SIfront4_AMD2_f4",
          "PNN_SIfront4_AMD2_f8",
          "PNN_CIfront1_AMD3",
          "PNN_CIfront1_AMD3_f2",
          "PNN_CIfront1_AMD3_f4",
          "PNN_CIfront1_AMD3_f8",
          "PNN_SIfront4_AMD3",
          "PNN_SIfront4_AMD3_f2",
          "PNN_SIfront4_AMD3_f4",
          "PNN_SIfront4_AMD3_f8",
          ]
#----

#++++ Options
plot = False
#----


dslist = []
for sname in snames:
    #++++ Open datasets avg and vid
    print(f"Opening {sname}")
    grid_out, out = open_simulation(path+f"out.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
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



    #+++++ Get proper times and window
    #+++++ Growth cycle method
    q_hat_r = 1 + vid.Ro_r - 1/vid.Ri_r
    T_growth = 1/np.sqrt(-vid.f_0**2 * q_hat_r)
    T_gr_hat = T_growth/vid.T_inertial
    #-----

    #+++++ Dissipation max method
    ε_avg = avg.ε.pnmean('z')
    t_εmax = ε_avg.time[ε_avg.argmax('time').values]
    #-----

    #+++++ Choose actual time for lin growth
    t_lin = 15*T_gr_hat
    #t_lin = 3*t_εmax/4

    vid = vid.sel(time=slice(0, t_lin))
    avg = avg.sel(time=slice(0, t_lin))
    #-----
    #-----
    
    #+++++ Get HSP and LSP
    if "shearprod_y" in avg.variables:
        SPy = avg.shearprod_y
        SPz = avg.shearprod_z

    else:

        #+++++ Calculate prelim variables for LSP and HSP
        print("Calculating shear production rates manually")
        U_bg = out.U.sel(xF=0)
        dUdz = grid_out.derivative(U_bg, 'z', boundary="extend")
        dUdy = grid_out.derivative(U_bg, 'y', boundary="extend")

        u_acc = vid.u_tot - U_bg
        u_afc = grid_vid.interp(u_acc, 'y', boundary="extrapolate")
        u_acf = grid_vid.interp(u_acc, 'z', boundary="extrapolate")
        v = vid.v
        w = vid.w
        #-----

        #+++++ Calculate LSP and HSP
        SPy = -(u_afc * v * dUdy).mean(('yF'))
        SPz = -(u_acf * w * dUdz).mean(('yC'))
        #-----
    #-----
    
    #+++++ Save
    dsout = xr.Dataset(dict(SPy=SPy, SPz=SPz))
    dsout["SPy_avg"] = dsout.SPy.pnmean(('z'))
    dsout["SPz_avg"] = dsout.SPz.pnmean(('z'))
    dsout["SPy_lin"] = dsout.SPy_avg.isel(time=-1)
    dsout["SPz_lin"] = dsout.SPz_avg.isel(time=-1)

    dsout.to_netcdf(f"data_post/shearprod_{sname}.nc")

    dsout = dsout.expand_dims("simulation").assign_coords(simulation=[sname])
    dslist.append(dsout)
    #-----
   
    
    #+++++ Plot
    if plot:
        print("Plotting...")
        from matplotlib import pyplot as plt
        size=5
        ncols=1; nrows=2
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*size, 0.7*nrows*size),
                                 squeeze=False, sharey=False)
        axesf = axes.flatten()

        def plot_SP(ds, axesf):
            (ds.SPy_avg/abs(ds.SPy_avg+ds.SPz_avg)).plot.line(x='time', 
                                              ax=axesf[0], 
                                              label='horizontal SP / abs(total SP)')
            ds.ω_x.plot(ax=axesf[1], rasterized=True, label="vorticity")
            return

        dsout["ω_x"] = vid.ω_x.isel(time=-1)
        plot_SP(dsout, axesf)
    
    
        #++++ Polish and Save plot
        for ax in axesf:
            ax.grid(True)
            ax.legend()
        fig.tight_layout()
        fig.savefig(f"figures_check/SP_{sname}.png")
        #----
    #-----


#++++ Create and save dataset
print("Saving allshearprod...")
allshearprod = xr.concat(dslist, dim="simulation")
allshearprod.to_netcdf("data_post/allshearprod.nc")
#----
