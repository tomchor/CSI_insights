import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from aux01_physfuncs import adjust_variables, get_IC_front
from dask.diagnostics import ProgressBar
π = np.pi

#++++ Define directory and simulation name
path = f"simulations/data/"
snames = ["PNN_CIfront1",
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
          #"PNN_CIintjet01",
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


IClist = []
for sname in snames:
    #++++ Open datasets avg and vid
    print(sname)
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    vid = adjust_variables(vid)
    #-----

    #+++++ Get parameters and create simulation dimension
    vid0 = vid.sel(time=0)
    IC = get_IC_front(y=vid.yF, z=vid.zF, 
                      y_0=vid.y_0, z_0=vid.z_0, 
                      σ_y=vid.σ_y, σ_z=vid.σ_z, 
                      u_0=vid.u_0, N2_inf=vid.N2_inf, f_0=
                      vid.f_0, Lz=vid.Lz)

    if False:
        from matplotlib import pyplot as plt
        IC.b.pncontourf(y='z',)
        plt.figure()
        vid0.b.pncontourf(y='z',)

    if False:
        from matplotlib import pyplot as plt
        opts = dict(vmin=-2.5, vmax=2.5, cmap="RdBu_r")
        IC.PV.pnimshow(y='z', **opts)
        PV_c = 0*float(IC.PV.min())
        IC.PV.pncontour(y='z', levels=[PV_c], c='k')
        #plt.figure()
        #(vid0.PV/(vid.f_0*vid.N2_inf)).pnimshow(y='z', **opts)
    #----

    #++++ Interpolate to get masks
    IC = IC.rename(q_hat="PV_aff")
    IC["PV_acc"] = grid_vid.interp(IC.PV_aff, ['y', 'z'])
    IC["PV_acf"] = grid_vid.interp(IC.PV_aff, 'y')
    IC["PV_afc"] = grid_vid.interp(IC.PV_aff, 'z')

    IC["mask_acc"] = IC.PV_acc < 0
    IC["mask_acf"] = IC.PV_acf < 0
    IC["mask_afc"] = IC.PV_afc < 0
    #----

    #++++ Expand simulation coord and append to list
    IC = IC.expand_dims("simulation").assign_coords(simulation=[sname])
    IClist.append(IC)
    #----
    
#++++ Create and save dataset
allICs = xr.concat(IClist, dim="simulation")
with ProgressBar():
    allICs.to_netcdf("data_post/allICs.nc")
#----

