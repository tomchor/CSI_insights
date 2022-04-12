import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from dask.diagnostics import ProgressBar
from matplotlib import pyplot as plt
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIfront1",
#          "PNN_CIfront2",
#          "PNN_CIfront3",
#          "PNN_SIfront2",
#          "PNN_SIfront3",
#          "PNN_SIfront4",
#          "PNN_SIfront5",
#          "PNN_SIfront6",
#          "PNN_CIintjet01",
          ]
#----

#++++ Options
variables = ["u", "v", "w"]
#----


for i, sname in enumerate(snames):
    print(f"Opening {sname}")

    #++++ Open datasets avg and vid
    print("Calculating Reb")
    grid_out, out = open_simulation(path+f"out.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #----

    #++++ Get subdomain
    out = out.sel(time=5, method="nearest")

    if "front" in sname:
        y_r = 1/4 * np.sqrt(2) * out.σ_y + out.y_0
    else:
        y_r = 1/2 * np.sqrt(2) * out.σ_y + out.y_0
    side = 2*out.σ_y
    yF_r = out.sel(yF=y_r+side, method="nearest").yF
    yF_l = out.sel(yF=y_r-side, method="nearest").yF
    yslice = slice(yF_l, yF_r) # Make sure len(yF) = len(yC)+1

    depth=np.inf
    zF_r = out.sel(zF=out.z_0+depth, method="nearest").zF
    zF_l = out.sel(zF=out.z_0-depth, method="nearest").zF
    zslice = slice(zF_l, zF_r) # Make sure len(zF) = len(zC)+1

    out = out.sel(yC=yslice, yF=yslice, zC=zslice, zF=zslice).reset_coords()
    subgrid_out = pn.get_grid(out, topology="PNN")

    ix = len(out.xC)//2; x0 = dict(x=out.xC[ix], method="nearest")
    kz = len(out.zC)//2; z0 = dict(z=out.zC[kz], method="nearest")
    #----


    for var in variables:
        print(var)
        print()

        #++++ Calculations
        out[f"{var}_mean"] = out[var].biject('x').mean("x")
        out[f"{var}_fluc"] = out[var] - out[f"{var}_mean"]
        #----
    
        #+++++ Start figure
        nrows=2; ncols=3
        size = 4
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.3*ncols*size, nrows*size),
                                         constrained_layout=True, squeeze=False)
        axesf = axes.flatten()
        #-----

        #+++++ Plot
        with ProgressBar():
            out[var].biject().sel(**x0).plot.imshow(ax=axesf[0], x='y', rasterized=True)
            out[f"{var}_fluc"].biject().sel(**x0).plot.imshow(ax=axesf[1], x='y', rasterized=True)
            out[f"{var}_mean"].biject().plot.imshow(ax=axesf[2], x='y', rasterized=True)

            out[var].biject().sel(**z0).plot.imshow(ax=axesf[3], x='y', rasterized=True)
            out[f"{var}_fluc"].biject().sel(**z0).plot.imshow(ax=axesf[4], x='y', rasterized=True)
            out.ω_x.isel(xC=ix).pnimshow(ax=axesf[5], x='y', rasterized=True)
        #-----


        #++++ Save plot
        fig.savefig(f"figures_check/xavg_{var}_{sname}.pdf")
        #----
