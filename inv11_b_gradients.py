import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, prettify_names
from matplotlib import pyplot as plt
from dask.diagnostics import ProgressBar
from cmocean import cm


#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIsurfjet1",
          "PNN_SIsurfjet4",
          ]
#----


for i, sname in enumerate(snames):
    #++++ Open datasets avg and vid
    print(sname)
    grid_out, out = open_simulation(path+f"out.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #-----


    #++++ Get subdomain
    print("Getting a σ_y-wide slice of the domain")
    if "surfjet" in sname:
        y_r = 1/4 * np.sqrt(2) * out.σ_y + out.y_0
    else:
        y_r = 1/2 * np.sqrt(2) * out.σ_y + out.y_0
    side = 1/2*out.σ_y
    yF_r = out.sel(yF=y_r+side, method="nearest").yF
    yF_l = out.sel(yF=y_r-side, method="nearest").yF
    yslice = slice(yF_l, yF_r) # Make sure len(yF) = len(yC)+1

    depth=np.inf
    zF_r = out.sel(zF=out.z_0+depth, method="nearest").zF
    zF_l = out.sel(zF=out.z_0-depth, method="nearest").zF
    zslice = slice(zF_l, zF_r) # Make sure len(zF) = len(zC)+1

    tslice = slice(None, 3)

    out = out.sel(yC=yslice, yF=yslice, zC=zslice, zF=zslice, time=tslice)
    subgrid_out = pn.get_grid(out, topology="PNN")
    #-----

    #++++ Calculations
    def get_gradients(ds, grid=None):
        if grid is None: 
            grid = pn.get_grid(ds, topology="PNN")

        ds["dbdx"] = grid.interp(grid.derivative(ds.b, 'x'), 'x')
        ds["dbdy"] = grid.interp(grid.derivative(ds.b, 'y', boundary="extrapolate"), 'y')
        ds["dbdz"] = grid.interp(ds.dbdz, 'z')

        ds["dbds"] = np.sqrt(ds.dbdx**2 + ds.dbdy**2 + ds.dbdz**2)
        ds["dbds_mean"] = ds.dbds.pnmean(('x', 'y', 'z'))

        return ds

    out = get_gradients(out)
    avg["ε_mean"] = avg.ε.pnmean(('z'))
    
    with ProgressBar():
        out.dbds_mean.load()
        out.ε_mean.load()

    out["dbds_mean"] = (out.dbds_mean - out.dbds_mean.min())/out.dbds_mean.max()
    avg["ε_mean"] /= avg.ε_mean.max()
    #----

    #++++ Plot
    out.dbds_mean.plot(label="dbds")
    avg.ε_mean.plot(label="ε")
    plt.legend(); plt.grid()
    plt.show()
    #----

#++++ Make it pretty
for ax in axesf:
    ax.set_xlim(0, 12)
    ax.set_ylabel("")

    ax.legend()
    ax.grid(True)
#----

#-----
fig.savefig(f"figures_paper/dissip_evolution.pdf")
#-----

