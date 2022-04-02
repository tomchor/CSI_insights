import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from aux01_physfuncs import adjust_variables
from dask.diagnostics import ProgressBar
from matplotlib import pyplot as plt
π = np.pi

#++++ Computation options
import dask
if __name__ == '__main__':
    compute_flags = dict(num_workers=18, memory_limit='5GB') # Processes doesn't work very well
else:
    compute_flags = dict()

extra = ""
#----

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = [#"PNN_CIsurfjet1",
          #"PNN_CIsurfjet2",
          #"PNN_CIsurfjet3",
          #"PNN_CIsurfjet4",
          #"PNN_CIsurfjet5",
          "PNN_SIsurfjet1",
          #"PNN_SIsurfjet2",
          #"PNN_SIsurfjet3",
          #"PNN_SIsurfjet4",
          #"PNN_SIsurfjet5",
          "PNN_SIsurfjet6",
          #"PNN_CIintjet01",
          ]
#----

nrows=2; ncols=2
size = 4

Rolim = 1.5
Rilim = 1.
ν_m = 1e-6
allparams = xr.load_dataset("data/allparams.nc")
for i, sname in enumerate(snames):
    print(f"Opening {sname}")
    print(f"Extra: {extra}")

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
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #----

    #++++ Get subdomain
    print("Getting a σ_y-wide slice of the domain")
    out = adjust_variables(out)
    if "surfjet" in sname:
        y_r = allparams.sel(simulation=sname).y_qmin
    else:
        y_r = 1/2 * np.sqrt(2) * out.σ_y + out.y_0
    side_orig = 1/2*out.σ_y
    side = out.σ_y
    yF_r = out.sel(yF=y_r+side, method="nearest").yF
    yF_l = out.sel(yF=y_r-side, method="nearest").yF
    yslice = slice(yF_l, yF_r) # Make sure len(yF) = len(yC)+1

    if "halfz" in extra:
        depth= 40
    else:
        depth=np.inf
    zF_r = out.sel(zF=out.z_0+depth, method="nearest").zF
    zF_l = out.sel(zF=out.z_0-depth, method="nearest").zF
    zslice = slice(zF_l, zF_r) # Make sure len(zF) = len(zC)+1

    tslice = 0

    out = out.sel(yC=yslice, yF=yslice, zC=zslice, zF=zslice, time=tslice)
    subgrid_out = pn.get_grid(out, topology="PNN")
    #-----

    #++++ Get conditional avg mask
    if "cond" in extra:
        print("Getting masks for conditional averaging")
        ε_c = 1e-10
        mask_ccc = out.ε > ε_c
        mask_ccf = subgrid_out.interp(out.ε, 'z', boundary="extrapolate") > ε_c
        mask_cfc = subgrid_out.interp(out.ε, 'y', boundary="extrapolate") > ε_c
        mask_fcc = subgrid_out.interp(out.ε, 'x', boundary="extrapolate") > ε_c

        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            out["ε"] = out.ε.where(mask_ccc, drop=True)
            out["ν_e"] = out.ν_e.where(mask_ccc, drop=True)
            out["dbdz"] = out.dbdz.where(mask_ccf, drop=True)
            out["u"] = out.u.where(mask_fcc, drop=True)
            out["v"] = out.v.where(mask_cfc, drop=True)
            out["w"] = out.w.where(mask_ccf, drop=True)
    #----

    #++++ Calculations
    print("Calculating...")
    def calc_variables(ds, grid=None):
        if grid is None: 
            grid = pn.get_grid(ds, topology="PNN")
        ds["ε_mean"] = ds.ε.pnmean(('x',))

        ds["dvdx"] = grid.derivative(ds.v, "x")
        ds["dudy"] = grid.derivative(ds.u, "y", boundary="extrapolate")
        ds["Ro"] = (ds.dvdx - ds.dudy)/ds.f_0
        ds["Ro_mean"] = ds.Ro.pnmean(('x',))

        ds["dudz2"] = grid.interp(grid.derivative(ds.u, "z", boundary="extrapolate")**2, "x")
        ds["dvdz2"] = grid.interp(grid.derivative(ds.v, "z", boundary="extrapolate")**2, "y")
        ds["S2"] = (ds.dudz2 + ds.dvdz2)
        ds["Ri"] = (ds.dbdz / ds.S2)
        ds["Ri_mean"] = ds.Ri.pnmean(('x',))

        ds["u_mean"] = ds.u.pnmean(('x',))
        ds["PV_mean"] = ds.PV.pnmean(('x',))

        return ds


    out = calc_variables(out, subgrid_out)
    #----
    
    #++++ Create figure
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.3*ncols*size, nrows*size),
                                     constrained_layout=True, squeeze=False)
    axesf = axes.flatten()
    #-----

    #+++++ Plot
    with ProgressBar():
        out.Ri_mean.plot.imshow(ax=axesf[0], rasterized=True, vmin=-Rilim, vmax=+Rolim, cmap="RdBu_r")
        out.Ro_mean.plot.imshow(ax=axesf[1], rasterized=True, vmin=-Rolim, vmax=+Rolim, cmap="RdBu_r")
        out.u_mean.plot.imshow(ax=axesf[2], rasterized=True, vmin=out.u_0, vmax=-out.u_0, cmap="RdBu_r")
        out.PV_mean.plot.imshow(ax=axesf[3], rasterized=True, cmap="RdBu_r")
    #----

    #++++ Get relevant lines
    for ax in axesf:
        ax.axvline(x=y_r, ls=":", c="w")
        ax.axvline(x=y_r-side_orig, ls="--", c="k")
        ax.axvline(x=y_r+side_orig, ls="--", c="k")
    #----

    #++++ Save figure
    fig.savefig(f"figures_check/inv13_Ro_{sname}.png")
    #----

