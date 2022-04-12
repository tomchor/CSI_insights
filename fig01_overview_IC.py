import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, regularize_ds
from aux02_plotting import letterize
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse
from cmocean import cm


#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIfront1",
          "PNN_SIfront4",
          ]
#----

#++++ Figure parames (ncols, nrows, vars)
variables = ["q_hat", "ω_x", "ε"]
times = [0, 5, 5]
vmins = [-1.5, -8e-3, 1e-10]
vmaxs = [+1.5, +8e-3, 1e-8]
cmaps = ["RdBu_r", cm.balance, "inferno"]
norms = [None, None, LogNorm()]
bconts = [True, True, False]

ncols=len(snames)
nrows=len(variables)
size=4.5

fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*size, nrows*size/2.5),
                         squeeze=False, constrained_layout=True,
                         sharex=True, sharey=True)
axesf = axes.flatten()
#-----

def fmt(x):
    s = f"{x:.2f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return f"{s} m/s"


allparams = xr.load_dataset("data/allparams.nc")
for i, sname in enumerate(snames):
    #++++ Open datasets avg and vid
    print(sname)
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    vid = regularize_ds(vid, check_PV=True, check_utot=False, check_btot=False)
    #-----

    #+++++ Select relevant times and subdomain
    vid = vid.sel(time=times, method="nearest")
    yslice = slice(1000, 7000)
    vid = vid.sel(yC=yslice, yF=yslice)

    #vid = pn.downsample(vid, yC=200, yF=200, round_func=np.ceil)
    #-----

    #+++++ Get variables ready
    if "q_hat" in variables:
        vid["q_hat"] = vid.PV / (vid.f_0 * vid.N2_inf)
        vid.q_hat.attrs = dict(long_name=r"Normalized PV ($\hat q_b$)")
    vid.ω_x.attrs = dict(long_name=r"x-vorticity", units="1/s")
    vid.ε.attrs = dict(long_name=r"$\varepsilon_k$", units="m$^2$/s$^3$")
    #-----

    #++++ Colorbar?
    if sname == snames[-1]:
        add_cbar=True
    else:
        add_cbar=False
    #----

    #++++ Plot reference point
    y_r = allparams.sel(simulation=sname).y_qmin
    z_r = 0
    axes[0, i].scatter(y_r, z_r, zorder=50, c='w', clip_on=False, edgecolor='k', s=40)
    #----

    #+++++ Plot panels
    for j, (time, var) in enumerate(zip(times, variables)):
        vmin, vmax = vmins[j], vmaxs[j]
        cmap = cmaps[j]
        bcont = bconts[j]
        norm = norms[j]

        vidj = vid.isel(time=j)

        vidj[var].pnimshow(ax=axes[j, i], cmap=cmap, 
                           vmin=vmin, vmax=vmax, norm=norm,
                           add_colorbar=add_cbar, rasterized=True)

        if time==0:
            lim = abs(vid.u_0)
            levels=np.linspace(-lim, lim, 10)
            CS = vidj.u.pncontour(ax=axes[j, i], levels=levels, colors="g", linewidths=0.4, linestyles="solid")
            axes[j, i].clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=6)

        if bcont:
            levels = np.linspace(vidj.b.min(), vidj.b.max(), 20)
            vidj.b.pncontour(ax=axes[j, i], levels=levels, colors="k", linewidths=0.4, linestyles="dashed")

        axes[j, i].set_title(f"Inertial time = {time:.1f}")
    #----

#-----
letterize(axesf, 0.05, 0.9, fontsize=14)
fig.savefig(f"figures_paper/overview_cross_section.pdf")
#-----

