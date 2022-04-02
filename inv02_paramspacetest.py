import numpy as np
from matplotlib import pyplot as plt
import pynanigans as pn
import xarray as xr
from aux00_utils import filter_sims
from aux02_plotting import plot_scatter

#++++ Open and curate datasets
allparams = xr.load_dataset("data/allparamstest.nc")

allparams = filter_sims(allparams, prettify=True, only_main=False)

allparams["Ri_qmininv"] = 1/allparams.Ri_qmin
#----

#++++ Create figure
size=5
ncols=1; nrows=1
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*size, size*nrows*.8),
                         constrained_layout=True,
                         squeeze=False, sharey=False, sharex="row")

axesf = axes[:, ncols-1].flatten()
#-----

#++++ Standardize
def standardize_ax(ax, min_ri_inv=-1, max_ri_inv=6,
                   min_ro=-2, max_ro=1/2):

    #++++ Plot heatmap
    Ri_inv = xr.DataArray(np.linspace(-1, 15), dims=['Ri_inv'])
    Ri_inv = Ri_inv.assign_coords(Ri_inv=Ri_inv)

    Ro_tot = xr.DataArray(np.linspace(min_ro, max_ro), dims=['Ro_tot'])
    Ro_tot = Ro_tot.assign_coords(Ro_tot=Ro_tot)

    q_hat = 1 + Ro_tot - Ri_inv
    q_hat.attrs = dict(long_name=r"$\hat q$", units="-")
    q_hat = q_hat.sel(Ri_inv=slice(min_ri_inv, 1.1*max_ri_inv))
    q_hat.plot.contourf(cmap="Blues_r", zorder=0, vmax=0, 
                        vmin=np.floor(q_hat.min()),
                        add_colorbar=True,
                        cbar_kwargs=dict(shrink=0.9, orientation="horizontal"))
    #----

    #++++ Plot relevant "barriers"
    Ri_inv = xr.DataArray(np.linspace(0, max_ri_inv), dims=['Ri_inv'])
    Ri_inv = Ri_inv.assign_coords(Ri_inv=Ri_inv)

    #ax.plot(Ri_inv, Ri_inv-1, c="k", ls="--", label="$q=0$")
    Ri_inv2 = Ri_inv.sel(Ri_inv=slice(0.5, None))
    ax.plot(Ri_inv2, -Ri_inv2, c="k", ls="-.", label="$Ro_r Ri_r=-1$")

    ax.fill_between(Ri_inv, Ri_inv-1, 10, color="lightgray", zorder=0)

    Ro = xr.DataArray(np.linspace(-3, -1), dims=['Ro'])
    Ro = Ro.assign_coords(Ro=Ro)

    ax.fill_betweenx(Ro_tot, min_ri_inv, 0, color="darkgray", zorder=0)
    #----

    ax.grid(True)
    ax.set_ylim(min_ro, max_ro)
    ax.set_xlim(min_ri_inv, max_ri_inv)

    return

standardize_ax(axesf[0], min_ro=-5., max_ro=0, min_ri_inv=-1, max_ri_inv=7)
#----


#+++ Plot reference parameter data
opts = dict(clip_on=False, zorder=50, s=100)
def plot_refs(ax, easy=False):
    plot_scatter(allparams, ax=ax, x="Ri_qmininv", y="Ro_qmin", hue="simulation", **opts)
    ax.set_xlabel(r"$1/Ri_r$")
    ax.set_ylabel(r"$Ro_r$")
    return
plot_refs(axesf[0])
#----


#++++ Save figure
axesf[0].legend(loc="upper left", bbox_to_anchor=(1., 1), fontsize='x-small', labelspacing=1.1)
#fig.savefig(f"figures_check/paramspace_check.pdf")
#----
