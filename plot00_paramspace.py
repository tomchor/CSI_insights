import numpy as np
from matplotlib import pyplot as plt
import pynanigans as pn
import xarray as xr

#++++ Open and curate datasets
allparams = xr.load_dataset("data/allparams.nc")
alleffs = xr.load_dataset("data/alleffs.nc")

allparams["Ri_rinv"] = 1/allparams.Ri_r
allparams["Ri_bulk_inv"] = 1/allparams.Ri_bulk

alleffs = alleffs.reindex(simulation=allparams.simulation)
allparams = allparams.combine_first(alleffs[["Γ_last", "γ_last"]])

varlist = ["N2_inf", "σ_y", "PV_r", "Γ_last"]
#----

#++++ Create figure
size=4.5
ncols=3; nrows=len(varlist)
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*size, size),
                         constrained_layout=True,
                         squeeze=False, sharey=False, sharex="row")

#++++ Merge left panels into one
def merge_rows(fig, axes, col=0):
    gs = axes[0, col].get_gridspec()
    for ax in axes[:, col]:
        ax.remove()
    ax = fig.add_subplot(gs[:, col])
    return ax
ax0 = merge_rows(fig, axes, col=0)
ax1 = merge_rows(fig, axes, col=1)
#----

axesf = axes[:, ncols-1].flatten()
#-----

#++++ Standardize
def standardize_ax(ax, min_ri_inv=-1, max_ri_inv=6,
                   min_ro=-2, max_ro=1/2):

    Ri_inv = xr.DataArray(np.linspace(-1, 15), dims=['Ri_inv'])
    Ri_inv.assign_coords(Ri_inv=Ri_inv)

    Ro = xr.DataArray(np.linspace(-3, -1), dims=['Ro'])
    Ro.assign_coords(Ro=Ro)

    ax.plot(Ri_inv, Ri_inv-1, c="k", ls="--", label="$q=0$")
    ax.plot(Ri_inv, -Ri_inv, c="k", ls="-.", label="$Ro_r Ri_r=-1$")
    ax.axhline(y=0, c="k", ls=":")
    ax.axvline(x=0, c="k", ls=":")
    ax.fill_between(Ri_inv, Ri_inv-1, 10, color="lightgray", zorder=0)
    ax.fill_betweenx(Ro, Ro+1, 0, color="darkgray", zorder=0)

    ax.grid(True)
    ax.set_ylim(min_ro, max_ro)
    ax.set_xlim(min_ri_inv, max_ri_inv)

    return

standardize_ax(ax0, min_ro=-3, max_ro=0.5, min_ri_inv=-1, max_ri_inv=7)
standardize_ax(ax1, min_ro=-3, max_ro=0.5, min_ri_inv=-1, max_ri_inv=7)
#----


#+++ Plot reference parameter data
def plot_refs(ax):
    allparams.plot.scatter(ax=ax, x="Ri_rinv", y="Ro_r", hue="simulation", clip_on=False)

    #++++ Make it pretty
    ax.set_xlabel(r"$1/Ri_r$")
    ax.set_ylabel(r"$Ro_r$")
    #----
    return
plot_refs(ax0)
#----

#+++ Plot bulk parameter data
def plot_refs(ax, **kwargs):
    allparams.plot.scatter(ax=ax, x="Ri_bulk_inv", y="Ro_bulk", hue="simulation", clip_on=False, **kwargs)

    #++++ Make it pretty
    ax.set_xlabel(r"$1/Ri_{bulk}$")
    ax.set_ylabel(r"$Ro_{bulk}$")
    #----
    return
plot_refs(ax1, add_guide=False)
#----


#+++ Plot domain set-up information
xsim = range(len(allparams.simulation))
for ax, var in zip(axesf, varlist):
    ax.bar(xsim, allparams[var])
    ax.set_xticks(xsim)
    if ax != axesf[-1]:
        ax.set_xticklabels([])
    else:
        ax.set_xticklabels(allparams.simulation.values, rotation=20)

    if var=="N2_inf":
        ax.set_yscale("log")
    elif var=="PV_r":
        ax.set_ylim(None, 0)
    elif var=="Γ":
        ax.set_ylim(0, 0.8)
    ax.set_title(var)
    ax.grid()
#----


#++++ Save figure
ax0.legend(loc="upper left")
fig.savefig(f"figures_check/paramspace.pdf")
#----
