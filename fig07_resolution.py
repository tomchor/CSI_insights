import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, filter_sims
from matplotlib import pyplot as plt
from cmocean import cm
from aux02_plotting import letterize, color_base, marker_base



#++++ Define directory and simulation name
Cnames = ["PNN_CIsurfjet1",
          "PNN_CIsurfjet1_f2",
          "PNN_CIsurfjet1_f4",
          "PNN_CIsurfjet1_f8",
          ]
Snames = ["PNN_SIsurfjet4",
          "PNN_SIsurfjet4_f2",
          "PNN_SIsurfjet4_f4",
          "PNN_SIsurfjet4_f8",
          ]
snames = [Cnames, Snames]
alleffs = xr.load_dataset(f"data_post/alleffs.nc", decode_times=False)
alleffs.Δz.attrs = dict(units="m", long_name=r"$\Delta z$")

#alleffs = filter_sims(alleffs, prettify=True, only_main=True)
#----

#++++ Figure parames (ncols, nrows, vars)
ncols=2
nrows=1
size=4.5

fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*size, nrows*size/1.2),
                         squeeze=False, constrained_layout=True,
                         sharex="row", sharey=False)
axesf = axes.flatten()
#-----



#+++ Plot!
for i, lname in enumerate(snames):
    dsall = alleffs.sel(simulation = lname)
    dsall["Δz/Lo"] = alleffs.Δz / alleffs.Lo
    dsall["Δz/Lo_avg"] = alleffs.Δz / alleffs.Lo_avg
    dsall["Δz/Lo_inf"] = alleffs.Δz / alleffs.Lo_inf

    #dsall.plot.scatter(ax=axesf[0], x="Δz", y="Δz/Lo", color=color_base[i], marker="x", label=lname[0])
    dsall.plot.scatter(ax=axesf[0], x="Δz", y="Δz/Lo_inf", color=color_base[i], marker="o", label=lname[0])
    #dsall.plot.scatter(ax=ax[0], x="Δz", y="Δz/Lo_avg",)

    dsall.plot.scatter(ax=axesf[1], x="Δz", y="Γ_last", color=color_base[i], marker="o", label=lname[0])
#----

#++++ Make it pretty
for ax in axesf:
    ax.set_ylim(0, None)
    ax.grid(True)
    ax.legend(fontsize='x-small', labelspacing=1.1)

axesf[0].set_title("Spacing normalized by the Ozmidov scale")
axesf[1].set_title("Cumulative mixing effieciency")
#----

#-----
letterize(axesf, 0.1, 0.9, fontsize=14)
fig.savefig(f"figures_paper/resolution.pdf")
#-----

