import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, prettify_names
from aux02_plotting import letterize
from matplotlib import pyplot as plt
from cmocean import cm


#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIfront1",
          "PNN_SIfront4",
          ]
#----

#++++ Figure parames (ncols, nrows, vars)
ncols=len(snames)
nrows=2
size=4.5

fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*size, nrows*size/2),
                         squeeze=False, constrained_layout=True,
                         sharex=True, sharey="row")
axesf = axes.flatten()
#-----

for i, sname in enumerate(snames):
    #++++ Open datasets avg and vid
    print(sname)
    #grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
    #                                use_inertial_periods=True,
    #                                topology=sname[:3],
    #                                squeeze=True,
    #                                load=False,
    #                                open_dataset_kwargs=dict(chunks=dict(time=1)),
    #                                )
    ds_eff = xr.load_dataset(f"data/efficiencies_{sname}.nc", decode_times=False)
    #-----

    #+++++ Prepare for plotting
    ds_eff = prettify_names(ds_eff)
    ds_eff = ds_eff.squeeze()

    ε_max = ds_eff.ε.max()
    ds_eff["ε"] = ds_eff.ε / ε_max
    ds_eff["dBPEdt"] = ds_eff.dBPEdt / ε_max
    #-----

    #+++++ Plot panels
    ds_eff.ε.plot.line(ax=axes[0, i], c='r', ls="-", label=r"$\langle\epsilon_k\rangle$ [normalized]")
    ds_eff.dBPEdt.plot.line(ax=axes[0, i], c='b', ls="-", label=r"$\langle\epsilon_p\rangle$ [normalized]")

    ds_eff.γ.plot.line(ax=axes[1, i], c="k", ls="--", label=r"$\gamma$")
    ds_eff.Γ.plot.line(ax=axes[1, i], c="k", ls="-", label=r"$\Gamma$")
    #----

#++++ Make it pretty
for ax in axesf:
    ax.set_xlim(0, 12)
    ax.set_ylabel("")

    ax.legend()
    ax.grid(True)
#----

#-----
letterize(axesf, 0.03, 0.9, fontsize=14)
fig.savefig(f"figures_paper/dissip_evolution.pdf")
#-----

