import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, prettify_names
from aux02_plotting import letterize
from matplotlib import pyplot as plt
from cmocean import cm


#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data_post/"
snames = ["PNN_CIfront1",
          "PNN_SIfront5",
          ]
#----


for i, sname in enumerate(snames):
    #++++ Figure parames (ncols, nrows, vars)
    ncols=1
    nrows=1
    size=4.5
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*size, nrows*size/2),
                             squeeze=False, constrained_layout=True,
                             sharex=True, sharey="row")
    #-----

    #++++ Open datasets avg and vid
    print(sname)
    ds_eff = xr.load_dataset(f"data_post/efficiencies_{sname}.nc", decode_times=False)
    #-----

    #+++++ Prepare for plotting
    ds_eff = prettify_names(ds_eff)
    ds_eff = ds_eff.squeeze()

    ε_max = ds_eff.ε.max()
    ds_eff["ε"] = ds_eff.ε / ε_max
    ds_eff["dBPEdt"] = ds_eff.dBPEdt / ε_max
    #-----

    #+++++ Plot panels
    ds_eff.γ.plot.line(ax=axes[0, 0], c="k", ls="--", label=r"Mixing efficiency ($\gamma$)")
    #----

    #++++ Make it pretty
    ax = axes[0, 0]
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 0.4)
    ax.set_ylabel("")
    ax.set_title("")

    ax.legend()
    ax.grid(True)

    fig.savefig(f"figures_talks/ME_evolution_{sname}.pdf")
#-----

