import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, pnames
from aux02_plotting import letterize
from matplotlib import pyplot as plt
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = [#"FNN_CIintjet01",
          "FNN_CIsurfjet1",
          "FNN_CIsurfjet3",
          "FNN_SIsurfjet4",
          ]
#----

#++++ Options
debug = False
#----



#++++ Adjust colormap to shade unstable regions
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

RdBu_r = cm.get_cmap('RdBu_r', 256)
newcolors = RdBu_r(np.linspace(0, 1, 256))
gray = np.array([180/256, 180/256, 180/256, 1])
newcolors[ 256//2 : int(256*(1/2+1/8)), :] = gray
newcmp = ListedColormap(newcolors)
#----


#++++ Plotting
print("Starting to plot")
if debug:
    jump=3
else:
    jump=5

vidwind = xr.load_dataset(f"data/shearprod_secondary_{snames[0]}.nc")
ncols = len(vidwind.time.values[::jump])
if debug:
    nrows = 6
else:
    nrows = len(snames)+1
size = 2.
fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(1.2*ncols*size, nrows*size),
                         constrained_layout=True, sharex=False, sharey=True)

gs_avg = axes[-1, 0].get_gridspec()
ax_avgs = []
for joint in range(-3, 0):
    for ax in axes[joint, :].flatten():
        ax.remove()
    ax_avgs.append(fig.add_subplot(gs_avg[joint, :]))
letterize(np.array(fig.axes), 0.05, 0.9, fontsize=14, zorder=10)
#----




for ns, sname in enumerate(snames):
    #++++ Open datasets avg and vid
    print(f"Opening {sname}")
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )

    vidwind = xr.load_dataset(f"data/shearprod_secondary_{sname}.nc")
    vidwind = vidwind.reset_coords()
    #----

    #++++ Plot line plots
    ax_avg = ax_avgs[ns]
    vidwind.wb_2nd_mean.plot(ax=ax_avg, x="time", 
                             label=r"$\langle{w'b'}\rangle_s$", lw=4, marker="o")
    vidwind.SPy_2nd_mean.plot(ax=ax_avg, x="time", 
                             label=r"$\langle SP_h^\mathrm{second}\rangle_s$", lw=4, marker="o")
    vidwind.SPz_2nd_mean.plot(ax=ax_avg, x="time", 
                             label=r"$\langle SP_v^\mathrm{second}\rangle_s$", lw=4, marker="o")
    #----

    #++++ Helpful labels
    pname = pnames[sname]
    if ns==0: ax_avg.legend(fontsize=13, loc="upper left", bbox_to_anchor=(0.1, 0, 1, 1))
    ax_avg.grid(True)
    ax_avg.set_ylabel(f"TKE equation terms\n{pname}", fontsize=12)
    #----

    #++++ Plots snapshots
    if ns == 0:
        print("Plotting snapshots")
        for col, time in enumerate(vidwind.time.values[::jump]):
            if col<ncols-1:
                cbar=False
                cbar_kwargs=None
            else:
                cbar=True
                cbar_kwargs=dict(shrink=0.8)
            vidsnap0 = vidwind.sel(time=time).squeeze()
    
            lim=1
            vidsnap0.Ri.pnimshow(ax=axes[0,col], vmin=-lim, vmax=lim, 
                                 cmap=newcmp, add_colorbar=cbar,
                                 cbar_kwargs=cbar_kwargs,
                                 rasterized=True)
    
            axes[0,col].set_title(f"{time:.2f} inertial periods")
    
            if debug:
                lim=8e-3
                vidsnap0.ω_x.plot.imshow(ax=axes[1,col], vmin=-lim, vmax=lim, 
                                      cmap="RdBu_r", add_colorbar=cbar, rasterized=True)
    
                levs = np.linspace(vidsnap0.b_tot.min(), vidsnap0.b_tot.max(), 20)
                vidsnap0.b_tot.plot.contour(ax=axes[1,col], levels=levs, linewidths=0.1, linestyles="--", colors="k")
    
                lim=1e-7
                vidsnap0.SPy_2nd.plot.imshow(ax=axes[2,col], vmin=-lim, vmax=lim, 
                                        cmap="RdBu_r", add_colorbar=cbar, rasterized=True)
    
                vidsnap0.SPz_2nd.plot.imshow(ax=axes[3,col], vmin=-lim, vmax=lim, 
                                        cmap="RdBu_r", add_colorbar=cbar, rasterized=True)
    
                vidsnap0.wb_2nd.plot.imshow(ax=axes[4,col], vmin=-lim, vmax=lim, 
                                        cmap="RdBu_r", add_colorbar=cbar, rasterized=True)
    
    #----
    
    #++++ Prettify and save
    axes.flatten()[-1].grid()
    if debug:
        fig.savefig(f"figures_paper/2Dinst_snapshot_{snames[0]}_debug.pdf")
    else:
        fig.savefig(f"figures_paper/2Dinst_snapshot_{snames[0]}.pdf")
    #----
