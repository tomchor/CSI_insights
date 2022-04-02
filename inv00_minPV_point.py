import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from cmocean import cm


#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIsurfjet1",
          "PNN_CIsurfjet2",
          "PNN_CIsurfjet3",
          "PNN_CIsurfjet4",
          "PNN_CIsurfjet5",
          "PNN_SIsurfjet1",
          "PNN_SIsurfjet2",
          "PNN_SIsurfjet3",
          "PNN_SIsurfjet4",
          "PNN_SIsurfjet5",
          "PNN_SIsurfjet6",
          ]
#----

#++++ Figure parames (ncols, nrows, vars)
ncols=3
nrows=int(np.ceil(len(snames)/ncols))
size=4.5

fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*size, nrows*size/2.5),
                         squeeze=False, constrained_layout=True,
                         sharex=True, sharey=True)
axesf = axes.flatten()
#-----

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
    #-----

    #+++++ Select relevant times and subdomain
    vid = vid.sel(time=0, method="nearest")
    yslice = slice(2000, 6000)
    zslice = slice(-20, 0)
    vid = vid.sel(yC=yslice, yF=yslice, zC=zslice, zF=zslice)
    vid = pn.downsample(vid, yC=200, yF=200, round_func=np.ceil)
    #-----

    #+++++ Get variables ready
    if "PV" not in vid.variables:
        vid["q_hat"] = (vid.PV_ver + vid.PV_hor) / (vid.f_0 * vid.N2_inf)
    else:
        vid["q_hat"] = vid.PV / (vid.f_0 * vid.N2_inf)
    vid.q_hat.attrs = dict(long_name=r"Normalized PV ($\hat q$)")
    #-----

    #++++ Plot reference point
    axesf[i].scatter(vid.y_r, vid.z_r, zorder=50, c='w', clip_on=False, edgecolor='k', s=40)

    y_qmin = allparams.sel(simulation=sname).y_qmin
    axesf[i].scatter(y_qmin, 0, zorder=50, c='y', clip_on=False, edgecolor='k', s=40)
    #----

    #+++++ Plot panels
    vid.q_hat.pnimshow(ax=axesf[i], cmap="RdBu_r", vmin=None, vmax=None,
                       rasterized=True)
    vid.q_hat.pncontour(ax=axesf[i], colors='k', levels=20)
    axesf[i].set_title(sname)
    #----

#-----
fig.savefig(f"figures_check/overview_min_PV.png")
#-----

