import numpy as np
import xarray as xr
from cmocean import cm
from matplotlib import pyplot as plt
from aux00_utils import open_simulation
import pynanigans as pn
plt.rcParams['figure.constrained_layout.use'] = True


# Define dir and file names
#++++
dirname = "ISI_jet"
#dirname = "testing_ISI_jet"
#dirname = "testing_ISI_jet2"
#dirname = "test_profile"
#path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
path = f"/glade/u/home/tomasc/work_cheyenne/{dirname}/data/"
snames = [#"FNN_CIintjet01",
          "FNN_CIsurfjet1",
          "FNN_CIsurfjet3",
          "FNN_SIsurfjet4",
          #"FNN_TEST3SIsurfjet4",
          #"FNN_sloshsurfjet1",
          #"FNN_stabsurfjet1",
          #"PNN_CIintjet01",
          #"PNN_CIintjet2",
          #"PNN_CIsurfjet1",
          #"PNN_CIsurfjet2",
          #"PNN_CIsurfjet3",
          #"PNN_CIsurfjet4",
          #"PNN_CIsurfjet5",
          #"PNN_SIsurfjet1",
          #"PNN_SIsurfjet2",
          #"PNN_SIsurfjet3",
          #"PNN_SIsurfjet4",
          #"PNN_SIsurfjet5",
          #"PNN_SIsurfjet6",
          ]
#----



#+++++ Definte plotting function
def plot_horvort(vid, fig, tt, framedim="time", avg=None, **kwargs):

    #++++ Create subplots and get correct time
    plt.rcParams['font.size'] = 8
    axes = fig.subplots(ncols=2, nrows=3, 
                        gridspec_kw=dict(height_ratios=[0.4, 0.4, 0.2]),
                        )
    axesf = axes.flatten()

    vid0 = vid.isel(time=tt)
    #----

    #++++ First panel
    lim = 0.02
    vid0.tke.pnimshow(ax=axesf[0], cmap="inferno", y='z',
                                    vmin=0, vmax=lim, 
                                    **kwargs)
    lim = abs(vid0.u_0)
    vid0.u_tot.plot.contour(ax=axesf[0], colors="w", linewidths=0.2,
                                         y='zC', levels=np.linspace(-lim, lim, 30),
                                         )
    #---- 

    #++++ Second panel
    lim = 5e-9
    vid0.PV.pnimshow(ax=axesf[1], cmap=cm.balance, 
                              y='z',
                            vmin=-lim, vmax=lim, 
                            **kwargs)
    vid0.PV.plot.contour(ax=axesf[1], colors="k", linewidths=0.1,
                                      levels=[0])
    #---- 

    #++++ Third panel
    lim = 1e-8
    vid0.ε.pnimshow(ax=axesf[2], cmap="afmhot",
                                   y='z',
                                   vmin=0, vmax=lim, 
                                   **kwargs)
    #---- 

    #++++ Fourth panel
    lim = 1e-2
    levels = np.linspace(vid0.b_tot.min(), vid0.b_tot.max(), 20)
    vid0.ω_x.plot.imshow(ax=axesf[3], cmap=cm.balance,
                                     y='zF',
                                     vmin=-lim, vmax=lim, 
                                     **kwargs)
    vid0.b_tot.plot.contour(ax=axesf[3], colors="k", linewidths=0.2,
                                    y='zC', levels=levels,
                                    )
    #---- 

    #++++ Fifth panel
    if avg is not None:
        avg0 = avg.sel(time=slice(0, vid0.time))
        avg0.ε.plot(ax=axesf[4], color='r', label="ε")
        avg0.dBPEdt.plot(ax=axesf[4], color='b', label=r"$\epsilon_p$ (dBPE/dt)")
        avg0.sponge_dissip.plot(ax=axesf[4], color='g', label="$\epsilon_s$ (sponge)", ls='--')
    #----

    #++++ Sixth panel
    γ = avg0.dBPEdt / (avg0.dBPEdt + avg0.ε)
    γ_sp = avg0.dBPEdt / (avg0.dBPEdt + avg0.ε + avg0.sponge_dissip)
    γ_sp.plot(ax=axesf[5], color='g', label="$\gamma_{sp}$ (w/ sponge)", ls='--')
    γ.plot(ax=axesf[5], color='k', label="$\gamma$", ls='-')
    axesf[5].set_ylim(-0.05, 1.1)
    #----

    #++++ Fix axes
    for i, ax in enumerate(axesf):
        if i==0:
            time_sec = float(vid0.time*vid.T_inertial)
            ax.set_title(f"Time = {time_sec/3600:.2f} hours, {vid0.time.item():.2f} Inertial periods")
        else:
            ax.set_title('')
        for ax in axes[:,1]:
            ax.set_ylabel('')
        for ax in axes[0]:
            ax.set_xlabel('')

        ax.set_xlim(vid0.yF.min(), vid0.yC.max())
    for ax in axes[-1]:
        ax.legend()
        ax.grid()
    #----

    fig.set_constrained_layout_pads(w_pad=0, h_pad=0, hspace=0, wspace=0)
    return axes, fig
#------



for sname in snames:

    #++++ Open dataset
    if __name__ == "__main__": print(f"\nOpening {sname}\n")
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(),
                                    )
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    grid_out, out = open_simulation(path+f"out.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    load=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #------
    
    #++++ Downsample and chunk to save space
    vid = vid.sel(time=slice(None, 10*vid.T_inertial))
    vid = pn.downsample(vid, yC=1000, yF=1000, zC=500, zF=500, time=150, round_func=np.ceil)
    vid = vid.pnchunk(maxsize_4d=1000**2, round_func=np.ceil)
    #----
    
    
    #++++ Include avg quantities in vid for plotting
    avg["BPE"] = -avg.b_sorted * avg.zC
    avg_0d = avg.mean("zC")
    avg_0d["dBPEdt"] = avg_0d.BPE.differentiate("time") / avg.T_inertial
    #----
    
    
    #++++ Begin plotting
    kwargs = dict(avg=avg_0d, interpolation="none", 
                  add_colorbar=True, 
                  cbar_kwargs=dict(extend='both',
                  shrink=0.8),
                  rasterized=True,
                  )
    
    if 0:
        fig = plt.figure()
        plot_horvort(vid, fig, tt=20, **kwargs 
                     )
    else:
        from xmovie import Movie
        anim_horvort = Movie(vid, plotfunc=plot_horvort,
                             framedim='time',
                             pixelwidth=1800,
                             pixelheight=1000,
                             dpi=200,
                             frame_pattern='frame_%05d.png',
                             fieldname=None,
                             input_check=False,
                             **kwargs
                     )
    
        # The condition below is necessary for processes scheduler: 
        # https://stackoverflow.com/questions/18204782/runtimeerror-on-windows-trying-python-multiprocessing
        if __name__ == '__main__':
            from dask.distributed import Client
    
            #client = Client(n_workers=18, memory_limit='0.5GB', processes=False)
            #print("Client :", client)
            print("Start frame saving")
            anim_horvort.save(f"anims2/{sname}_horvort.mp4",
                              remove_frames=True,
                              remove_movie=False,
                              progress=True,
                              verbose=False,
                              overwrite_existing=True,
                              framerate=16,
                              parallel=True,
                              #parallel_compute_kwargs=dict(num_workers=18, memory_limit='5GB', processes=False), # 24 min
                              parallel_compute_kwargs=dict(num_workers=18, memory_limit='5GB', scheduler="processes"), # 2 min
                              )
            plt.close("all")
    #-----
    
