import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
π = np.pi


#++++ Define directory and simulation name
dirname = "ISI_jet"
#dirname = "testing_ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIsurfjet1",
          #"PNN_CIsurfjet1_f2",
          #"PNN_CIsurfjet1_f4",
          #"PNN_CIsurfjet1_f8",
          #"PNN_SIsurfjet4",
          #"PNN_SIsurfjet4_f2",
          #"PNN_SIsurfjet4_f4",
          #"PNN_SIsurfjet4_f8",
          ]
#----

for sname in snames:
    #++++ Open datasets avg and vid
    print(f"Opening {sname}")
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
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
    #----


    #+++++ Exclude sponge layer if it isn't excluded already
    yslice = slice(2/8*vid.Ly, 6/8*vid.Ly)
    avg = avg.sel(yC=yslice, yF=yslice)
    vid = vid.sel(yC=yslice, yF=yslice)

    avg = avg.pnchunk(maxsize_4d=1000**2)
    vid = vid.pnchunk(maxsize_4d=1000**2, round_func=np.ceil)
    #-----


    #++++ Get times for minimum and maximum in dissipation
    avg_0d = avg.mean(("zC", "zF"))
    it_max = int(avg_0d.ε.argmax("time"))
    t_max = avg_0d.isel(time=it_max).time
    #----

    # Get important times and perform average
    #++++
    times_avg = slice(None, None, None)
    avg = avg.sel(time=times_avg)

    def_times = [0, 0.5, 2, 5, 10]
    times_vid = np.sort(np.append(t_max,  [el for el in def_times if el > t_max ]))
    vid = vid.sel(time=times_vid, method="nearest")
    #-----


    #++++++ Max(<ε>) calculation
    if False:
        print("Calculating ε_max by percentile")
        ε_max = np.percentile(vid.ε.sel(time=t_max, method="nearest"), 99.)

    else:
        print("Calculating ε_max by coarsening")
        y_window = int(100//vid.ΔyF.max()) # 100 m resolution in horizontal
        z_window = int(10//vid.ΔzF.max()) # 10 m resolution in vertical

        ε_max = vid.ε.coarsen(yC=y_window, zC=z_window, boundary="trim").mean().max()
    #------

    #++++ Filter out places with low dissipation
    vid_filt = vid.where(vid.ε >= 1e-10, drop=True)
    #----

    #+++++ Ozmidov scale calculation
    vid_filt["N"] = np.sqrt(vid_filt.N2_inf)
    Lo_avg = 2*π*np.sqrt(vid_filt.ε / vid_filt.N**3)
    Lb_avg = 2*π*np.sqrt(vid_filt.tke) / vid_filt.N
    #-----

