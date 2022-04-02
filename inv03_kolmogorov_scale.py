import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
π = np.pi


#++++ Define directory and simulation name
dirname = "ISI_jet"
#dirname = "testing_ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = [#"FNN_CIintjet01",
          "FNN_CIsurfjet1",
          "FNN_CIsurfjet3",
          "FNN_SIsurfjet4",
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
    yslice = slice(1/8*vid.Ly, 7/8*vid.Ly)
    avg = avg.sel(yC=yslice, yF=yslice)
    vid = vid.sel(yC=yslice, yF=yslice)

    avg = avg.pnchunk(maxsize_4d=1000**2)
    vid = vid.pnchunk(maxsize_4d=1000**2, round_func=np.ceil)
    #-----


    #++++ Get times for minimum and maximum in dissipation
    avg_0d = avg.mean(("zC", "zF"))
    it_min = int(avg_0d.ε.sel(time=slice(0, 5)).argmin("time"))
    it_max = int(avg_0d.ε.argmax("time"))
    t_min, t_max = avg_0d.isel(time=[it_min, it_max]).time
    #----

    # Get important times and perform average
    #++++
    times_vid = np.sort(np.append(t_max, [0, 0.5, 2, 5, 10]))
    times_avg = slice(None, None, None)

    avg = avg.sel(time=times_avg)
    vid = vid.sel(time=times_vid, method="nearest")
    #-----


    # Kolmogorov scale calculation
    #++++++++
    if False:
        print("Calculating ε_max by percentile")
        ε_max = np.percentile(vid.ε.sel(time=t_max, method="nearest"), 99.)
    else:
        print("Calculating ε_max by coarsening")
        y_window = int(100//vid.ΔyF.max()) # 100 m resolution in horizontal
        z_window = int(10//vid.ΔzF.max()) # 10 m resolution in vertical

        ε_max = vid.ε.coarsen(yC=y_window, zC=z_window, boundary="trim").mean().max()
    ν = np.array([vid.νh, vid.νz])
    η = (ν**3/float(ε_max))**(1/4)
    Δξ = [float(vid.ΔyC.max()), float(vid.ΔzF.max())]
    print(f"Should be >1: {η/Δξ}\n")
    #------


    #+++++ Ozmidov scale calculation
    if False:
        N2 = grid_vid.interp(vid.dbdz, "z")
        N2 = np.maximum(N2, 0) # Get rid of negative values
        L_o = np.sqrt(vid.ε / N2**(3/2))
    #-----

