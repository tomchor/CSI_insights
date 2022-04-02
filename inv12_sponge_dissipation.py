import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIsurfjet1",
          "PNN_CIsurfjet2",
          "PNN_CIsurfjet3",
          "PNN_SIsurfjet1",
          "PNN_SIsurfjet2",
          "PNN_SIsurfjet3",
          "PNN_SIsurfjet4",
          "PNN_SIsurfjet5",
          "PNN_SIsurfjet6",
          "PNN_CIintjet01",
          ]
#----


for sname in snames:
    #++++ Open datasets avg and vid
    print(f"Opening {sname}")
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #----


 
    #++++ Get important times and perform average
    times_vid = [0, 0.5, 2, 3, 5, 10]
    times_avg = slice(None, None, None)
    
    avg = avg.sel(time=times_avg)
    avg_0d = avg.mean(("zC"))
    #-----
    
    #++++ Integrate dissipations in time
    ε_k_int = avg_0d.ε.integrate("time")
    ε_s_int = avg_0d.sponge_dissip.integrate("time")
    #----

    #++++ Ratio between dissipations
    print(f"Ratio between dissipations for simulation {sname} is ", float(ε_s_int / ε_k_int))
    #----
