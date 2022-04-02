import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from aux01_physfuncs import get_ref_quantities, adjust_variables, get_bulk_parameters, get_relevant_attrs
π = np.pi

#++++ Define directory and simulation name
path = f"simulations/data/"
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
          #"PNN_CIintjet01",
          "FNN_CIsurfjet1",
          "FNN_CIsurfjet3",
          "FNN_SIsurfjet4",
          ]
#----

paramlist = []
for sname in snames:
    #++++ Open datasets avg and vid
    print(sname)
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=True,
                                    open_dataset_kwargs=dict(),
                                    )
    #-----

    #+++++ Get parameters and create simulation dimension
    avg = adjust_variables(avg)
    params_sim = get_relevant_attrs(avg)
    params_bulk = get_bulk_parameters(avg)
    if "surfjet" in sname:
        setup="front"
    elif "intjet" in sname:
        setup="jet"
    else:
        raise NameError
    params_ref = get_ref_quantities(avg, setup=setup)

    params = xr.Dataset(params_sim | params_ref | params_bulk)
    params = params.expand_dims("simulation").assign_coords(simulation=[sname])
    #-----

    #++++ Append reference quantities to list
    paramlist.append(params)
    #----
    
#++++ Create and save dataset
allparams = xr.concat(paramlist, dim="simulation")
allparams.to_netcdf("data_post/allparams.nc")
#----

