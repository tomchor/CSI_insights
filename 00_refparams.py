import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from aux01_physfuncs import get_ref_quantities, adjust_variables, get_bulk_parameters, get_relevant_attrs
π = np.pi

#++++ Define directory and simulation name
path = f"simulations/data/"
snames = ["PNN_CIfront1",
          "PNN_CIfront2",
          "PNN_CIfront3",
          "PNN_CIfront4",
          "PNN_CIfront5",
          "PNN_SIfront1",
          "PNN_SIfront2",
          "PNN_SIfront3",
          "PNN_SIfront4",
          "PNN_SIfront5",
          "PNN_SIfront6",
          #"PNN_CIintjet01",
          "FNN_CIfront1",
          "FNN_CIfront3",
          "FNN_SIfront4",
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
    if "front" in sname:
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

