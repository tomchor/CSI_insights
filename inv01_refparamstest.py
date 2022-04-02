import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from aux01_physfuncs import (get_ref_quantities, adjust_variables, 
                             get_bulk_parameters, get_relevant_attrs,
                             get_IC_jet, get_IC_front,)
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIsurfjet1",
          #"PNN_CIsurfjet2",
          #"PNN_CIsurfjet3",
          #"PNN_CIsurfjet4",
          #"PNN_CIsurfjet5",
          #"PNN_SIsurfjet1",
          #"PNN_SIsurfjet2",
          "PNN_SIsurfjet3",
          #"PNN_SIsurfjet4",
          #"PNN_SIsurfjet5",
          #"PNN_SIsurfjet6",
          "PNN_CIintjet01",
          #"FNN_CIsurfjet1",
          #"FNN_SIsurfjet4",
          ]
#----

paramlist = []
for sname in snames:
    #++++ Open datasets avg and avg
    print(sname)
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(),
                                    )
    #-----

    #+++++ Get parameters and create simulation dimension
    avg = adjust_variables(avg)
    if sname=="PNN_CIintjet01":
        print(f"N2_inf was {float(avg.N2_inf)}")
        #avg.attrs["N2_inf"] = 4.5e-5
        print(f"N2_inf is {float(avg.N2_inf)}")

        print(f"σ_z was {float(avg.σ_z)}")
        avg.attrs["σ_z"] = 120
        print(f"σ_z is {float(avg.σ_z)}")

        print(f"σ_y was {float(avg.σ_y)}")
        #avg.attrs["σ_y"] = 600
        print(f"σ_y is {float(avg.σ_y)}")

        print(f"u_0 was {float(avg.u_0)}")
        avg.attrs["u_0"] = -0.3
        print(f"u_0 is {float(avg.u_0)}")

        print(f"f_0 was {float(avg.f_0)}")
        #avg.attrs["f_0"] = .7e-4
        print(f"f_0 is {float(avg.f_0)}")

    params_sim = get_relevant_attrs(avg)
    params_bulk = get_bulk_parameters(avg)
    if "surfjet" in sname:
        setup="front"
        IC = get_IC_front(y=avg.yF, z=avg.zF, 
                          y_0=avg.y_0, z_0=avg.z_0, 
                          σ_y=avg.σ_y, σ_z=avg.σ_z, 
                          u_0=avg.u_0, N2_inf=avg.N2_inf, f_0=
                          avg.f_0, Lz=avg.Lz)
    elif "intjet" in sname:
        setup="jet"
        IC = get_IC_jet(y=avg.yF, z=avg.zF, 
                        y_0=avg.y_0, z_0=avg.z_0, 
                        σ_y=avg.σ_y, σ_z=avg.σ_z, 
                        u_0=avg.u_0, N2_inf=avg.N2_inf, f_0=
                        avg.f_0, Lz=avg.Lz)
    else:
        raise NameError
    params_ref = get_ref_quantities(avg, setup=setup)
    RoRi_neg = -params_ref["Ro_qmin"] * params_ref["Ri_qmin"]
    print(f"-Ri Ro = {RoRi_neg}")

    #++++ Check stratification
    if (IC.b.differentiate('zF')<0).any():
        raise ValueError("Negative N2")
    #----

    params = xr.Dataset(params_sim | params_ref | params_bulk)
    params = params.expand_dims("simulation").assign_coords(simulation=[sname])
    #-----

    #++++ Append reference quantities to list
    paramlist.append(params)
    #----
    
#++++ Create and save dataset
allparams = xr.concat(paramlist, dim="simulation")
allparams.to_netcdf("data/allparamstest.nc")
#----

