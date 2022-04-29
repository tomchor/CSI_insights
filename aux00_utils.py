import xarray as xr
import pynanigans as pn

#++++ Create PV, u_tot and b_tot if they don't already exist
def regularize_ds(ds, 
                  check_utot=True, 
                  check_btot=True,
                  check_PV=True,
                  ):
    """ Make datasets more regular (i.e., have the same variables, etc.) """
    if check_utot:
        if "u_tot" not in ds.variables.keys():
            ds["u_tot"] = ds.u

    if check_btot:
        if "b_tot" not in ds.variables.keys():
            ds["b_tot"] = ds.b

    if check_PV:
        if "PV" not in ds.variables.keys():
            ds["PV"] = ds.PV_ver + ds.PV_hor
    
    return ds
#----

#++++ Open simulation following the standard way
def open_simulation(fname, use_inertial_periods=True, 
                    open_dataset_kwargs=dict(),
                    regularize_ds_kwargs=dict(),
                    load=False,
                    squeeze=True,
                    unique=True,
                    verbose=False,
                    topology="PPN", **kwargs):
    if verbose: print(sname, "\n")
    
    #++++ Open dataset and create grid before squeezing
    if load:
        ds = xr.load_dataset(fname, decode_times=False, **open_dataset_kwargs)
    else:
        ds = xr.open_dataset(fname, decode_times=False, **open_dataset_kwargs)
    grid_ds = pn.get_grid(ds, topology=topology, **kwargs)
    #----

    #++++ Squeeze?
    if squeeze: ds = ds.squeeze()
    #----

    #++++ Normalize units and regularize
    if use_inertial_periods:
        ds = pn.normalize_time_by(ds, seconds=ds.T_inertial, new_units="Inertial period")
    ds = regularize_ds(ds, **regularize_ds_kwargs)
    #----

    #++++ Returning only unique times:
    if unique:
        import numpy as np
        _, index = np.unique(ds['time'], return_index=True)
        if verbose and (len(index)!=len(ds.time)): print("Cleaning non-unique indices")
        ds = ds.isel(time=index)
    #----

    return grid_ds, ds
#----

#++++ Simulation names 
pnames = dict(PNN_CIfront1="CIfront1",
              PNN_CIfront2="CIfront2",
              PNN_CIfront3="CIfront3",
              PNN_CIfront4="CIfront4",
              PNN_CIfront5="CIfront5",
              PNN_CIfront6="CIfront6",
              PNN_CIfront7="CIfront7",
              PNN_SIfront1="SIfront1",
              PNN_SIfront2="SIfront2",
              PNN_SIfront3="SIfront3",
              PNN_SIfront4="SIfront4",
              PNN_SIfront5="SIfront5",
              PNN_SIfront6="SIfront6",
              PNN_SIfront7="SIfront7",
              PNN_SIfront8="SIfront8",
              PNN_CIintjet01="CIintjet1",
              FNN_CIfront1="2D_CIfront1",
              FNN_CIfront3="2D_CIfront3",
              FNN_SIfront4="2D_SIfront4",
              FNN_CIintjet01="2D_CIintjet1",
              PNN_CIfront1_AMD="CIfront1_AMD",
              PNN_SIfront4_AMD="SIfront4_AMD",
              )

main_sims = ["PNN_CIfront1",
             "PNN_CIfront2",
             "PNN_CIfront3",
             "PNN_CIfront4",
             "PNN_CIfront5",
             "PNN_CIfront6",
             "PNN_CIfront7",
             "PNN_SIfront1",
             "PNN_SIfront2",
             "PNN_SIfront3",
             "PNN_SIfront4",
             "PNN_SIfront5",
             "PNN_SIfront6",
             "PNN_SIfront7",
             "PNN_SIfront8",
             ]
#-----


#++++ Make names match the paper
def prettify_names(ds):
    return ds.assign_coords(simulation=[ pnames[sim] for sim in ds.simulation.values ])
#----

#++++ Filter main
def filter_sims(ds, prettify=True, only_main=False):
    if only_main:
        simulations = list(filter(lambda x: x in main_sims, ds.simulation.values))
        ds = ds.reindex(simulation=simulations)

    if prettify:
        ds = prettify_names(ds)

    return ds
#----
