import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from dask.diagnostics import ProgressBar
π = np.pi

#++++ Define directory and simulation name
path = f"simulations/data/"
snames = [#"FNN_CIintjet01",
          #"FNN_CIsurfjet1",
          "FNN_CIsurfjet3",
          #"FNN_SIsurfjet4",
          ]
#----


for sname in snames:
    #++++ Open datasets avg and vid
    print(f"Opening {sname}")
    grid_out, out = open_simulation(path+f"out.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #----

 
    #++++ Get important times and perform average
    if sname == "FNN_CIintjet01":
        t_max1 = 2.9
        t_max2 = 3.1
        window = dict(y=slice(6200, 7000), z=slice(-300, -200))

    elif sname == "FNN_CIsurfjet1":
        t_max1 = 1.70
        t_max2 = 2.00
        window = dict(y=slice(5200, 5400), z=slice(-40, 0))

    elif sname == "FNN_CIsurfjet3":
        t_max1 = 1.40
        t_max2 = 1.65
        window = dict(y=slice(4000, 4500), z=slice(-20, 0))

    elif sname == "FNN_SIsurfjet4":
        t_max1 = 0.90
        t_max2 = 1.25
        window = dict(y=slice(3700, 3900), z=slice(None, -20))

    else:
        raise NameError(f"Define vaues for simulation {sname}")

    newindow = dict()
    ykeys = ["yF", "yC"]
    zkeys = ["zF", "zC"]
    for k, v in window.items():
        if k=="y":
            newindow.update({ykey : v for ykey in ykeys})
        if k=="z":
            newindow.update({zkey : v for zkey in zkeys})
    window = newindow
    #-----
    
    #++++ Selecting only our window
    vidsnaps = vid.sel(time=slice(t_max1, t_max2))
    vidwind = vidsnaps.sel(**window)

    gridwind = pn.get_grid(vidwind, topology="FNN")
    #----


    #++++ Define calculations as in-place functions
    def add_Ri(ds, grid=None):
        if grid is None: 
            grid = pn.get_grid(ds, topology="FNN")

        ds["dbdz_int"] = grid.interp(ds.dbdz, 'z')
        ds["v_int"] = grid.interp(ds.v, 'y')
        ds["Ri"] = ds.dbdz_int / (ds.u.differentiate('zC')**2 + ds.v_int.differentiate('zC')**2)

        return ds

    def add_background(ds, t_init):
        ds["U_bg"] = ds.sel(time=t_init, method="bfill").u_tot
        ds["V_bg"] = ds.sel(time=t_init, method="bfill").v
        ds["W_bg"] = ds.sel(time=t_init, method="bfill").w
        ds["B_bg"] = ds.sel(time=t_init, method="bfill").b_tot
        return ds

    def add_fluctuations(ds, grid=None):
        if grid is None: 
            grid = pn.get_grid(ds, topology="FNN")

        bd = dict(boundary="extend")
        ds["u_fcc"] = ds.u_tot - ds.U_bg
        ds["v_cfc"] = ds.v - ds.V_bg
        ds["w_ccf"] = ds.w - ds.W_bg
        ds["b_acc"] = ds.b_tot - ds.B_bg

        ds["u_ccc"] = grid.interp(ds.u_fcc, 'x', **bd)
        ds["u_cfc"] = grid.interp(ds.u_ccc, 'y', boundary="extend")
        ds["u_ccf"] = grid.interp(ds.u_ccc, 'z', boundary="extend")

        ds["v_ccf"] = grid.interp(grid.interp(ds.v_cfc, 'y', **bd), 'z', **bd)

        ds["w_cfc"] = grid.interp(grid.interp(ds.w_ccf, 'y', **bd), 'z', **bd)

        ds["dUdy"] = grid.derivative(ds.U_bg, 'y', boundary="extend")
        ds["dVdy"] = grid.derivative(ds.V_bg, 'y', boundary="extend")
        ds["dWdy"] = grid.derivative(ds.W_bg, 'y', boundary="extend")

        ds["dUdy_cfc"] = grid.interp(ds.dUdy, 'x', **bd)
        ds["dVdy_cfc"] = grid.interp(ds.dVdy, 'y', **bd)
        ds["dWdy_cfc"] = grid.interp(ds.dWdy, 'z', **bd)

        ds["dUdz"] = grid.derivative(ds.U_bg, 'z', **bd)
        ds["dVdz"] = grid.derivative(ds.V_bg, 'z', **bd)
        ds["dWdz"] = grid.derivative(ds.W_bg, 'z', **bd)

        ds["dUdz_ccf"] = grid.interp(ds.dUdz, 'x', **bd)
        ds["dVdz_ccf"] = grid.interp(ds.dVdz, 'y', **bd)
        ds["dWdz_ccf"] = grid.interp(ds.dWdz, 'z', **bd)

        ds["b_acf"] = grid.interp(ds.b_acc, 'z', **bd)

        return ds

    def add_shear_prod(ds):
        ds["SPy_2nd"] = -ds.u_cfc * ds.v_cfc * ds.dUdy_cfc \
                    -ds.v_cfc**2 * ds.dVdy_cfc \
                    -ds.w_cfc * ds.v_cfc * ds.dWdy_cfc

        ds["SPz_2nd"] = -ds.u_ccf * ds.w_ccf * ds.dUdz_ccf \
                    -ds.v_ccf * ds.w_ccf * ds.dVdz_ccf \
                    -ds.w_ccf**2 * ds.dWdz_ccf

        ds["wb_2nd"] = ds.w_ccf * ds.b_acf
        return ds

    def add_windowmean(ds):
        ds["wb_2nd_mean"] = ds.wb_2nd.pnmean(('x', 'y', 'z'))
        ds["SPy_2nd_mean"] = ds.SPy_2nd.pnmean(('x', 'y', 'z'))
        ds["SPz_2nd_mean"] = ds.SPz_2nd.pnmean(('x', 'y', 'z'))
        return ds
    #----

    #++++ Apply calculations
    vidwind = add_Ri(vidwind, grid=gridwind)
    vidwind = add_background(vidwind, t_max1)
    vidwind = add_fluctuations(vidwind, grid=gridwind)
    vidwind = add_shear_prod(vidwind)
    vidwind = add_windowmean(vidwind)
    #-----

    #+++++ Save to disk
    with ProgressBar():
        vidwind.to_netcdf(f"data_post/shearprod_secondary_{sname}.nc")
    #-----
