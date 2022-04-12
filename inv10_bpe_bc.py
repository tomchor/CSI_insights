import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from dask.diagnostics import ProgressBar
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIfront1",
#          "PNN_CIfront2",
#          "PNN_CIfront3",
#          "PNN_SIfront2",
#          "PNN_SIfront3",
#          "PNN_SIfront4",
#          "PNN_SIfront5",
#          "PNN_SIfront6",
#          "PNN_CIintjet01",
          ]
#----



allrefs = xr.open_dataset("data/allrefs.nc")
for i, sname in enumerate(snames):
    print(f"Opening {sname}")

    #++++ Open datasets avg and vid
    grid_out, out = open_simulation(path+f"out.{sname}.nc", 
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
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    )
    #----


    #----
    out = out.sel(time=slice(0, 0.9))
    avg = avg.sel(time=slice(0, 0.9))
    avg["BPE"] = -avg.zC * avg.b_sorted
    avg["dBPEdt"] = avg.BPE.pnmean('z').differentiate('time')
    out = out.isel(time=1)
    out.load()
    #----

    #----
    out.attrs = out.attrs | dict(Lx=out.ΔxC.sum())
    A = out.Lx * out.Ly

    def get_H_func_int(b, b_prime, **kwargs):
        """ b needs to be a float
            b_prime needs to be a DataArray but can't depend on time
        """
        H_func = np.heaviside(b_prime-b, 1/2)
        if type(H_func) is xr.DataArray:
            return H_func.integrate(('xC', 'yC', 'zC'))
        else:
            return np.trapz(np.trapz(np.trapz(H_func, axis=0), axis=1), axis=2)

    def get_z_star(b, b_prime, **kwargs):
        """ b here can be an dataArray but can't depend on time
            b_prime needs to be a DataArray but can't depend on time
        """
        A = out.Lx * out.Ly
        z_star = xr.apply_ufunc(get_H_func_int, b, b_prime, **kwargs) / A
        return z_star
    #----

    #----
    # We're using Dask ararys so we need to use map_blocks
    b = out.b_tot.isel(yC=-1, zC=-1)#, time=1)
    b_prime = out.b_tot#.isel(time=1)
    pause
    z_star = xr.map_blocks(get_z_star, b, args=[b_prime])
    #----
