import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from aux01_physfuncs import adjust_variables
from dask.diagnostics import ProgressBar
π = np.pi

#++++ Computation options
import dask
if __name__ == '__main__':
    compute_flags = dict(num_workers=18, memory_limit='5GB') # Processes doesn't work very well
else:
    compute_flags = dict()

extra = "_mask"
#----

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
          "PNN_CIfront1_f2",
          "PNN_CIfront1_f4",
          "PNN_CIfront1_f8",
          "PNN_SIfront4_f2",
          "PNN_SIfront4_f4",
          "PNN_SIfront4_f8",
          "PNN_CIfront1_AMD",
          "PNN_CIfront1_AMD_f2",
          "PNN_CIfront1_AMD_f4",
          "PNN_CIfront1_AMD_f8",
          "PNN_SIfront4_AMD",
          "PNN_SIfront4_AMD_f2",
          "PNN_SIfront4_AMD_f4",
          "PNN_SIfront4_AMD_f8",
          "PNN_CIfront1_AMD2",
          "PNN_CIfront1_AMD2_f2",
          "PNN_CIfront1_AMD2_f4",
          "PNN_CIfront1_AMD2_f8",
          "PNN_SIfront4_AMD2",
          "PNN_SIfront4_AMD2_f2",
          "PNN_SIfront4_AMD2_f4",
          "PNN_SIfront4_AMD2_f8",
          "PNN_CIfront1_AMD3",
          "PNN_CIfront1_AMD3_f2",
          "PNN_CIfront1_AMD3_f4",
          "PNN_CIfront1_AMD3_f8",
          "PNN_SIfront4_AMD3",
          "PNN_SIfront4_AMD3_f2",
          "PNN_SIfront4_AMD3_f4",
          "PNN_SIfront4_AMD3_f8",
          ]
#----


ν_m = 1e-6
allparams = xr.load_dataset("data_post/allparams.nc")
for i, sname in enumerate(snames):
    print(f"Opening {sname}")
    print(f"Extra: {extra}")

    #++++ Open datasets avg and vid
    print("Calculating Reb")
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
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #----

    #++++ Get subdomain
    print("Getting a σ_y-wide slice of the domain")
    out = adjust_variables(out)
    if "front" in sname:
        y_r = allparams.sel(simulation=sname).y_qmin
    else:
        y_r = 1/2 * np.sqrt(2) * out.σ_y + out.y_0

    if "mask" in extra:
        side = 2*out.σ_y
    else:
        side = 1/2*out.σ_y
    yF_r = out.sel(yF=y_r+side, method="nearest").yF
    yF_l = out.sel(yF=y_r-side, method="nearest").yF
    yslice = slice(yF_l, yF_r) # Make sure len(yF) = len(yC)+1

    if "halfz" in extra:
        depth= 40
    else:
        depth=np.inf
    zF_r = out.sel(zF=out.z_0+depth, method="nearest").zF
    zF_l = out.sel(zF=out.z_0-depth, method="nearest").zF
    zslice = slice(zF_l, zF_r) # Make sure len(zF) = len(zC)+1

    tslice = slice(1, None)

    out = out.sel(yC=yslice, yF=yslice, zC=zslice, zF=zslice, time=tslice)
    subgrid_out = pn.get_grid(out, topology="PNN")
    #-----

    #++++ Get conditional avg mask
    if "cond" in extra:
        print("Getting masks for conditional averaging")
        ε_c = 1e-10
        mask_ccc = out.ε > ε_c
        mask_ccf = subgrid_out.interp(out.ε, 'z', boundary="extend") > ε_c
        mask_cfc = subgrid_out.interp(out.ε, 'y', boundary="extend") > ε_c
        mask_fcc = subgrid_out.interp(out.ε, 'x', boundary="extend") > ε_c

        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            out["ε"] = out.ε.where(mask_ccc, drop=True)
            out["ν_e"] = out.ν_e.where(mask_ccc, drop=True)
            out["dbdz"] = out.dbdz.where(mask_ccf, drop=True)
            out["u"] = out.u.where(mask_fcc, drop=True)
            out["v"] = out.v.where(mask_cfc, drop=True)
            out["w"] = out.w.where(mask_ccf, drop=True)

    if "mask" in extra:
        print("Getting masks for IC averaging")
        allICs = xr.open_dataset("data_post/allICs.nc").chunk(dict(simulation=1))
        mask_acc = allICs.sel(simulation=sname).mask_acc
        mask_acf = allICs.sel(simulation=sname).mask_acf
        mask_afc = allICs.sel(simulation=sname).mask_afc
    #----

    #++++ Calculations
    print("Calculating...")
    def calc_variables(ds, grid=None):
        if grid is None: 
            grid = pn.get_grid(ds, topology="PNN")

        if "mask" in extra:
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                ε_mask = ds.ε.where(mask_acc, drop=True)
                νe_mask = ds.ν_e.where(mask_acc, drop=True)
                dbdz_mask = ds.dbdz.where(mask_acf, drop=True)
                ds["ε_mean"] = ε_mask.pnmean(('x', 'y', 'z'))
                ds["νe_mean"] = νe_mask.pnmean(('x', 'y', 'z'))
                ds["dbdz_mean"] = dbdz_mask.pnmean(('x', 'y', 'z'))
                ds["Re_b_strain"] = (ds.ε / ds.ν_e).where(mask_acc, drop=True).pnmean(('x', 'y', 'z')) / ds.N2_inf
        else:
            ds["ε_mean"] = ds.ε.pnmean(('x', 'y', 'z'))
            ds["νe_mean"] = ds.ν_e.pnmean(('x', 'y', 'z'))
            ds["dbdz_mean"] = ds.dbdz.pnmean(('x', 'y', 'z'))
            ds["Re_b_strain"] = (ds.ε / ds.ν_e).pnmean(('x', 'y', 'z')) / ds.N2_inf

        ds["Re_b_avg_sgs"] = ds.ε_mean / (ds.νe_mean * ds.N2_inf)
        ds["Re_b_avg_molec"] = ds.ε_mean / (ν_m * ds.N2_inf)
        ds["Re_b_point_sgs"] = ds.ε_mean / (ds.νe_mean * ds.dbdz_mean)
        ds["Re_b_point_molec"] = ds.ε_mean / (ν_m * ds.dbdz_mean)

        ds["dvdx"] = grid.derivative(ds.v, "x")
        ds["dudy"] = grid.derivative(ds.u, "y", boundary="extend")
        ds["Ro"] = (ds.dvdx - ds.dudy)/ds.f_0
        if "mask" in extra:
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                Ro_mask = ds.Ro.where(mask_afc, drop=True)
                ds["Ro_mean"] = Ro_mask.pnmean(('x', 'y', 'z'))
        else:
            ds["Ro_mean"] = ds.Ro.pnmean(('x', 'y', 'z'))

        ds["dudz2"] = grid.interp(grid.derivative(ds.u, "z", boundary="extend")**2, "x")
        ds["dvdz2"] = grid.interp(grid.derivative(ds.v, "z", boundary="extend")**2, "y")
        ds["S2"] = (ds.dudz2 + ds.dvdz2)
        ds["Ri"] = (ds.dbdz / ds.S2)
        if "mask" in extra:
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                Ri_mask = ds.Ri.where(mask_acf, drop=True)
                ds["Ri_mean"] = Ri_mask.pnmean(('x', 'y', 'z'))
        else:
            ds["Ri_mean"] = ds.Ri.pnmean(('x', 'y', 'z'))

        ds["u2"] = grid.interp(ds.u**2, "x")
        ds["v2"] = grid.interp(ds.v**2, "y")
        ds["w2"] = grid.interp(ds.w**2, "z")
        ds["ke"] = ds.u2 + ds.v2 + ds.w2
        if "mask" in extra:
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                ke_mask = ds.ke.where(mask_acc, drop=True)
                ds["ke_mean"] = ke_mask.pnmean(("x", "y", "z"))
        else:
            ds["ke_mean"] = ds.ke.pnmean(("x", "y", "z"))

        ds["Fr_avg_mean"] = ds.ε_mean / (ds.ke_mean * np.sqrt(ds.N2_inf))
        ds["Fr_point_mean"] = ds.ε_mean / (ds.ke_mean * np.sqrt(ds.dbdz_mean))

        ds["u_p"] = ds.u - ds.u.mean('xF')
        ds["v_p"] = ds.v - ds.v.mean('xC')
        ds["w_p"] = ds.w - ds.w.mean('xC')
        ds["u_p2"] = grid.interp(ds.u_p**2, "x")
        ds["v_p2"] = grid.interp(ds.v_p**2, "y")
        #ds["w_p2"] = grid.interp(ds.w_p**2, "z")
        #ds["ke_p"] = ds.u_p2 + ds.v2 + ds.w2
        #ds["ke_p"] = ds.u_p2 + ds.v_p2 + ds.w_p2
        ds["ke_p"] = ds.u_p2 + ds.v_p2
        if "mask" in extra:
            with dask.config.set(**{'array.slicing.split_large_chunks': True}):
                ke_p_mask = ds.ke_p.where(mask_acc, drop=True)
                ds["ke_p_mean"] = ke_p_mask.pnmean(("x", "y", "z"))
        else:
            ds["ke_p_mean"] = ds.ke_p.pnmean(("x", "y", "z"))

        ds["Fr_p_avg_mean"] = ds.ε_mean / (ds.ke_p_mean * np.sqrt(ds.N2_inf))
        ds["Fr_p_point_mean"] = ds.ε_mean / (ds.ke_p_mean * np.sqrt(ds.dbdz_mean))

        return ds

    out = calc_variables(out, subgrid_out)
    #----

    #++++ Get ε peak
    ε_avg = avg.ε.sel(time=slice(None, 5, 10)).pnmean('z')
    t_εmax = ε_avg.time[ε_avg.argmax('time').values]
    #----


    #++++ Save results
    print("Starting to save")
    out.Re_b_avg_sgs.attrs = dict(long_name=r"$Re_b = mean(\epsilon) / mean(\nu_e) N^2_\infty$")
    out.Re_b_avg_molec.attrs = dict(long_name=r"$Re_b = mean(\epsilon) / \nu_m N^2_\infty$")

    out.Re_b_point_sgs.attrs = dict(long_name=r"$Re_b = mean(\epsilon) / mean(\nu_e) mean(db/dz)$")
    out.Re_b_point_molec.attrs = dict(long_name=r"$Re_b = mean(\epsilon) / \nu_m mean(db/dz)$")

    out.Fr_avg_mean.attrs = dict(long_name=r"$Fr = mean(\epsilon) / mean(KE) N_\infty$")
    out.Fr_point_mean.attrs = dict(long_name=r"$Fr = mean(\epsilon) / mean(KE) \sqrt{mean(db/dz)}$")

    out.Fr_p_avg_mean.attrs = dict(long_name=r"$Fr = mean(\epsilon) / mean(KE') N_\infty$")
    out.Fr_p_point_mean.attrs = dict(long_name=r"$Fr = mean(\epsilon) / mean(KE') \sqrt{mean(db/dz)}$")

    dsout = xr.Dataset(dict(Re_b_avg_sgs=out.Re_b_avg_sgs, 
                            Re_b_avg_molec=out.Re_b_avg_molec,
                            Re_b_point_sgs=out.Re_b_point_sgs,
                            Re_b_point_molec=out.Re_b_point_molec,
                            Ri_mean=out.Ri_mean,
                            Ro_mean=out.Ro_mean,
                            ε_mean=out.ε_mean,
                            Fr_avg_mean=out.Fr_avg_mean,
                            Fr_point_mean=out.Fr_point_mean,
                            Fr_p_avg_mean=out.Fr_p_avg_mean,
                            Fr_p_point_mean=out.Fr_p_point_mean,
                            νe_mean=out.νe_mean,
                            t_εmax=t_εmax,
                            ))

    fname = f"data_post/Reb_{sname}{extra}.nc"
    print(f"Saving results to {fname}")
    with ProgressBar():
        dsout.to_netcdf(fname)
    #----
