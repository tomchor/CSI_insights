import numpy as np
import xarray as xr
π = np.pi


#++++ Initial condition
def get_IC_front(y, z, y_0=4e3, z_0=0, σ_y=1600, σ_z=80, u_0=-.2, N2_inf=1e-6, f_0=1e-4, Lz=80):
    """ Gets initial condition for a front set-up """
    from scipy.special import erf

    f_y = np.exp(-(y-y_0)**2/σ_y**2)
    F_y = 1/2 * np.sqrt(π) * σ_y * (erf((y-y_0)/σ_y) + 1)
    f_z = (z-z_0)/σ_z + 1
    df_z = 1/σ_z

    u = u_0 * f_y * f_z
    b = -u_0 * f_0 * F_y * df_z + N2_inf*(z+Lz)

    Ro = -u.differentiate(y.name) / f_0
    Ri = N2_inf / u.differentiate(z.name)**2
    q_hat = 1 + Ro - 1/Ri

    dsout = xr.Dataset(dict(u=u, b=b,
                            Ro=Ro, Ri=Ri,
                            q_hat=q_hat,))
    return dsout


def get_IC_jet(y, z, y_0=4e3, z_0=-250, σ_y=1600, σ_z=80, u_0=-.2, N2_inf=1e-6, f_0=1e-4, Lz=80): 
    """ Gets initial condition for an interior jet set-up """
    from scipy.special import erf

    f_y = np.exp(-(y-y_0)**2/σ_y**2)
    F_y = 1/2 * np.sqrt(π) * σ_y * (erf((y-y_0)/σ_y) + 1)
    f_z = np.exp(-(z-z_0)**2/σ_z**2)
    df_z = -2 * (z-z_0)/σ_z**2 * f_z

    u = u_0 * f_y * f_z
    b = -u_0 * f_0 * F_y * df_z + N2_inf * (z + Lz)

    Ro = -u.differentiate(y.name) / f_0
    Ri = b.differentiate(z.name) / u.differentiate(z.name)**2
    q_hat = 1 + Ro - 1/Ri

    dsout = xr.Dataset(dict(u=u, b=b,
                            Ro=Ro, Ri=Ri,
                            q_hat=q_hat,))
    return dsout
#----


#+++++ Convenience functions
def adjust_variables(ds):
    """
    Adjusts 
        σy to σ_y
        σz to σ_z
        f0 to f_0
        b0 to b_0
    """
    if "σy" in ds.attrs:
        ds.attrs["σ_y"] = ds.σy
    if "σz" in ds.attrs:
        ds.attrs["σ_z"] = ds.σz
    if "f0" in ds.attrs:
        ds.attrs["f_0"] = ds.f0
    if "b0" in ds.attrs:
        ds.attrs["b_0"] = ds.b0
    return ds


def get_relevant_attrs(ds):
    """
    Get what are deemed relevant attrs from ds.attrs and return them
    as a dataset
    """
    wantedkeys = ['N2_inf', 'σ_y', 'z_0', 'name', 'interval', 'Lz', 
                  'T_inertial', 'σ_z', 'u_0', 'Nz', 'Ny', 'date', 'Ly', 'b_0', 
                  'Oceananigans', 'y_0', 'νz', 'f_0', 'LES']
    attrs = dict((k, v) for k,v in ds.attrs.items() if k in wantedkeys)
    return attrs
#-----


#+++++ Functions for refence quantities
def get_ref_quantities(ds, setup="front"):
    """ 
    Returns a DataArray with reference quantities from dataset `ds` 
    including everything already in its attributes
    """
    params = dict()

    y_qmin, z_qmin, Ro_qmin, Ri_qmin = get_qmin(ds, setup=setup)
    params["y_qmin"] = y_qmin
    params["z_qmin"] = z_qmin
    params["Ro_qmin"] = Ro_qmin
    params["Ri_qmin"] = Ri_qmin
    return params


def get_qmin(ds, setup="front"):

    vardict = dict(y=ds.yF, z=ds.zF, 
                   y_0=ds.y_0, z_0=ds.z_0, 
                   σ_y=ds.σ_y, σ_z=ds.σ_z, 
                   u_0=ds.u_0, N2_inf=ds.N2_inf, 
                   f_0=ds.f_0, Lz=ds.Lz)

    if setup=="front":
        IC = get_IC_front(**vardict)
    elif setup=="jet":
        IC = get_IC_jet(**vardict)
    else:
        raise NameError
    q_hat = IC.q_hat

    q_hat_min = q_hat.where(q_hat==q_hat.min(), drop=True).squeeze()
    y_qmin = float(q_hat_min.yF)
    if setup=="front":
        z_qmin = float(q_hat_min.zF)
    elif setup=="jet":
        try:
            z_qmin = float(q_hat_min.zF[-1])
        except IndexError:
            z_qmin = float(q_hat_min.zF)
    else:
        raise NameError
    Ro_qmin = float(IC.Ro.sel(yF=y_qmin, zF=z_qmin))
    Ri_qmin = float(IC.Ri.sel(yF=y_qmin, zF=z_qmin))
    return y_qmin, z_qmin, Ro_qmin, Ri_qmin
#-----


#+++++ Functions for bulk quantities
def get_bulk_parameters(ds, Pr=1, ν_m=1.05e-6):
    """ Gets bulk parameters based on flow geometry """
    params = dict()

    params["Ro_bulk"] = ds.u_0 / (ds.σ_y * ds.f_0)
    params["Ri_bulk"] = ds.σ_z**2 * ds.N2_inf / ds.u_0**2
    if ds.LES:
        params["Re_bulk"] = ds.u_0 * ds.σ_y / ν_m
    else:
        params["Re_bulk"] = ds.u_0 * ds.σ_y / ds.νz
    params["PV_bulk"] = ds.N2_inf * ds.f_0 * (1 + params["Ro_bulk"] - 1/params["Ri_bulk"])
    params["Ratio"] = ds.σ_y / ds.σ_z
    params["delta_ar"] = ds.σ_z / ds.σ_y
    if "Pr" not in params:
        params["Pr"] = Pr
    else:
        params["Pr"] = ds.Pr

    return params
#-----


#++++ Functions for Kloosterziel's adjustment hypothesis
def zero_crossings(arr):
    import numpy as np
    return np.where(np.diff(np.sign(arr)))[0]

def find_right_bound(m_0, r_l, dim='r', verbose=False):
    #m_c = m_0.sel(**{dim : r_l, "method" : "nearest"})
    m_c = m_0.interp({dim : r_l})
    crossings = m_0[dim].isel(**{dim : zero_crossings(m_0 - m_c)})
    if verbose:
        print(f"rₗ= {r_l}, m_c = {m_c.item()}")
        print(f"Crossings: {crossings.values}")

    if len(crossings)==0:
        r_r = m_0[dim][-1]
        print(f"No crossings for rₗ= {float(r_l)}. Assuming rᵣ= {r_r.item()}")
    else:
        r_r = crossings[-1]
    return r_r

def get_Δm(m_0, r_l, dim='r', verbose=False):
    """
    Calculate ∫ₗʳ(m_c - m(y)) dy, which then needs to be minimized
    according to equation (1.1) of doi:10.1017/S0022112007006593
    """
    r_l = float(r_l)
    r_r = find_right_bound(m_0, r_l=r_l, dim=dim, verbose=verbose)
    m0_chunk = m_0.sel(**{dim : slice(r_l, r_r)})
    if verbose: print(f"Length of chunk is {len(m0_chunk)}")
    if len(m0_chunk)==0:
        Δm_int = 99
    else:
        mc = m_0.interp({dim : r_l})
        Δm = m0_chunk - mc
        Δm_int = Δm.integrate(dim)
    if verbose: print(f"rₗ= {r_l},    rᵣ = {r_r.item()},     Δm_int = {float(Δm_int)}\n")
    return Δm_int


def predict_final_state(m_0, dim='r',
                        left_bound=None, right_bound=None, 
                        verbose=False, **kwargs):
    """
    Predicts the final state of a DataArray m_0 for total momentum.

    Example:
    m_final = predict_final_state(m_dim, verbose=True, dim='y', method="COBYLA")
    """

    if left_bound is None:
        left_bound = m_0.where(m_0<m_0.sel(**{dim : 0, "method" : "nearest"}), drop=True)[dim][0]

    if right_bound is None:
        right_bound = m_0[dim].isel({dim : zero_crossings(m_0.differentiate(dim))[0]})

    initial_guess = (left_bound + right_bound)/2
    if verbose: 
        print(f"Left and right bounds for minimization are {float(left_bound)} and {float(right_bound)}")
        print(f"Initial guess: {initial_guess.values}\n")

    func = lambda x: get_Δm(m_0, x, verbose=verbose)**2
    result = minimize(func, initial_guess, bounds=[(left_bound, right_bound)], **kwargs)

    r_l = float(result.x)
    if verbose: print("Finding final r_l")
    r_r = find_right_bound(m_0, r_l=r_l, dim=dim, verbose=False)

    m_c = m_0.interp({dim : r_l})
    m_f = m_0.where((m_0[dim]<r_l) | (m_0[dim]>r_r), m_c)
    return m_f
#----


#++++ Calculation of background (unavailable) potential energy
def order_xyz(darray):
    """ 
    Orders values of `darray` by flattening every dimension except for time.
    Useful for calculating available potential energy
    """
    import numpy as np
    import dask as da

    stacked_dims = [ dim for dim in darray.dims if dim != "time" ]
    da_aux = darray.stack(aux_dim=stacked_dims)
    if "time" in darray.dims:
        da_aux = da_aux.transpose("time", "aux_dim")
    if darray.chunks:
        sort = lambda x: da.array.map_blocks(np.sort, x, axis=-1)
    else:
        sort = lambda x: np.sort(x, axis=-1)
    da_ord = xr.apply_ufunc(sort, da_aux, dask="allowed")
    return da_ord.unstack()


def zipsort(da_a, da_b, unsorted_dim="time",
            return_indices=False, indices=None):
    """
    Sort two DataArrays based on the first one only.

    Only works if both `da_a` and `da_b` are chunked in `unsorted_dim`
    with size 1 chunks
    """
    from dask.array import map_blocks
    assert da_a.dims == da_b.dims, "Dimensions aren't the same"
    for dim in da_a.dims:
        assert np.allclose(da_a[dim], da_b[dim]), f"Coordinates of {dim} aren't the same"

    sorted_dims = [ dim for dim in da_a.dims if dim != unsorted_dim ]
    daa_aux = da_a.stack(aux_dim=sorted_dims).transpose(unsorted_dim, "aux_dim") # stack all dims to be sorted into one
    dab_aux = da_b.stack(aux_dim=sorted_dims).transpose(unsorted_dim, "aux_dim") # stack all dims to be sorted into one

    if indices is None:
        indices = map_blocks(np.argsort, daa_aux.data, axis=-1, dtype=np.int64)
    if return_indices:
        return indices

    def reorder(A, ind): return A[0,ind]
    daa_aux.data = map_blocks(reorder, daa_aux.data, indices, dtype=np.float64)
    dab_aux.data = map_blocks(reorder, dab_aux.data, indices, dtype=np.float64)
    return daa_aux.unstack(), dab_aux.unstack()
#----


#+++++ Instability angles
def get_growth_rate(α, M2=0, F2=0, N2=0):
    T1 = N2 * np.sin(α)**2
    T2 = F2 * np.cos(α)**2
    T3 = M2 * np.sin(2*α)
    return -T1 -T2 +T3

def get_angles(M2=0, F2=0, N2=0):
    α0 = 1/2 * np.arctan2(2*M2, (N2-F2))
    α1 = α0 + π/2
    return [α0, α1]
#-----


#++++ Other stuff
def get_Ri(ds, grid, zderiv_kw=dict(boundary="extrapolate")):
    if "dbdz" not in ds.variables.keys():
        ds["dbdz"] = grid.derivative(ds.b, "z", **zderiv_kw)

    if "dvdz" not in ds.variables.keys():
        ds["dvdz"] = grid.interp(grid.derivative(ds.v, "z", **zderiv_kw), "y")

    if "dudz" not in ds.variables.keys():
        ds["dudz"] = grid.derivative(ds.u_tot, "z", **zderiv_kw)

    Ri = ds.dbdz /(ds.dudz**2 + ds.dvdz**2)
    return Ri
#----

