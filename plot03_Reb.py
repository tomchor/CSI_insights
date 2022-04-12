import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, pnames
from aux01_physfuncs import adjust_variables
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
π = np.pi


#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIfront1",
          "PNN_CIfront2",
          "PNN_CIfront3",
          "PNN_CIfront4",
          "PNN_CIfront5",
          "PNN_SIfront1",
          "PNN_SIfront2",
          "PNN_SIfront3",
          #"PNN_SIfront4",
          "PNN_SIfront5",
          "PNN_SIfront6",
#          "PNN_CIintjet01",
          ]
extra = "_mask"
#----

#+++++ Start figure and markers
nrows=3; ncols=2
size = 4
fig_all, axes_all = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.3*ncols*size, nrows*size),
                                 constrained_layout=True, squeeze=False)
axesf_all = axes_all.flatten()

markers = ["o", "^", "s", "p", "D", "P", "x", "*", ">", "<", "8"]
#-----

#++++ Plotting functions
def scatter_points(axes, add_guide=True, **kwargs):
    vmin0=-5; vmax0=0; hue0="RoRi_r"
    vmin1=1; vmax1=10; hue1="time"
    vmin2=-5; vmax2=0; hue2="RoRi_outer"
    vmin3=0; vmax3=0.25; hue3="γ_coeff"
    vmin4=0; vmax4=20; hue4="Re_b_avg_sgs"
    vmin5=0; vmax5=50; hue5="Re_b_point_molec"
    vmin6=0; vmax6=1e-2; hue6="Fr_p_avg_mean"
    vmin7=0; vmax7=1e-2; hue7="Fr_p_point_mean"
    vmin8=0; vmax8=0.25; hue8="γ"
    vmin9=0; vmax9=50; hue9="Re_b_avg_molec"
    if "cond" in extra:
        vmin10=5e-11; vmax10=5e-9; hue10="ε_mean"
    else:
        vmin10=1e-11; vmax10=1e-9; hue10="ε_mean"
    vmin11=0; vmax11=5; hue11="-RoRi_r"
    vmin12=0; vmax12=20; hue12="Re_b_point_sgs"

    dsplot.plot.scatter(ax=axes[0], x="Re_b_avg_sgs", y="γ", hue=hue10,
                        vmin=vmin10, vmax=vmax10, norm=LogNorm(),
                        add_guide=add_guide, marker=markers[i], label=pnames[sname],
                        **kwargs)

    dsplot.plot.scatter(ax=axes[1], x="Re_b_avg_molec", y="γ", hue=hue10,
                        vmin=vmin10, vmax=vmax10, norm=LogNorm(),
                        add_guide=add_guide, marker=markers[i], label=pnames[sname],
                        **kwargs)

    dsplot.plot.scatter(ax=axes[2], x="Re_b_avg_sgs", y="-RoRi_outer", hue=hue11,
                        vmin=vmin11, vmax=vmax11,
                        add_guide=add_guide, marker=markers[i], label=pnames[sname],
                        **kwargs)

    dsplot.plot.scatter(ax=axes[3], x="-RoRi_outer", y="γ", hue=hue11,
                        vmin=vmin11, vmax=vmax11,
                        add_guide=add_guide, marker=markers[i], label=pnames[sname],
                        **kwargs)

    dsplot.plot.scatter(ax=axes[5], x="Ri_mean", y="γ", hue=hue11,
                        vmin=vmin11, vmax=vmax11,
                        add_guide=add_guide, marker=markers[i], label=pnames[sname],
                        **kwargs)

    return None


def prettify_ax(ax, slope=True, xlog=True, ylog=True, 
                ysymlog=False, xsymlog=False, legend=False):
    if slope:
        add_slope(ax, slope=slope)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    if ysymlog:
        ax.set_yscale("symlog", linthresh=0.01)
    if xsymlog:
        ax.set_xscale("symlog", linthresh=0.01)
    if legend: ax.legend()
    ax.grid(True)
    return

def add_slope(ax, coeff=1e-1, slope=-1/2, log_xlim=(-1, 0)):
    xcoor = np.logspace(*log_xlim)
    ax.plot(xcoor, coeff*xcoor**(slope), c="k")
    return
#----


allparams = xr.open_dataset("data/allparams.nc")
for i, sname in enumerate(snames):
    print(f"Opening {sname}")

    #+++++ Load datasets
    if extra:
        ds_Reb = xr.load_dataset(f"data/Reb_{sname}{extra}.nc")
    else:
        ds_Reb = xr.load_dataset(f"data/Reb_{sname}.nc")
    dseff = xr.load_dataset(f"data/efficiencies_{sname}.nc").squeeze()

    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    avg = adjust_variables(avg)
    #-----


    #++++ Process avg
    ε_avg = avg.ε.sel(time=slice(None, 5, 10)).pnmean('z')
    #----

    #+++++ Make plotting easier
    dseff = dseff.sel(time=ds_Reb.time)
    dsplot = xr.Dataset(dict(Re_b_avg_sgs=ds_Reb.Re_b_avg_sgs, 
                             Re_b_avg_molec=ds_Reb.Re_b_avg_molec,
                             Re_b_point_sgs=ds_Reb.Re_b_point_sgs,
                             Re_b_point_molec=ds_Reb.Re_b_point_molec,
                             Ro_mean=ds_Reb.Ro_mean, Ri_mean=ds_Reb.Ri_mean, ε_mean=ds_Reb.ε_mean,
                             Fr_p_avg_mean=ds_Reb.Fr_p_avg_mean, Fr_p_point_mean=ds_Reb.Fr_p_point_mean,
                             νe_mean=ds_Reb.νe_mean, 
                             γ=dseff.γ, Γ=dseff.Γ,
                             γ_coeff=dseff.γ_coeff, Γ_coeff=dseff.Γ_coeff,
                             ))
    dsplot["Ri_r"] = allparams.sel(simulation=sname).Ri_r + 0*dsplot.time
    dsplot["RoRi_r"] = (allparams.Ri_r * allparams.Ro_r).sel(simulation=sname) + 0*dsplot.time
    dsplot["-RoRi_r"] = -(allparams.Ri_r * allparams.Ro_r).sel(simulation=sname) + 0*dsplot.time
    dsplot["RoRi_outer"] = dsplot.Ro_mean * dsplot.Ri_mean
    dsplot["-RoRi_outer"] = -dsplot.Ro_mean * dsplot.Ri_mean

    t_εmax = ε_avg.time[ε_avg.argmax('time').values]
    dsplot = dsplot.sel(time=slice(t_εmax,None)) # Remove laminar flow
    dsplot = dsplot.where(dsplot.ε_mean>1e-10) # Remove low-turbulence flows
    #-----

    #++++ Prepare names for plotting
    dsplot.RoRi_outer.attrs = dict(long_name=r"$\langle{Ro}\rangle_s\,\langle{Ri}\rangle_s$")
    dsplot["-RoRi_outer"].attrs = dict(long_name=r"$-\langle{Ro}\rangle_s\,\langle{Ri}\rangle_s$")
    dsplot["-RoRi_r"].attrs = dict(long_name=r"$-Ro_r \, Ri_r$")
    dsplot.Re_b_avg_sgs.attrs = dict(long_name=r"$Re_b^\mathrm{sgs}$")
    dsplot.Re_b_avg_molec.attrs = dict(long_name=r"$Re_b^\mathrm{mol}$")
    dsplot.ε_mean.attrs = dict(long_name=r"$\langle{\epsilon}\rangle_s$", units=r"m$^2$/s$^3$")
    #----

    #+++++ Plot all simulations
    add_guide = True if i==0 else False
    scatter_points(axesf_all, add_guide=add_guide, 
                   cbar_kwargs=dict(shrink=0.7, aspect=50))
    #-----

#++++ Make plot pretty
for i, ax in enumerate(axes_all[0,:]):
    prettify_ax(ax, slope=False)
add_slope(axes_all[0,0], log_xlim=(0, 1.2), slope=-1/2, coeff=1.2e-1)
add_slope(axes_all[0,1], log_xlim=(1.5, 3), slope=-1/2, coeff=7e-1)
for i, ax in enumerate(axes_all[1,:]):
    prettify_ax(ax, slope=False, xsymlog=True, ylog=True)

for ax in axesf_all:
    ax.set_title("")
    ax.set_ylim(1e-2, None)

#if "cond" in extra:
#    for i, ax in enumerate(axes_all[:,2]):
#        ax.set_xlim(None, -1e-2)
axesf_all[2].legend()
#----

#++++ Save plot
fig_all.savefig(f"figures_check/eff_Reb_all{extra}_debug.pdf")
#----
