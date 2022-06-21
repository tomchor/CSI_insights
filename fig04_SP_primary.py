import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation, filter_sims
from matplotlib import pyplot as plt
from cmocean import cm
from aux02_plotting import plot_scatter, letterize



#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIfront1",
          "PNN_SIfront4",
          ]
alleffs = xr.load_dataset(f"data_post/alleffs.nc", decode_times=False)
allprods = xr.load_dataset("data_post/allshearprod.nc", decode_times=False)
allparams = xr.load_dataset("data_post/allparams.nc", decode_times=False)

alleffs = filter_sims(alleffs, prettify=True, only_main=True)
allprods = filter_sims(allprods, prettify=True, only_main=True)
allparams = filter_sims(allparams, prettify=True, only_main=True)
#----

#++++ Prepare for plotting
dsall = xr.merge([alleffs, allprods])
dsall["Ri_qmin"] = allparams.Ri_qmin
dsall["N2_inf"] = allparams.N2_inf
dsall["M2_0"] = - allparams.f_0 * allparams.u_0 / allparams.σ_z
dsall["N2_norm"] = allparams.N2_inf / allparams.f_0**2
dsall["N2_norm2"] = allparams.N2_inf / dsall.M2_0

dsall.Ri_qmin.attrs = dict(long_name=r"$Ri_r$")
dsall.N2_inf.attrs = dict(long_name=r"$N^2_\infty$")
dsall.M2_0.attrs = dict(long_name=r"$M^2_0$")
dsall.N2_norm.attrs = dict(long_name=r"$N^2_\infty / f^2$")
dsall.N2_norm2.attrs = dict(long_name=r"$N^2_\infty / M_0^2$")

dsall["ratio_measured"] = abs(dsall.SPy_lin / dsall.SPz_lin)
dsall.ratio_measured.attrs = dict(long_name=r"$R_\mathrm{SP}^\mathrm{prim} = \langle SP_h^\mathrm{prim} \rangle/\langle SP_v^\mathrm{prim}\rangle$")

dsall.Γ_last.attrs = dict(long_name=r"$\Gamma_\infty$")

F2_r = allparams.f_0**2 * (1 - allparams.Ro_qmin)
ratio_approx = allparams.Ri_qmin * allparams.Ro_qmin
ratio_complete = -ratio_approx*(1 - F2_r/allparams.N2_inf)

ratio_complete.attrs = dict(long_name=r"$R_\mathrm{SP}^\mathrm{prim} = -Ri_r Ro_r\left[1 - \frac{f^2}{N^2}(1+Ro)\right]$")
dsall["ratio_complete"] = ratio_complete

ratio_approx.attrs = dict(long_name=r"$R_\mathrm{SP}^\mathrm{prim} \approx -Ri_r Ro_r$")
dsall["ratio_approx"] = -ratio_approx

γ_max = dsall.γ.max("time")
γ_max.attrs = dict(long_name=r"$\mathrm{max}(\gamma)$")
dsall["γ_max"] = γ_max
#----

#++++ Figure parames (ncols, nrows, vars)
ncols=2
nrows=1
size=4.5

fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(ncols*size, nrows*size/1.2),
                         squeeze=False, constrained_layout=True,
                         sharex=False, sharey=True)
axesf = axes.flatten()
#-----

#+++ Plot!
for ax in axesf:
    ax.set_xscale("log")
    ax.axhline(y=0.17, c="gray", ls="--", zorder=0)
plot_scatter(dsall, ax=axesf[0], x="ratio_measured", y="Γ_last", hue="simulation", add_guide=True, s=80, edgecolors='k')
plot_scatter(dsall, ax=axesf[1], x="ratio_approx", y="Γ_last", hue="simulation", add_guide=False, s=80, edgecolors='k')
#plot_scatter(dsall, ax=axesf[2], x="N2_norm", y="Γ_last", hue="simulation", add_guide=False, s=80, edgecolors='k')
#plot_scatter(dsall, ax=axesf[3], x="N2_norm2", y="Γ_last", hue="simulation", add_guide=False, s=80, edgecolors='k')

handles, labels = axesf[0].get_legend_handles_labels()

plot_scatter(dsall, ax=axesf[0], x="ratio_measured", y="γ_max", hue="simulation", add_guide=False, s=30, alpha=0.5)
plot_scatter(dsall, ax=axesf[1], x="ratio_approx", y="γ_max", hue="simulation", add_guide=False, s=30, alpha=0.5)
#plot_scatter(dsall, ax=axesf[2], x="N2_norm", y="γ_max", hue="simulation", add_guide=False, s=30, alpha=0.5)
#plot_scatter(dsall, ax=axesf[3], x="N2_norm2", y="γ_max", hue="simulation", add_guide=False, s=30, alpha=0.5)
#----

#++++ Make it pretty
for ax in axesf:
    ax.set_ylim(0, None)
    ax.grid(True)
    ax.set_ylabel("Cumulative mixing efficiency ($\Gamma_\infty$)")
    ax.set_title("")
axesf[1].legend(handles, labels, loc="upper left", bbox_to_anchor=(1., 1), fontsize='x-small', labelspacing=1.1)

axesf[0].set_xlabel(dsall.ratio_measured.attrs["long_name"])
axesf[1].set_xlabel(dsall.ratio_approx.attrs["long_name"])
#----

#-----
letterize(axesf, 0.1, 0.9, fontsize=14)
fig.suptitle("Cumulative mixing efficiency as a function of the ratio\nbetween horizontal and vertical shear production rates")
fig.savefig(f"figures_paper/effs_vs_SP.pdf")
#-----

