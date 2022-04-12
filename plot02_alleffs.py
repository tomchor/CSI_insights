import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["PNN_CIfront1",
          "PNN_CIfront2",
          "PNN_CIfront3",
          "PNN_SIfront2",
          "PNN_SIfront3",
          "PNN_SIfront4",
          "PNN_SIfront5",
          "PNN_SIfront6",
          ]
#----


#+++++ Plot efficiency metrics
from matplotlib import pyplot as plt
fig, axes = plt.subplots(ncols=1, nrows=3, figsize=(9, 9),
                         constrained_layout=True,
                         sharex=False, squeeze=False)
#-----
 
for sname in snames:
    #++++ Open datasets avg
    print(f"Opening {sname}")
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    load=False,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    ds_eff = xr.load_dataset(f"data/efficiencies_{sname}.nc")
    shearprod = xr.load_dataset(f"data/shearprod_{sname}.nc")
    #----

    #++++ Offset calculations
    Γ = ds_eff.Γ.isel(time=-1)
    #----

    #++++ Auxiliary variables
    φ_c = np.degrees(np.arctan(-1-avg.Ro_r))
    φ_Ri = np.degrees(np.arctan(-1/avg.Ri_r))

    Ro_b = avg.u_0/(avg.σy * avg.f_0)
    Ri_b = avg.σz**2 * avg.N2_inf / avg.u_0**2
    F2_b = avg.f_0**2 * (1 - Ro_b)

    RiRo = avg.Ro_r * avg.Ri_r
    F2_r = avg.f_0**2 * (1 - avg.Ro_r)
    RiRo_complete = -RiRo*(1 - F2_r/avg.N2_inf)
    SPy_ratio = (shearprod.SPy_lin/(shearprod.SPy_lin + shearprod.SPz_lin)).values
    #-----

    axes[0,0].bar(φ_Ri, Γ, width=5, alpha=0.8,
                  label=f"ϕ = {φ_Ri:.2f}, {sname}")
    axes[1,0].bar(RiRo_complete, Γ, width=0.2, alpha=0.8,
                  label=f"$-Ri Ro (1-F^2/N^2)$ = {RiRo_complete:.2f}, {sname}")
    axes[2,0].bar(SPy_ratio, Γ, width=0.1, alpha=0.8, 
                  label=f"{sname}")

axes[0,0].set_xlabel(r"$\phi_{Ri_b}$")
axes[1,0].set_xlabel(r"$-Ri Ro (1 - F^2/N^2)$ (Horizontal shear / Vertical shear)")
axes[2,0].set_xlabel(r"Measured horizontal shear / Measured total shear")

for ax in axes.flatten():
    ax.legend()
    ax.grid()
    ax.set_ylabel(f"$\Gamma$ (mixing efficiency)")

fig.savefig(f"figures_check/efficiencies_barplot.pdf")
