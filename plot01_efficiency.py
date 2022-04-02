import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
π = np.pi

#++++ Define directory and simulation name
dirname = "ISI_jet"
path = f"/glade/u/home/tomasc/scratch_cheyenne/{dirname}/data/"
snames = ["FNN_CIsurfjet1",
          #"FNN_SIsurfjet4",
          "PNN_CIintjet01",
          #"PNN_CIsurfjet1",
          #"PNN_SIsurfjet2",
          #"PNN_SIsurfjet3",
          #"PNN_SIsurfjet4",
          #"PNN_SIsurfjet5",
          #"PNN_SIsurfjet6",
          ]
#----


for sname in snames:
    #++++ Open datasets avg and vid
    print(f"Opening {sname}")
    grid_vid, vid = open_simulation(path+f"vid.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    squeeze=True,
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    ds_eff = xr.load_dataset(f"data/efficiencies_{sname}.nc")
    #----

    #+++++ Plot efficiency metrics
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12, 5),
                             constrained_layout=True,
                             sharex=True,)
    
    ds_eff.dBPEdt.plot(ax=axes[0,0], label=r"$dBPE/dt $")
    ds_eff.ε.plot(ax=axes[0,0], label="Dissipation")
    ds_eff.γ.plot(ax=axes[0,1], label="Bounded efficiency (no sponge dissipation)")
    ds_eff.γ2.plot(ax=axes[0,1], label="Bounded efficiency (with sponge dissipation)")
    
    ds_eff.nomin_int.plot(ax=axes[1,0], label=r"$\int dBPE/dt  dt$")
    ds_eff.denom_int.plot(ax=axes[1,0], label=r"$\int( dBPE/dz  + \epsilon)dt$")
    ds_eff.denom_int2.plot(ax=axes[1,0], label=r"$\int( dBPE/dz  + \epsilon + \epsilon_s)dt$")
    ds_eff.Γ.plot(ax=axes[1,1], label="Integrated Bounded efficiency (no sponge dissipation)")
    ds_eff.Γ2.plot(ax=axes[1,1], label="Integrated Bounded efficiency (with sponge dissipation)")
    #ds_eff.γ_weighted.plot(ax=axes[1,1], label="Weighted (by ε) integral of Bounded efficiency")
    
    for ax in axes.flatten():
        ax.grid()
        ax.legend()
    for ax in axes[:,1]:
        ax.set_ylim(-0.1, .6)
    fig.savefig(f"figures_check/efficiency_{sname}.pdf")
    #-----
        

