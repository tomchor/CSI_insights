import numpy as np
import pynanigans as pn
import xarray as xr
from aux00_utils import open_simulation
from matplotlib import pyplot as plt
π = np.pi

#++++ Define directory and simulation name
path = f"simulations/data/"
snames = ["PNN_CIfront1",
          #"PNN_CIfront2",
          #"PNN_CIfront3",
          #"PNN_SIfront1",
          #"PNN_SIfront2",
          #"PNN_SIfront3",
          "PNN_SIfront4",
          #"PNN_SIfront5",
          #"PNN_SIfront6",
          #"PNN_CIintjet01",
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
    grid_avg, avg = open_simulation(path+f"avg.{sname}.nc", 
                                    use_inertial_periods=True,
                                    topology=sname[:3],
                                    open_dataset_kwargs=dict(chunks=dict(time=1)),
                                    )
    #----

    #+++++++ BPE calculation
    avg["BPE"] = -(avg.b_sorted * avg.zC)
    avg["dBPEdt"] = avg.BPE.load().differentiate("time")/avg.T_inertial
    avg["ε_p"] = avg.dBPEdt
    #-------

 
    #++++ Get important times and perform average
    times_avg = slice(None, None, 5)
    
    avg = avg.sel(time=times_avg)
    avg_0d = avg.mean(("zC", "zF"))
    #-----

    #+++++ Get SGS buoyancy flux
    qb = (- vid.κ_e * grid_vid.interp(vid.dbdz, 'z')).load()
    qb_abs = abs(qb)
    wb_abs = abs(vid.wb_res)
    avg_0d["qb"] = qb.pnmean(('y', 'z'))
    avg_0d["qb_abs"] = qb_abs.pnmean(('y', 'z'))
    avg_0d["wb_abs"] = wb_abs.pnmean(('y', 'z'))
    #-----
    
    #++++ Integrate dissipations in time
    εk_int = avg_0d.ε.integrate("time")
    εp_int = avg_0d.ε_p.integrate("time")
    εs_int = avg_0d.sponge_dissip.integrate("time")
    wb_int = avg_0d.wb_res.integrate("time")
    qb_int = avg_0d.qb.integrate("time")
    #----

    #++++ Ratio between dissipations
    print(f"Ratio between ε_sponge and ε_k for simulation {sname} is ", float(εs_int / εk_int))
    print(f"Ratio between w'b' and ε_k for simulation {sname} is ", float(wb_int / εk_int))
    print(f"Ratio between w'b' and ε_p for simulation {sname} is ", float(wb_int / εp_int))
    print(f"Ratio between w'b' and q_b for simulation {sname} is ", float(wb_int / qb_int))
    #----

    #++++ Plot results
    plt.figure()
    #avg_0d.ε.plot(label=r"Mean KE dissipation rate $\varepsilon_k$")
    avg_0d.ε_p.plot(label=r"Mean buoyancy mixing rate $\varepsilon_p$")
    (-avg_0d.wb_res).plot(label=r"$-\langle w'b'\rangle$")

    fig, ax = plt.gcf(), plt.gca()
    ax.set_title(sname)
    ax.legend(); ax.grid(True)
    fig.savefig(f"figures_check/{sname}_wb_εp.png")
    #----
