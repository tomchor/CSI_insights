import numpy as np
import xarray as xr
from aux01_grid import get_grid
from cmocean import cm
import holoviews as hv
hv.extension('bokeh')
import hvplot.xarray
π = np.pi

# Preamble
#++++++

# Define dir and file names
#++++
dirname = "inst_2dbounded"
dirname = "CIinst_2D"
path = f"../{dirname}/"
fname = "instFNN_SIsurfjet4"
#----

# Open dataset
#++++
vid = xr.open_dataset(path+f"vid.{fname}.nc", chunks={"time": 50})
vid = vid.squeeze(drop=True)
if 'yC' in vid.coords:
    horcoord = 'yC'
else:
    horcoord = 'xC'
#------
#------


#+++++ Change time units
tf = [0, 0.5, 1, 1.5, 2, 5, 10]
tf = [0, 2, 5]

vid = vid.assign_coords(time = vid.time/vid.T_inertial); vid.time.attrs = dict(untis="inertial times")
vid = vid.sel(time=tf, method="nearest")
#------


#+++++ Define GRID
if "FPN" in fname:
    grid = get_grid(vid, periodic=['x', 'y'])
    boundary = None
elif "FNN" in fname:
    grid = get_grid(vid, periodic=['x'], topology="FNN")
    boundary = "extend"
#------


#++++ Calculations
dudy = vid.u_tot.differentiate("yC")
dbdz = grid.interp(vid.dbdz, 'z')
barot = (vid.f_0 - dudy) * dbdz

dudz = vid.u_tot.differentiate("zC")
baroc = -vid.f_0*dudz**2
#----


#+++++ Plotting options
geom = dict(width=550, height=400, 
            xlim=(vid.yF.min().item(), vid.yF.max().item()), 
            ylim=(vid.zF.min().item(), vid.zF.max().item()),
            shared_axes=False)
unit = dict(clim=(-1.5, 1.5), cmap=cm.balance)
dissip = dict(clim=(-3e-7, 3e-7), cmap=cm.balance)

ncols = 3
#-----

#++++ PLOT
Ro_plot = vid.Ro.hvplot(title="Rossby number", x='yF', y='zC', **(geom | unit))
Ri_plot = vid.Ri.hvplot(title="Richardson number", x='yC', y='zF', **(geom | unit))

SPy_plot = vid.SP_y.hvplot(title="Shear Prod y", x='yC', y='zC', **(geom | dissip))
SPz_plot = vid.SP_z.hvplot(title="Shear Prod z", x='yC', y='zC', **(geom | dissip))

PV_plot = vid.PV.hvplot(title="Potential Vorticity", x='yF', y='zF', 
                        clim=(-5e-8, 5e-8), cmap=cm.balance, **geom)
ωx_plot = vid.ω_x.hvplot(title="ω_x", x='yF', y='zF',
                         clim=(-1e-2, 1e-2), cmap=cm.balance, **geom)

BT_plot = barot.hvplot(title=r"$(f+\zeta) N^2$", x='yC', y='zC',
                       **geom)
BC_plot = baroc.hvplot(title=r"$f (\partial u / \partial z)^2$", x='yC', y='zC',
                       **geom)

PV_zero = vid.PV.hvplot.contour(x='yF', y='zF', levels=[0], line_width=0.5, **geom)
b_cont = vid.b_tot.hvplot.contour(x='yC', y='zC', levels=30, line_width=0.5, **geom)
#----

#+++++ ASSEMBLE AND SAVE
plot = (Ro_plot +\
        Ri_plot +\
        PV_plot * PV_zero +
        SPy_plot +\
        SPz_plot +\
        ωx_plot * b_cont +\
        BT_plot +\
        BC_plot 
        ).cols(ncols)

hv.save(plot, f"figures_check/{fname}_snapshots.html")
#-----
