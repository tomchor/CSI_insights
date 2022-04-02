import numpy as np
from itertools import chain


#++++ Instability angles
def plot_angles(angles, ax):
    vlim = ax.get_ylim()
    hlim = ax.get_xlim()
    length = np.max([np.diff(vlim), np.diff(hlim)])
    for α in np.array(angles):
        v0 = np.mean(vlim)
        h0 = np.mean(hlim)
        h1 = h0 + length*np.cos(α)
        v1 = v0 + length*np.sin(α)
        ax.plot([h0, h1], [v0, v1], c='k', ls='--')
    return
#----


#++++ Define colors and markers
color_base = ["b", "C1", "C2", "C3", "C5", "C8"]
marker_base = ["o", "v", "P"]

colors = color_base*len(marker_base)
markers = list(chain(*[ [m]*len(color_base) for m in marker_base ]))
#markers = marker_base*len(color_base)
#colors = list(chain(*[ [m]*len(marker_base) for m in color_base ]))
#----


#++++ Standardized plotting
def plot_scatter(ds, ax=None, x=None, y=None, hue="simulation", add_guide=True, cycle_kwargs=dict(), **kwargs):
    for i, sim in enumerate(ds[hue].values):
        #++++ Getting values for specific point
        xval = ds.sel(**{hue:sim})[x]
        yval = ds.sel(**{hue:sim})[y]
        marker = markers[i]
        color = colors[i]
        cycle_kw_i = { k:v[i] for k, v in cycle_kwargs.items() }
        #----

        #++++ Define label (or not)
        if add_guide:
            label=sim
        else:
            label=""
        #----

        #++++ Plot
        ax.scatter(xval, yval, c=color, marker=marker, label=sim, **(kwargs | cycle_kw_i))
        #----

        #++++ Include labels
        try:
            ax.set_xlabel(xval.attrs["long_name"])
        except:
            ax.set_xlabel(xval.name)

        try:
            ax.set_ylabel(yval.attrs["long_name"])
        except:
            ax.set_ylabel(yval.name)
        #----

    return 
#----


#++++ Letterize plot axes
def letterize(axes, x, y, coords=True, bbox=dict(boxstyle='round', 
                                                        facecolor='white', alpha=0.8),
                     **kwargs):
    from string import ascii_lowercase
    for ax, c in zip(axes.flatten(), ascii_lowercase*2):
        ax.text(x, y, c, transform=ax.transAxes, bbox=bbox, **kwargs)
    return
#----

