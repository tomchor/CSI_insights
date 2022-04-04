import numpy as np
import pynanigans as pn
import xarray as xr
import pandas as pd

#++++ Set rosetta stone and variables of interest
variables_main = ["σ_y", "f_0", "u_0", "N2_inf", "Ro_qmin", "Ri_qmin", "delta_ar", "Γ_last"]
variables_aux = ["σ_y", "f_0", "u_0", "N2_inf", "νz", "Ro_qmin", "Ri_qmin", "delta_ar"]

rosetta = dict(σ_y = r"$\sigma_y$ (m)",
               σ_z = r"$\sigma_z$ (m)",
               f_0 = r"$f$ (1/s)",
               u_0 = r"$u_0$ (m/s)",
               N2_inf = r"$N^2_0$ (1/s$^2$)",
               νz = r"$\nu_e$ (m$^2$/s)",
               Ro_qmin = r"$Ro_r$",
               Ri_qmin = r"$Ri_r$",
               delta_ar = r"$\delta$",
               simulation = "Simulation",
               Γ_last = r"$\Gamma_\infty$",
               )
simstone = dict(PNN_ = "",
                FNN_ = r"2D\_",
                surfjet = r"front",
                )
#----


#++++ Choose which simulations are main and which are aux
main_sims = ["PNN_CIsurfjet1",
             "PNN_CIsurfjet2",
             "PNN_CIsurfjet3",
             "PNN_CIsurfjet4",
             "PNN_CIsurfjet5",
             "PNN_SIsurfjet1",
             "PNN_SIsurfjet2",
             "PNN_SIsurfjet3",
             "PNN_SIsurfjet4",
             "PNN_SIsurfjet5",
             "PNN_SIsurfjet6",
             ]

aux_sims = ["FNN_CIsurfjet1",
            "FNN_CIsurfjet3",
            "FNN_SIsurfjet4",
            #"PNN_CIintjet01",
            ]
#----

#++++ Open dataset
allparams = xr.load_dataset("data_post/allparams.nc")
alleffs = xr.load_dataset("data_post/alleffs.nc")

allparams["νz"] = allparams.νz.where(np.logical_not(allparams.LES), other=0)
allparams = xr.merge([allparams, alleffs.Γ_last])
#----

#++++ Separate into two
params_main = allparams[variables_main]
params_main = params_main.reindex(simulation=main_sims)
params_main = params_main.rename({k:v for k, v in rosetta.items() if k in params_main.variables.keys()})

params_aux = allparams[variables_aux]
params_aux = params_aux.reindex(simulation=aux_sims)
params_aux = params_aux.rename({k:v for k, v in rosetta.items() if k in params_aux.variables.keys()})
#----

#++++ Replace inside indices and turn params into 
def replace_indices(df):
    df = df.to_dataframe()
    df.index = pd.Series(df.index).replace(simstone, regex=True)
    return df

params_main = replace_indices(params_main)
params_aux = replace_indices(params_aux)
#----

#++++
# Change format of columns that require less precision
def change_format(df):
    if rosetta["νz"] in df.columns:
        df[rosetta["νz"]] = df[rosetta["νz"]].map(lambda x: '\SI{%.1e}{}' % x)

    if rosetta["Γ_last"] in df.columns:
        df[rosetta["Γ_last"]] = df[rosetta["Γ_last"]].map(lambda x: '%3.2f' % x)

    df[rosetta["delta_ar"]] = df[rosetta["delta_ar"]].map(lambda x: '\SI{%.2e}{}' % x)
    df[rosetta["f_0"]] = df[rosetta["f_0"]].map(lambda x: '\SI{%.1e}{}' % x)
    df[rosetta["u_0"]] = df[rosetta["u_0"]].map(lambda x: '\SI{%.1e}{}' % x)
    df[rosetta["N2_inf"]] = df[rosetta["N2_inf"]].map(lambda x: '\SI{%.1e}{}' % x)
    df[rosetta["Ro_qmin"]] = df[rosetta["Ro_qmin"]].map(lambda x: '$%2.1f$' % x)
    df[rosetta["Ri_qmin"]] = df[rosetta["Ri_qmin"]].map(lambda x: '%2.1f' % x)
    return df

params_main = change_format(params_main)
params_aux = change_format(params_aux)
#-----

#++++ Output!
options = dict(escape=False, multicolumn_format="l")
print(params_main.to_latex(**options))
print(params_aux.to_latex(**options))

params_main.to_latex("figures_paper/params_main.tex", encoding='utf-8', **options)
params_aux.to_latex("figures_paper/params_aux.tex", encoding='utf-8', **options)
#----
