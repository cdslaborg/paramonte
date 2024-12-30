#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 17

marker ={ "IK" : "."
        , "RK" : "-"
        , "CK" : "-"
        }
ylab =  { "IK" : "Discrete Integer Uniform CDF"
        , "RK" : "Continuous Real Uniform CDF"
        , "CK" : "Continuous Complex Uniform CDF"
        }

for kind in ["IK", "RK", "CK"]:

    df = pd.read_csv("main.unif.cdf."+kind+".txt", delimiter = " ")

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    if kind == "CK":
        plt.plot( df.values[:, 0]
                , df.values[:,2]
                , marker[kind]
                , color = "r"
                )
        plt.plot( df.values[:, 1]
                , df.values[:,3]
                , marker[kind]
                , color = "blue"
                )
        ax.legend   ( ["real", "imaginary"]
                    , fontsize = fontsize
                    )
    else:
        plt.plot( df.values[:, 0]
                , df.values[:, 1]
                , marker[kind]
                , color = "r"
                )

    ax.set_xlabel("X", fontsize = 17)
    ax.set_ylabel(ylab[kind], fontsize = 17)

    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    plt.savefig("getUnifCDF."+kind+".png")
