#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

linewidth = 2
fontsize = 17

marker ={ "CK" : "-"
        , "IK" : "."
        , "RK" : "-"
        }
xlab =  { "CK" : "redshift: z ( real/imaginary components )"
        , "IK" : "redshift: z ( integer-valued )"
        , "RK" : "redshift: z ( real-valued )"
        }
legends =   [ r"$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) = (0.3, 0.7, 0.0, 0.0)$"
            , r"$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) = (1.0, 0.0, 0.0, 0.0)$"
            , r"$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) = (0.0, 1.0, 0.0, 0.0)$"
            , r"$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) = (.30, .30, .40, 0.0)$"
            , r"$(\Omega_M, \Omega_\Lambda, \Omega_R, \Omega_K) = (.25, .25, .25, .25)$"
            ]

for kind in ["IK", "CK", "RK"]:

    pattern = "*." + kind + ".txt"
    fileList = glob.glob(pattern)
    if len(fileList) == 1:

        df = pd.read_csv(fileList[0], delimiter = ", ")

        fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
        ax = plt.subplot()

        if kind == "CK":
            plt.plot( df.values[:, 0] - 1
                    , df.values[:,1:6]
                    , marker[kind]
                    , linewidth = linewidth
                   #, color = "r"
                    )
            plt.plot( df.values[:, 1] - 1
                    , df.values[:,1:6]
                    , marker[kind]
                    , linewidth = linewidth
                   #, color = "blue"
                    )
        else:
            plt.plot( df.values[:, 0] - 1
                    , df.values[:,1:6]
                    , marker[kind]
                    , linewidth = linewidth
                   #, color = "r"
                    )
            ax.legend   ( legends
                        , fontsize = fontsize
                        )

        plt.xticks(fontsize = fontsize - 2)
        plt.yticks(fontsize = fontsize - 2)
        ax.set_xlabel(xlab[kind], fontsize = 17)
        ax.set_ylabel("The Hubble Parameter [km / s / Mpc]", fontsize = 17)
        #ax.set_xlim([0.0, 5])
        #ax.set_ylim([0.05, 0.5])
        #ax.set_xscale("log")
        ax.set_yscale("log")

        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")

        plt.tight_layout()
        plt.savefig(fileList[0].replace(".txt",".png"))

    elif len(fileList) > 1:

        sys.exit("Ambiguous file list exists.")