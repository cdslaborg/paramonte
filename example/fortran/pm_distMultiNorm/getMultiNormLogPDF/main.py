#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import glob
import sys

linewidth = 2
fontsize = 17

marker ={ "CK" : "-"
        , "IK" : "."
        , "RK" : "-"
        }
xlab =  { "CK" : "X ( real/imaginary components )"
        , "IK" : "X ( integer-valued )"
        , "RK" : "X ( real-valued )"
        }
legends =   [ r"$\mu = 0.0,~\sigma = 3.0$"
            , r"$\mu = 0.0,~\sigma = 1.0$"
            , r"$\mu = 0.0,~\sigma = 0.3$"
            , r"$\mu = -2.,~\sigma = 1.0$"
            ]

for kind in ["IK", "CK", "RK"]:

        # Make 1d plot.

    pattern = "*.D1."+kind+".txt"
    fileList = glob.glob(pattern)
    if len(fileList) == 1:

        df = pd.read_csv(fileList[0], delimiter = ",")

        fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 300)
        ax = plt.subplot()

        if kind == "CK":
            plt.plot( df.values[:, 0]
                    , df.values[:,1:5]
                    , marker[kind]
                    , linewidth = linewidth
                   #, color = "r"
                    )
            plt.plot( df.values[:, 1]
                    , df.values[:,1:5]
                    , marker[kind]
                    , linewidth = linewidth
                   #, color = "blue"
                    )
        else:
            plt.plot( df.values[:, 0]
                    , df.values[:,1:5]
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
            ax.set_ylabel("Probability Density Function (PDF)", fontsize = 17)

            plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
            ax.tick_params(axis = "y", which = "minor")
            ax.tick_params(axis = "x", which = "minor")

            plt.tight_layout()
            plt.savefig(fileList[0].replace(".txt",".png"))

    elif len(fileList) > 1:

        sys.exit("Ambiguous file list exists.")
    # Make 2d plot.

    pattern = "*.D2."+kind+".txt"
    fileList = glob.glob(pattern)
    if len(fileList) == 1:

        df = pd.read_csv(fileList[0], delimiter = ",", header = None)
        df = df.values

        npnt = math.isqrt(len(df[:, 0]))
        gridx = np.reshape(df[:, 0], newshape = (npnt, npnt), order = 'F')
        gridy = np.reshape(df[:,1], newshape = (npnt, npnt), order = 'F')
        gridz = np.reshape(df[:,2], newshape = (npnt, npnt), order = 'C')

        fig, ax = plt.subplots(subplot_kw = {"projection": "3d"})
        fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 300)
        ax = fig.add_subplot(1, 1, 1, projection = '3d')
        
        ax.plot_surface(gridx, gridy, gridz, cmap = 'GnBu', linewidth = 0)
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel("PDF", fontsize = 17)
        
        plt.xticks(fontsize = fontsize - 2)
        plt.yticks(fontsize = fontsize - 2)
        
        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")
        
        plt.tight_layout()
        plt.savefig(fileList[0].replace(".txt",".png"))

    elif len(fileList) > 1:

        sys.exit("Ambiguous file list exists.")