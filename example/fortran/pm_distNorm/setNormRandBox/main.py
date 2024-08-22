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
xlab =  { "CK" : "Normal Random Number ( real/imaginary components )"
        , "IK" : "Normal Random Number ( integer-valued )"
        , "RK" : "Normal Random Number ( real-valued )"
        }
legends =   [ r"$\mu = -5.,~\sigma = 1.0$"
            , r"$\mu = 0.0,~\sigma = 1.0$"
            , r"$\mu = 2.0,~\sigma = 3.0$"
            ]

for kind in ["IK", "CK", "RK"]:

    pattern = "*." + kind + ".txt"
    fileList = glob.glob(pattern)
    if len(fileList) == 1:

        df = pd.read_csv(fileList[0], delimiter = " ", header = None)

        fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
        ax = plt.subplot()

        if kind == "CK":
            plt.hist( df.values[:,0:3]
                    , histtype = "stepfilled"
                    , alpha = 0.5
                    , bins = 75
                    )
        else:
            plt.hist( df.values[:,0:3]
                    , histtype = "stepfilled"
                    , alpha = 0.5
                    , bins = 75
                    )
            ax.legend   ( legends
                        , fontsize = fontsize
                        )
        plt.xticks(fontsize = fontsize - 2)
        plt.yticks(fontsize = fontsize - 2)
        ax.set_xlabel(xlab[kind], fontsize = 17)
        ax.set_ylabel("Count", fontsize = 17)
        ax.set_title("Histograms of {} Normal random numbers".format(len(df.values[:, 0])), fontsize = 17)

        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")

        plt.savefig(fileList[0].replace(".txt",".png"))

    elif len(fileList) > 1:

        sys.exit("Ambiguous file list exists.")