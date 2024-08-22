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
xlab =  { "CK" : "Random Value ( Real / Imaginary ))"
        , "IK" : "Random Value ( Integer )"
        , "RK" : "Random Value ( Real )"
        }
legends =   [ r"$\alpha = +.1, x_{min} = +3., x_{max} = +10$"
            , r"$\alpha = +2., x_{min} = +0., x_{max} = +5$"
            , r"$\alpha = +.5, x_{min} = +0., x_{max} = +8$"
            ]

for kind in ["IK", "CK", "RK"]:

    pattern = "*." + kind + ".txt"
    fileList = glob.glob(pattern)
    if len(fileList) == 1:

        df = pd.read_csv(fileList[0], delimiter = ", ")

        fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
        ax = plt.subplot()

        if kind == "CK":
            ax.hist ( df.values[:, 0]
                    , bins = 30
                    , histtype = "stepfilled"
                    , density = True
                    , alpha = 0.7
                    )
            ax.hist ( df.values[:,1]
                    , bins = 30
                    , histtype = "stepfilled"
                    , density = True
                    , alpha = 0.7
                    )
        else:
            ax.hist ( df.values[:,:]
                    , bins = 50
                    , histtype = "stepfilled"
                    , density = True
                    , alpha = 0.7
                    )
            ax.legend   ( legends[::-1]
                        , fontsize = fontsize
                        )

        plt.xticks(fontsize = fontsize - 2)
        plt.yticks(fontsize = fontsize - 2)
        ax.set_xlabel(xlab[kind], fontsize = 17)
        ax.set_ylabel("Density", fontsize = 17)
        ax.set_title("Histogram of {} randomly generated values".format(len(df.values)), fontsize = fontsize)
        #ax.set_yscale("log")
        #ax.set_xscale("log")

        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")
        ax.set_axisbelow(True)

        plt.tight_layout()
        plt.savefig(fileList[0].replace(".txt",".png"))

    elif len(fileList) > 1:

        sys.exit("Ambiguous file list exists.")