#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

fontsize = 17

marker ={ "CK" : "-"
        , "IK" : "."
        , "RK" : "-"
        }
xlab =  { "CK" : r"shape parameter: $\kappa$ ( real/imaginary )"
        , "IK" : r"shape parameter: $\kappa$ ( integer-valued )"
        , "RK" : r"shape parameter: $\kappa$ ( real-valued )"
        }

for kind in ["IK", "CK", "RK"]:

    pattern = "*." + kind + ".txt"
    fileList = glob.glob(pattern)
    if len(fileList) == 1:

        df = pd.read_csv(fileList[0], delimiter = " ")

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
        else:
            plt.plot( df.values[:, 0]
                    , df.values[:, 1]
                    , marker[kind]
                    , color = "r"
                    )

        plt.xticks(fontsize = fontsize - 2)
        plt.yticks(fontsize = fontsize - 2)
        ax.set_xlabel(xlab[kind], fontsize = fontsize)
        ax.set_ylabel("Normalization Factor of the PDF", fontsize = fontsize)
        ax.set_xlim([0, 8.])
        ax.set_ylim([0, 1.2])

        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")

        plt.savefig(fileList[0].replace(".txt",".png"))

    elif len(fileList) > 1:

        sys.exit("Ambiguous file list exists.")