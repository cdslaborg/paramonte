#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

linewidth = 2
fontsize = 17

marker ={ "CKG" : "-"
        , "IKG" : "."
        , "RKG" : "-"
        }
xlab =  { "CKG" : "Generalized Gamma Random Value ( real/imaginary components )"
        , "IKG" : "Generalized Gamma Random Value ( integer-valued )"
        , "RKG" : "Generalized Gamma Random Value ( real-valued )"
        }
legends =   [ "$\kappa, 1/\omega, 1/\sigma = 1.0, 0.5, 0.5$"
            , "$\kappa, 1/\omega, 1/\sigma = 2.0, 0.5, 1.0$"
            , "$\kappa, 1/\omega, 1/\sigma = 0.5, 2.0, 0.5$"
            , "$\kappa, 1/\omega, 1/\sigma = 0.2, 5.0, 0.2$"
            , "$\kappa, 1/\omega, 1/\sigma = .14, 7.0, .14$"
            , "$\kappa, 1/\omega, 1/\sigma = 2.0, 5.0, 0.3$"
            ]

for kind in ["IKG", "CKG", "RKG"]:

    pattern = "*." + kind + ".txt"
    fileList = glob.glob(pattern)
    if len(fileList) == 1:

        df = pd.read_csv(fileList[0], delimiter = ",", header = None)

        fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
        ax = plt.subplot()

        for j in range(len(df.values[0,:])):
            if kind == "CKG":
                plt.hist( df.values[:,j]
                        , histtype = "stepfilled"
                        , alpha = 0.5
                        , bins = 75
                        )
            else:
                plt.hist( df.values[:,j]
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
        ax.set_title("Histograms of {} Generalized Gamma random values".format(len(df.values[:, 0])), fontsize = 17)
        #ax.set_yscale("log")

        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")

        plt.savefig(fileList[0].replace(".txt",".png"))

    elif len(fileList) > 1:

        sys.exit("Ambiguous file list exists.")