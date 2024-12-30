#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys


fontsize = 17

prefixes =  { "setKmeansPP" : "Memberships based on k-means++."
            , "setKmeans"   : "k-means refinement of k-means++ memberships."
            }
for prefix in list(prefixes.keys()):

    legends = []
    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    fileList = glob.glob(prefix + ".*.txt")
    if len(fileList) == 2:

        for file in fileList:

            kind = file.split(".")[1]
            df = pd.read_csv(file, delimiter = ",", header = None)

            if kind == "center":
                ax.scatter  ( df.values[:, 1]
                            , df.values[:,2]
                            , zorder = 100
                            , marker = "*"
                            , c = "red"
                            , s = 50
                            )
                legends.append("center")
            elif kind == "sample":
                ax.scatter  ( df.values[:, 1]
                            , df.values[:,2]
                            , c = df.values[:, 0]
                            , s = 10
                            )
                legends.append("sample")
            else:
                sys.exit("Ambiguous file exists: {}".format(file))

        ax.legend(legends, fontsize = fontsize)
        plt.xticks(fontsize = fontsize - 2)
        plt.yticks(fontsize = fontsize - 2)
        ax.set_xlabel("X", fontsize = 17)
        ax.set_ylabel("Y", fontsize = 17)
        ax.set_title(prefixes[prefix], fontsize = fontsize)

        plt.axis('equal')
        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")
        ax.set_axisbelow(True)
        plt.tight_layout()

        plt.savefig(prefix + ".png")

    else:

        sys.exit("Ambiguous file list exists.")