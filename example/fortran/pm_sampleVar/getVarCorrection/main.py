#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

linewidth = 2
fontsize = 17

pattern = "*.RK.txt"
fileList = glob.glob(pattern)
if len(fileList) == 1:

    df = pd.read_csv(fileList[0], delimiter = ",")

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    plt.plot( df.values[:, 0]
            , df.values[:,1:]
            , linewidth = linewidth
           #, color = "r"
            )
    ax.legend   ( list(df.columns.values[1:])
                , fontsize = fontsize
                )

    plt.xticks(fontsize = fontsize - 2)
    plt.yticks(fontsize = fontsize - 2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(df.columns.values[0], fontsize = 17)
    ax.set_ylabel("Correction", fontsize = 17)

    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    plt.savefig(fileList[0].replace(".txt",".png"))

elif len(fileList) > 1:

    sys.exit("Ambiguous file list exists.")