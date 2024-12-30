#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys
import os

fontsize = 17

parent = os.path.basename(os.path.dirname(__file__))
pattern = parent + "*.txt"

fileList = glob.glob(pattern)
for file in fileList:

    df = pd.read_csv(file, delimiter = " ")

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    for i in range(1,len(df.values[0,:]+1)):

            plt.plot( df.values[:, 0]
                    , df.values[:,i]
                    , linewidth = 2
                    )

    plt.xticks(fontsize = fontsize - 2)
    plt.yticks(fontsize = fontsize - 2)
    ax.set_xlabel("x", fontsize = fontsize)
    if "abserr" in file:
        nbit = file.split(".")[1][2:]
        label = [r"x - erf(erfinv(x))"]
        ax.set_ylabel("{}-bits Absolute Error".format(nbit), fontsize = fontsize)
    else:
        label = [r"erfinv(x)"]
        ax.set_ylabel("Inverse Error Function", fontsize = fontsize)

    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    ax.legend   ( label
                , fontsize = fontsize
               #, loc = "center left"
               #, bbox_to_anchor = (1, 0.5)
                )

    plt.savefig(file.replace(".txt",".png"))