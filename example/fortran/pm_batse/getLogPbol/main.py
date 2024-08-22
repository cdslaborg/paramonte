#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

fontsize = 17

kind = "RK"

pattern = "*." + kind + ".txt"
fileList = glob.glob(pattern)
if len(fileList) == 1:

    df = pd.read_csv(fileList[0], delimiter = " ")

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    plt.plot( df.values[:,1]
            , df.values[:, 0]
            , linewidth = 2
            )

    plt.plot( df.values[:,2]
            , df.values[:, 0]
            , linewidth = 2
            )

    plt.xticks(fontsize = fontsize - 2)
    plt.yticks(fontsize = fontsize - 2)
    ax.set_xlabel("Bolometric energy flux range", fontsize = fontsize)
    ax.set_ylabel("Epeak [ keV ]", fontsize = fontsize)
    ax.set_xscale("log")
    ax.set_yscale("log")

    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    ax.legend   ( ["getLogPbol()", "setBandPhoton()"]
                , fontsize = fontsize
               #, loc = "center left"
               #, bbox_to_anchor = (1, 0.5)
                )

    plt.savefig(fileList[0].replace(".txt",".png"))

else:

    sys.exit("Ambiguous file list exists.")
