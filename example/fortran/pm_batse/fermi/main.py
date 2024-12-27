#!/USR/BIN/ENV python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

fontsize = 17

pattern = "fermi.z.out"
fileList = glob.glob(pattern)
if len(fileList) == 1:

    df = pd.read_csv(fileList[0], delimiter = ",")

    #### plot histogram

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    plt.hist( df["zmode"].values
            , histtype = "stepfilled"
            , alpha = 0.5
            , bins = 50
            )

    plt.xticks(fontsize = fontsize - 2)
    plt.yticks(fontsize = fontsize - 2)
    ax.set_xlabel("redshift", fontsize = fontsize)
    ax.set_ylabel("Count", fontsize = fontsize)
    #ax.set_xlim([0, 5])

    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    ax.legend   ( ["{} Fermi Redshifts".format(len(df["zmode"].values))]
                , fontsize = fontsize
               #, loc = "center left"
               #, bbox_to_anchor = (1, 0.5)
                )

    plt.savefig(fileList[0].replace(".out",".hist.png"))

    #### plot comparison

    dfz = df[df["zreal"] > 0]

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    ax.plot ( [0, 1] # [.1, .1]
            , [0, 1] # [np.max(dfz.values), np.max(dfz.values)]
            , transform = ax.transAxes
            , linewidth = 1.5
            , color = "black"
            , zorder = 1000
            )

    plt.scatter ( dfz["zreal"].values
                , dfz["zmode"].values
                , zorder = 1001
                , s = 20
                )

    plt.xticks(fontsize = fontsize - 2)
    plt.yticks(fontsize = fontsize - 2)
    ax.set_xlabel("Reference redshift", fontsize = fontsize)
    ax.set_ylabel("Predicted redshift", fontsize = fontsize)
    lb, ub = .9 * np.min(dfz.values), np.max(dfz.values) * 1.1
    ax.set_xlim([lb, ub])
    ax.set_ylim([lb, ub])
    ax.set_xscale("log")
    ax.set_yscale("log")

    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    ax.legend   ( ["Equality", "{} Fermi Redshifts".format(len(dfz["zmode"].values))]
                , fontsize = fontsize
               #, loc = "center left"
               #, bbox_to_anchor = (1, 0.5)
                )

    plt.savefig(fileList[0].replace(".out",".scatter.png"))

else:

    sys.exit("Ambiguous file list exists.")
