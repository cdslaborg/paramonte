#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

fontsize = 17

kind = "IK"
label = [ r"$\lambda = 0.1$"
        , r"$\lambda = 1.0$"
        , r"$\lambda = 4.0$"
        , r"$\lambda = 10.$"
        ]

pattern = "*." + kind + ".txt"
fileList = glob.glob(pattern)
if len(fileList) == 1:

    df = pd.read_csv(fileList[0], delimiter = " ", header = None)

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    for i in range(1,len(df.values[0,:]+1)):

            plt.plot( df.values[:, 0]
                    , df.values[:,i]
                    , marker = "."
                    )

    plt.xticks(fontsize = fontsize - 2)
    plt.yticks(fontsize = fontsize - 2)
    ax.set_xlim([0, None])
    ax.set_xlabel("Count", fontsize = fontsize)
    ax.set_ylabel("CDF", fontsize = fontsize)

    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    ax.legend   ( label
                , fontsize = fontsize
               #, loc = "center left"
               #, bbox_to_anchor = (1, 0.5)
                )

    plt.savefig(fileList[0].replace(".txt",".png"))

else:

    sys.exit("Ambiguous file list exists.")