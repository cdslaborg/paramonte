#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

linewidth = 2
fontsize = 17

for kind in ["RK"]:

    pattern = "*." + kind + ".txt"
    fileList = glob.glob(pattern)

    for file in fileList:

        df = pd.read_csv(file, delimiter = ",", header = None)

        # start with a square Figure
        fig = plt.figure(figsize=(8, 8))

        plt.rcParams.update({'font.size': fontsize - 2})
        ax = plt.subplot()

        ax.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "y", which = "minor")
        ax.tick_params(axis = "x", which = "minor")

        # the scatter plot:
        ax.plot     ( df.values[:, 0]
                    , df.values[:,1]
                    , linewidth = 2
                    , color = "black"
                    )
        ax.scatter  ( df.values[:, 0]
                    , df.values[:,1]
                    , s = 8
                    , zorder = 1000
                    , color = "red"
                    )

        ax.set_xlabel("X", fontsize = 17)
        ax.set_ylabel("Y", fontsize = 17)

        plt.savefig(file.replace(".txt",".png"))
