#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

linewidth = 2
fontsize = 17

for kind in ["neimean", "neinear", "neinext", "neiprev", "piwilin", "monopol", "rungeEffect"]:

    crd = pd.read_csv(glob.glob("*."+kind+".crd.txt")[0], delimiter = ",")
    pattern = "*."+kind+".interp.txt"
    fileList = glob.glob(pattern)

    for file in fileList:

        df = pd.read_csv(file, delimiter = ",")

        # start with a square Figure
        fig = plt.figure(figsize = (8, 6))
        ax = plt.subplot(1,1,1)

        ax.scatter  ( df.values[:, 0]
                    , df.values[:, 1]
                    , zorder = 1000
                    , c = "black"
                    , s = 8
                    )
        ax.scatter  ( crd.values[:, 0]
                    , crd.values[:,1]
                    , zorder = 1000
                    , c = "red"
                    , s = 20
                    )

        plt.minorticks_on()
        ax.set_xlabel("X", fontsize = 17)
        ax.set_ylabel("Y", fontsize = 17)
        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        ax.tick_params(axis = "x", which = "minor")
        ax.tick_params(axis = "y", which = "minor")
        ax.legend([file.split(".")[-3], "nodes"], fontsize = fontsize)

        plt.tight_layout()
        plt.savefig(file.replace(".txt",".png"))