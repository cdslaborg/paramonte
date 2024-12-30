#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys
import os

linewidth = 2
fontsize = 17

parent = os.path.basename(os.path.dirname(__file__))
pattern = parent + "*.txt"
files = glob.glob(pattern)

for file in files:

    df = pd.read_csv(file, delimiter = "\t")

    #print(df.values)
    fig = plt.figure(figsize = (8, 6))
    ax = plt.subplot(1,1,1)
    ax.plot ( df.values[:, 0]
            , df.values[:, 1]
            , zorder = 1000
            )
    #ax.scatter  ( df.values[:, 0]
    #            , df.values[:, 1]
    #            , zorder = 1000
    #            , s = 1
    #            )
    plt.minorticks_on()
    ax.set_xlabel(df.columns.values[0], fontsize = 17)
    ax.set_ylabel(df.columns.values[1], fontsize = 17)
    ax.tick_params(axis = "x", which = "minor")
    ax.tick_params(axis = "y", which = "minor")
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.legend(["Duluth, MN"], fontsize = fontsize)
    plt.tight_layout()
    plt.savefig(file.replace(".txt",".png"))
