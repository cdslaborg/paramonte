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

        # definitions for the axes
       #left, width = 0.1, 0.65
       #bottom, height = 0.1, 0.65
       #spacing = 0.015

        # start with a square Figure
        fig = plt.figure() # figsize = (8, 6)

        plt.rcParams.update({'font.size': fontsize - 2})
        ax = plt.subplot()
       #ax = fig.add_axes([left, bottom, width, height]) # scatter plot
       #ax_histx = fig.add_axes([left, bottom + height + spacing, width, 0.2], sharex = ax) # histx
       #ax_histy = fig.add_axes([left + width + spacing, bottom, 0.2, height], sharey = ax) # histy
       #
       #for axes in [ax, ax_histx, ax_histy]:
       #    axes.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
       #    axes.tick_params(axis = "y", which = "minor")
       #    axes.tick_params(axis = "x", which = "minor")
       #
       ## no labels
       #ax_histy.tick_params(axis = "y", labelleft = False)
       #ax_histx.tick_params(axis = "x", labelbottom = False)

        # the scatter plot:
        for i in range(0, len(df.values[0,:]), 2):
            ax.scatter  ( df.values[:,i]
                        , df.values[:,i+1]
                        , s = 8
                        , zorder = 1000
                        )

       #ax_histx.hist(df.values[:, 0], bins = 50, zorder = 1000)
       #ax_histy.hist(df.values[:, 1], bins = 50, orientation = "horizontal", zorder = 1000)

        #ax.set_aspect("equal")
        ax.set_xlim([0, df.values[-1,:].max()])
        ax.set_xlabel("ndim", fontsize = 17)
        ax.set_ylabel("Volume of ndim unit-ball", fontsize = 17)
        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        plt.minorticks_on()

        plt.tight_layout()
        plt.savefig(file.replace(".txt",".png"))
