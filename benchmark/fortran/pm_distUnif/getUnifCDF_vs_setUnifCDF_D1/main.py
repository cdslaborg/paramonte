#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 14

methods = ["setUnifCDF", "getUnifCDF"]

df = pd.read_csv("main.out")

####################################################################################################################################
#### Plot the runtimes.
####################################################################################################################################

ax = plt.figure(figsize = 1.25*np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

for method in methods:
    plt.plot( df["arraySize"].values
            , df[method].values
            , linewidth = 2
            )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size", fontsize = fontsize)
ax.set_ylabel("Runtime [ seconds ]", fontsize = fontsize)
ax.set_title("Rank-1 Runtime:\ngetUnifCDF vs. setUnifCDF.\nLower is better.", fontsize = fontsize)
ax.set_xscale("log")
ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( methods
           #, loc='center left'
           #, bbox_to_anchor=(1, 0.5)
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.getUnifCDF_vs_setUnifCDF_D1.runtime.png")

####################################################################################################################################
#### Plot the runtime ratios.
####################################################################################################################################

ax = plt.figure(figsize = 1.25*np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

plt.plot( df["arraySize"].values
        , np.ones( len(df["arraySize"].values) )
        , linestyle = "-"
        , linewidth = 2
       #, color = "black"
        )
plt.plot( df["arraySize"].values
        , df["getUnifCDF"].values / df["setUnifCDF"].values
        , linewidth = 2
       #, color = "r"
        )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size", fontsize = fontsize)
ax.set_ylabel("Runtime Ratio", fontsize = fontsize)
ax.set_title("""Rank-1 Runtime Ratio: getUnifCDF to setUnifCDF.
A value < 1 implies better performance of getUnifCDF.""", fontsize = fontsize)
ax.set_xscale("log")
#ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( ["setUnifCDF", "getUnifCDF"]
           #, loc='center left'
           #, bbox_to_anchor=(1, 0.5)
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.getUnifCDF_vs_setUnifCDF_D1.runtime.ratio.png")
