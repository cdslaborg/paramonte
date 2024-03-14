#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 14

methods = ["default", "iseqArg"]

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
ax.set_title("Finding pattern with external `iseq()` vs. default comparison.\nLower is better.", fontsize = fontsize)
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
plt.savefig("benchmark.setLoc-iseqArg_vs_default.runtime.png")

####################################################################################################################################
#### Plot the runtime ratios.
####################################################################################################################################

ax = plt.figure(figsize = 1.25*np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

plt.plot( df["arraySize"].values
        , np.ones(len(df["arraySize"].values))
       #, linestyle = "--"
       #, color = "black"
        , linewidth = 2
        )
plt.plot( df["arraySize"].values
        , df["iseqArg"].values / df["default"].values
        , linewidth = 2
        )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size", fontsize = fontsize)
ax.set_ylabel("Runtime compared to iseqArg", fontsize = fontsize)
ax.set_title("Runtime Ratio: Comparison with `iseq()` / default.\nLower means faster. Lower than 1 means faster than finding with `iseq()`.", fontsize = fontsize)
ax.set_xscale("log")
#ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( methods
           #, bbox_to_anchor = (1, 0.5)
           #, loc = "center left"
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.setLoc-iseqArg_vs_default.runtime.ratio.png")
