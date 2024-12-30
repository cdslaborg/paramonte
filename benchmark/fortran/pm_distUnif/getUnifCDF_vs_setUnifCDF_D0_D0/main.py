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

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

for method in methods:
    plt.hist( np.log10(df[method].values)
            , histtype = "stepfilled"
            , density = True
            , alpha = 0.7
            , bins = 30
            )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("$\log_{10}$ ( Runtime [ seconds ] )", fontsize = fontsize)
ax.set_ylabel("Count", fontsize = fontsize)
ax.set_title("Rank-0 Runtime:\ngetUnifCDF vs. setUnifCDF.\nLower is better.", fontsize = fontsize)
#ax.set_xscale("log")
#ax.set_yscale("log")
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
plt.savefig("benchmark.getUnifCDF_vs_setUnifCDF_D0_D0.runtime.png")

####################################################################################################################################
#### Plot the runtime ratios.
####################################################################################################################################

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

plt.hist( np.log10(df["getUnifCDF"].values / df["setUnifCDF"].values)
        , histtype = "stepfilled"
        , density = True
        , alpha = 0.7
        , bins = 30
        )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel(r"$\log_{10}$ ( Runtime Ratio [ seconds ] )", fontsize = fontsize)
ax.set_ylabel("Count", fontsize = fontsize)
ax.set_title(r"""$\log_{10}$ ( Rank-0 Runtime Ratio ): getUnifCDF to setUnifCDF.
A value < $\log_{10}(1)$ implies better performance of getUnifCDF.""", fontsize = fontsize)
#ax.set_xscale("log")
#ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")

plt.tight_layout()
plt.savefig("benchmark.getUnifCDF_vs_setUnifCDF_D0_D0.runtime.ratio.png")
