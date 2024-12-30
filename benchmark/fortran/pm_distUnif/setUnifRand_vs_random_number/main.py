#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 14

df = pd.read_csv("main.out")

####################################################################################################################################
#### Plot the runtimes.
####################################################################################################################################

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

plt.plot( df.values[:, 0]
        , df.values[:,1:]
        , linewidth = 2
        )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size", fontsize = fontsize)
ax.set_ylabel("Runtime [ seconds ]", fontsize = fontsize)
ax.set_title("Runtime:\nsetUnifRand() vs. random_number().\nLower is better.", fontsize = fontsize)
ax.set_xscale("log")
ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( list(df.columns.values[1:])
           #, loc='center left'
           #, bbox_to_anchor=(1, 0.5)
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.setUnifRand_vs_random_number.runtime.png")

####################################################################################################################################
#### Plot the runtime ratios.
####################################################################################################################################

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

# baseline

plt.plot( df.values[:, 0]
        , np.ones( len(df["arraySize"].values) )
        , linestyle = "-"
        , linewidth = 2
        )
for colname in df.columns.values[2:]:
    plt.plot( df.values[:, 0]
            , df[colname].values / df.values[:, 1]
            , linewidth = 2
            )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size", fontsize = fontsize)
ax.set_ylabel("Runtime Ratio", fontsize = fontsize)
ax.set_title("""Uniform RNG Runtime Ratio: setUnifRand() to random_number().
A value < 1 implies better performance of setUnifRand().""", fontsize = fontsize)
ax.set_xscale("log")
#ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( list(df.columns.values[1:])
           #, bbox_to_anchor=(1, 0.5)
           #, loc='center left'
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.setUnifRand_vs_random_number.runtime.ratio.png")