#!/USR/BIN/ENV python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 14

methods = ["getLogFactorial", "getFactorial"]
labels = [label+"(Whole Number)" for label in methods]

df = pd.read_csv("main.out")

####################################################################################################################################
#### Plot the runtimes.
####################################################################################################################################

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

for method in methods:
    plt.plot( df["Point"].values
            , df[method].values
            , linewidth = 2
            , linestyle = "-"
            , marker = "*"
            )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Whole Number", fontsize = fontsize)
ax.set_ylabel("Runtime [ seconds ]", fontsize = fontsize)
ax.set_title("getLogFactorial() vs. getFactorial()\nLower is better.", fontsize = fontsize)
ax.set_xscale("log")
ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( labels
           #, loc='center left'
           #, bbox_to_anchor=(1, 0.5)
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.getFactorial_vs_getLogFactorial.runtime.png")

####################################################################################################################################
#### Plot the runtime ratios.
####################################################################################################################################

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

plt.plot( df["Point"].values
        , np.ones(len(df["Point"].values))
        , linestyle = "-"
        , marker = "*"
       #, color = "black"
        , linewidth = 2
        )
plt.plot( df["Point"].values
        , df["getFactorial"].values / df["getLogFactorial"].values
        , linestyle = "-"
        , marker = "*"
        , linewidth = 2
        )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Whole Number", fontsize = fontsize)
ax.set_ylabel("Runtime Ratio", fontsize = fontsize)
ax.set_title("getFactorial() / getLogFactorial(), Lower means faster.\nLower than 1 means faster than getLogFactorial().", fontsize = fontsize)
ax.set_xscale("log")
#ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( labels
           #, bbox_to_anchor = (1, 0.5)
           #, loc = "center left"
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.getFactorial_vs_getLogFactorial.runtime.ratio.png")
