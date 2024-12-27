#!/USR/BIN/ENV python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import os
dirname = os.path.basename(os.getcwd()) 

fontsize = 14

####################################################################################################################################
#### Plot the runtime ratios for normal array.
####################################################################################################################################

df = pd.read_csv("main.normal.out")
colnames = list(df.columns.values)

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

for colname in colnames[1:]:
    plt.plot( df[colnames[0]].values
           #, np.ones(len(df[colnames[0]].values))
            , df[colname].values / df[colnames[1]].values
            , linewidth = 2
           #, linestyle = "-"
            )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size", fontsize = fontsize)
ax.set_ylabel("Runtime Ratio ( W.R.T. cenabled = .false. )", fontsize = fontsize)
ax.set_title("Performance when the input arrays\ncause no underflow instances.\nLower is better.", fontsize = fontsize)
ax.set_xscale("log")
#ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( colnames[1:]
           #, loc='center left'
           #, bbox_to_anchor=(1, 0.5)
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark." + dirname + ".normal.png")

####################################################################################################################################
#### Plot the runtime ratios for the array causing many underflows.
####################################################################################################################################

df = pd.read_csv("main.underflow.out")
colnames = list(df.columns.values)

ax = plt.figure(figsize = 1.25 * np.array([6.4,4.6]), dpi = 200)
ax = plt.subplot()

for colname in colnames[1:]:
    plt.plot( df[colnames[0]].values
           #, np.ones(len(df[colnames[0]].values))
            , df[colname].values / df[colnames[1]].values
            , linewidth = 2
           #, linestyle = "-"
            )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size", fontsize = fontsize)
ax.set_ylabel("Runtime Ratio W.R.T. {})".format(colnames[1]), fontsize = fontsize)
ax.set_title("Performance when the input arrays\ncause many underflow instances.\nLower is better.", fontsize = fontsize)
ax.set_xscale("log")
#ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( colnames[1:]
           #, loc='center left'
           #, bbox_to_anchor=(1, 0.5)
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark." + dirname + ".underflow.png")
