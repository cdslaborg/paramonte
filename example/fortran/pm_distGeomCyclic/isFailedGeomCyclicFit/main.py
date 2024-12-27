#!/USR/BIN/ENV python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

linewidth = 2
fontsize = 17

marker ={ "CK" : "-"
        , "IK" : "."
        , "RK" : "-"
        }
label = [ r"probSuccess = .2, period = 20"
        , r"Histogram Fit"
        ]

fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
ax = plt.subplot()

dfrnd = pd.read_csv("isFailedGeomCyclicFit.rnd.txt", delimiter = ",", header = None)
for j in range(len(dfrnd.values[0,:])):
    plt.hist( dfrnd.values[:,j]
            , histtype = "stepfilled"
            , density = False
            , alpha = 0.5
            , bins = 75
            )
dffit = pd.read_csv("isFailedGeomCyclicFit.fit.txt", delimiter = ",", header = None)
for j in range(1, len(dffit.values[0,:])):
    plt.plot( dffit.values[:, 0]
            , dffit.values[:,j]
            , linewidth = 2
            , c = "red"
            )
    ax.legend   ( label
                , fontsize = fontsize
                )
plt.xticks(fontsize = fontsize - 2)
plt.yticks(fontsize = fontsize - 2)
#ax.set_xlim(0, 10)
ax.set_xlabel("Cyclic Geometric Random Value ( integer-valued )", fontsize = 17)
ax.set_ylabel("Count", fontsize = 17)
ax.set_title("{} Cyclic Geometric random values and histogram fit".format(len(dfrnd.values[:, 0])), fontsize = 17)

plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")

plt.savefig("isFailedGeomCyclicFit.png")