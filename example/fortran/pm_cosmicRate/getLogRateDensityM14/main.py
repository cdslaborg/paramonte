#!/usr/bin/env python
#pip install paramonte

modelName = "M14"
filePrefix = "getLogRateDensity" + modelName

##################################################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 17

df = pd.read_csv(filePrefix + ".csv", delimiter = ",")

fig = plt.figure(figsize = 1.25*np.array([6.4,4.8]), dpi = 200)
ax = plt.subplot()

for colname in df.columns[1:]:
    plt.plot( df.values[:,0]
            , df[colname].values
            , linestyle = "-"
            , linewidth = 2
            )

labels = []
for colname in list(df.columns[1:]): labels.append(modelName + ": " + colname)
ax.legend   ( labels
            , fontsize = fontsize
           #, loc = "center left"
           #, bbox_to_anchor = (1, 0.5)
            )

ax.set_xlabel(df.columns[0], fontsize = fontsize)
ax.set_ylabel("Unnormalized Density", fontsize = fontsize)

plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)

plt.tight_layout()
plt.savefig(filePrefix + ".z.png")

##################################################################

if False:
    import paramonte as pm
    sim = pm.Paradram()
    sim.readSample()
    sim.sampleList[0].plot.histplot()
    sim.sampleList[0].plot.histplot.funcout.axes.set_xlim([0,3])
    sim.sampleList[0].plot.histplot.funcout.axes.set_xlabel("log( z + 1 )")
    import matplotlib.pyplot as plt
    plt.title("Histogram of the simulated redshifts from\nthe rate density model of " + modelName)
    plt.tight_layout()
    sim.sampleList[0].plot.histplot.savefig(fname = filePrefix + ".logzplus1.sample.png")

    # Plot the histogram of the corresponding redshift values.

    import numpy as np
    sim.sampleList[0].df["z"] = np.exp(sim.sampleList[0].contents.logzplus1) - 1
    sim.sampleList[0].plot.histplot(xcolumns = "z")
    sim.sampleList[0].plot.histplot.funcout.axes.set_xlim([0,12])
    sim.sampleList[0].plot.histplot.funcout.axes.set_xlabel("redshift: z")
    import matplotlib.pyplot as plt
    plt.title("Histogram of the simulated redshifts from\nthe rate density model of " + modelName)
    plt.tight_layout()
    sim.sampleList[0].plot.histplot.savefig(fname = filePrefix + ".z.sample.png")
