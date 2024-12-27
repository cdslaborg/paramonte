#!/USR/BIN/ENV python
#pip install paramonte

import os
examname = os.path.basename(os.getcwd())
modelName = examname[-3:]

##################################################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 17

df = pd.read_csv(examname + ".csv", delimiter = ",")

fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
ax = plt.subplot()

for colname in df.columns[1:]:
    plt.plot( df.values[:, 0]
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
plt.savefig(examname + ".z.png")

##################################################################

#### Visualize MCMC

import glob
linewidth = 2
fontsize = 17
files = glob.glob("./*_sample.txt")

for file in files:

    #basename = file.split("_")[0]
    df = pd.read_csv(file, delimiter = ",")
    if "_chain.txt" in file:
        sindex = 7 # start column index.
    elif "_sample.txt" in file:
        sindex = 1 # start column index.
    else:
        sys.exit("Unrecognized simulation output file: " + file)

    # histogram

    for histname in ["z", "logzplus1"]:

        fig = plt.figure(figsize = (8, 6))
        ax = plt.subplot(1,1,1)
        if histname == "z":
            ax.hist(df.values[:, sindex:])
        elif histname == "logzplus1":
            ax.hist(np.log(df.values[:, sindex:]))
        else:
            sys.exit("Unrecognized histogram name: " + histname)
        plt.minorticks_on()
        #ax.set_xscale("log")
        ax.set_ylabel("Count", fontsize = 17)
        ax.set_xlabel("log( z + 1 )", fontsize = 17)
        ax.tick_params(axis = "x", which = "minor")
        ax.tick_params(axis = "y", which = "minor")
        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        plt.title("Histogram of the simulated redshifts from\nthe rate density model of " + modelName)
        plt.tight_layout()
        plt.savefig(examname + "." + histname + ".sample.png")
