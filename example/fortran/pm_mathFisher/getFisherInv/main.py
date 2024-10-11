#!/usr/bin/env python
#pip install paramonte

import os
examname = os.path.basename(os.getcwd())

##################################################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 17

df = pd.read_csv(examname + ".csv", delimiter = ",")

fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
ax = plt.subplot()

plt.plot( df.values[:, 1]
        , df.values[:, 1]
        , linestyle = "-"
        , linewidth = 2
        )
for colname in df.columns[1:]:
    plt.plot( df.values[:, 0]
            , df[colname].values
            , linestyle = "-"
            , linewidth = 2
            )

labels = ["Equality"]
for colname in list(df.columns[1:]): labels.append(colname)
ax.legend   ( labels
            , fontsize = fontsize
           #, loc = "center left"
           #, bbox_to_anchor = (1, 0.5)
            )

ax.set_xlabel(df.columns[1], fontsize = fontsize)
ax.set_ylabel(df.columns[0], fontsize = fontsize)

plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)

plt.tight_layout()
plt.savefig(examname + ".png")