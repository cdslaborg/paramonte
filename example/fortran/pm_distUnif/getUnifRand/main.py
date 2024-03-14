#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

rand = pd.read_csv("main.unif.rand.txt")

fig = plt.figure(figsize = 1.25*np.array([6.4,4.8]), dpi = 200)
ax = plt.subplot()

plt.plot( rand.values[:,0]
        , rand.values[:,1]
        , "."
        , color = "r"
        )

ax.set_xlabel("Real Random Component", fontsize = 17)
ax.set_ylabel("Imaginary Random Component", fontsize = 17)
ax.set_title("A scatter plot of Complex-valued random numbers".format(len(rand)), fontsize = 17)

plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")

plt.savefig("getUnifRand.CK.png")
