#!/USR/BIN/ENV python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import sys

linewidth = 2
fontsize = 17

for kind in ["sin.RK"]:

    file = glob.glob("*crd*"+kind+".txt")[0]
    df = pd.read_csv(file, delimiter = ",")

    #print(df.values)
    fig = plt.figure(figsize = (8, 6))
    ax = plt.subplot(1,1,1)
    ax.plot ( df.values[:, 0]
            , df.values[:, 1]
            , zorder = 1000
            )
    plt.minorticks_on()
    ax.set_xlabel("x", fontsize = 17)
    ax.set_ylabel("f(x)", fontsize = 17)
    ax.tick_params(axis = "x", which = "minor")
    ax.tick_params(axis = "y", which = "minor")
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.legend([file.split(".")[-3] + "(x)"], fontsize = fontsize)
    plt.tight_layout()
    plt.savefig(file.replace(".txt",".png"))

    file = glob.glob("*acf*"+kind+".txt")[0]
    df = pd.read_csv(file, delimiter = ",")
    fig = plt.figure(figsize = (8, 6))
    ax = plt.subplot(1,1,1)
    ax.plot ( df.values[:, 0]
            , df.values[:, 1]
            , zorder = 1000
            )
    plt.minorticks_on()
    ax.set_xlabel("Lag", fontsize = 17)
    ax.set_ylabel("acf(f)", fontsize = 17)
    ax.tick_params(axis = "x", which = "minor")
    ax.tick_params(axis = "y", which = "minor")
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    plt.tight_layout()
    plt.savefig(file.replace(".txt",".png"))