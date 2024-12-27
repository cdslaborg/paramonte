#!/USR/BIN/ENV python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import re
fontsize = 15

# Sort files by digit.

files = sorted(glob.glob("main.norm.*.out"), key=lambda x:float(re.findall("(\d+)",x)[0]))

for file in files:

    df = pd.read_csv(file)

    fig = plt.figure(figsize = 1.25 * np.array([6.4, 4.8]), dpi = 200)
    ax = plt.subplot()

    plt.plot( df["sample"].values
            , df["ecdf"].values
            , linewidth = 3
            , color = "black"
            )

    plt.plot( df["sample"].values
            , df["lcdf"].values
            , linewidth = 2
            , color = "r"
            )

    plt.plot( df["sample"].values
            , df["ucdf"].values
            , linewidth = 2
            , color = "r"
            )

    ax.set_xlabel("X", fontsize = fontsize)
    ax.set_ylabel("Empirical Cumulative Distribution", fontsize = fontsize)
    ax.set_title("ECDF of a Normal sample of Size {}".format(len(df.index)), fontsize = fontsize)

    plt.minorticks_on()
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")

    ax.legend   ( ["ECDF", "CI 1%","CI 99%"]
                , bbox_to_anchor = (1, 0.5)
                , fontsize = fontsize
                , loc = "upper left"
                )

    plt.tight_layout()

    plt.savefig(file + ".png")
