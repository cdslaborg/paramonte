#!/USR/BIN/ENV python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 15

####################################################################################################################################

dtypes = ["random", "sorted"]

for dtype in dtypes:

    df = pd.read_csv("main.{}.out".format(dtype))

    #### Plot timings.

    ax = plt.figure(figsize = 1.25 * np.array([9.4,4.8]), dpi = 200)
    ax = plt.subplot()

    for colname in df.columns[1:]:
        plt.plot( df.values[:, 0]
                , df[colname].values
                , linewidth = 2
                )

    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    ax.set_xlabel("Array Size (# elements)", fontsize = fontsize)
    ax.set_ylabel("Time [ seconds ]", fontsize = fontsize)
    ax.set_title("timing of setSortedArray() vs. setSortedIndex(). Input data is {}.\n Lower is better.".format(dtype), fontsize = fontsize)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.minorticks_on()
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")
    ax.legend   ( df.columns[1:]
               #, loc='center left'
               #, bbox_to_anchor=(1, 0.5)
                , fontsize = fontsize
                )

    plt.tight_layout()
    plt.savefig("benchmark.sorting_vs_indexing.{}.png".format(dtype))

    #### Plot timing ratios.

    ax = plt.figure(figsize = 1.25 * np.array([9.4,4.8]), dpi = 200)
    ax = plt.subplot()

    start = 0

    plt.plot( df.values[:, 0]
            , np.ones(len(df.values[:, 0]))
            , linestyle = "--"
            , linewidth = 2
           #, color = "black"
            )
    for colname in df.columns[2:]:
        plt.plot( df.values[start:,0]
                , df[colname].values[start:] / df["setSortedArray"].values[start:]
                , linewidth = 2
                )

    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    ax.set_xlabel("Array Size (# elements)", fontsize = fontsize)
    ax.set_ylabel("Time Ratio ( setSortedIndex() / setSortedArray() )", fontsize = fontsize)
    ax.set_title("Ratio of timing of setSortedIndex() to setSortedArray(). Input data is {}.".format(dtype), fontsize = fontsize)
    ax.set_xscale("log")
    #ax.set_yscale("log")
    plt.minorticks_on()
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")
    ax.legend   ( ["setSortedArray"] + list(df.columns[2:])
               #, bbox_to_anchor=(1, 0.5)
               #, loc='center left'
                , fontsize = fontsize
                )

    plt.tight_layout()
    plt.savefig("benchmark.sorting_vs_indexing.{}.ratio.png".format(dtype))

####################################################################################################################################
