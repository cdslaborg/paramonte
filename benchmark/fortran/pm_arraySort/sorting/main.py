#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

fontsize = 15

colnames =  [ "setSortedQsorti"
            , "setSortedQsortr"
            , "setSortedQsortrdp"
            , "setSortedBubble"
            , "setSortedHeapi"
            , "setSortedInsertionl"
            , "setSortedInsertionb"
            , "setSortedMerge"
            , "setSortedSelection"
            , "setSortedShell"
            ]

dtypes = ["random", "sorted"]

####################################################################################################################################

df = {}

for dtype in dtypes:

    df[dtype] = pd.read_csv("main.{}.out".format(dtype))
    colnames = list(df[dtype].columns.values[1:])

    ax = plt.figure(figsize = 1.25 * np.array([9.4,4.8]), dpi = 200)
    ax = plt.subplot()

    for colname in colnames:
        plt.plot( df[dtype].values[:, 0]
                , df[dtype][colname].values / df[dtype]["setSortedQsorti"].values
                , linewidth = 2
                )

    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    ax.set_xlabel("Array Size (# elements)", fontsize = fontsize)
    ax.set_ylabel("Time normalized to qsorti timing", fontsize = fontsize)
    ax.set_title("Benchmarking of sorting algorithms compared to qsorti.\nInput data is " + r"$\mathrm{\bf{" + dtype + "}}$. Lower is better.", fontsize = fontsize)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.minorticks_on()
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    ax.tick_params(axis = "y", which = "minor")
    ax.tick_params(axis = "x", which = "minor")
    ax.legend   ( colnames
                , loc='center left'
                , bbox_to_anchor=(1, 0.5)
                , fontsize = fontsize
                )

    plt.tight_layout()
    plt.savefig("benchmark.sorting.{}.png".format(dtype))

####################################################################################################################################

ax = plt.figure(figsize = 1.25 * np.array([9.4,4.8]), dpi = 200)
ax = plt.subplot()

for colname in colnames:
    plt.plot( df["random"].values[:, 0]
            , df["sorted"][colname].values / df["random"][colname].values
            , linewidth = 2
            )

plt.xticks(fontsize = fontsize)
plt.yticks(fontsize = fontsize)
ax.set_xlabel("Array Size (# elements)", fontsize = fontsize)
ax.set_ylabel("Sorting Time Ratio ( Sorted / Random Input Array )", fontsize = fontsize)
ax.set_title("Sorted / Random input-arrays sorting-time ratio. Lower is better.", fontsize = fontsize)
ax.set_xscale("log")
ax.set_yscale("log")
plt.minorticks_on()
plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
ax.tick_params(axis = "y", which = "minor")
ax.tick_params(axis = "x", which = "minor")
ax.legend   ( colnames
            , loc='center left'
            , bbox_to_anchor=(1, 0.5)
            , fontsize = fontsize
            )

plt.tight_layout()
plt.savefig("benchmark.sorting.random.sorted.ratio.png")
