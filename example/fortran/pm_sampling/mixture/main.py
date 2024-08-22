#!/usr/bin/env python

import matplotlib.pyplot as plt
#import paramonte as pm
import pandas as pd
import numpy as np
import glob
import sys
import os

plt.rcParams['text.usetex'] = True
# Uncomment above line if LaTeX is installed on the system.
# On Ubuntu, this requires: sudo apt install dvipng texlive-latex-extra texlive-fonts-recommended cm-super
examname = os.path.basename(os.getcwd())

#### Visualize MCMC

linewidth = 2
fontsize = 17
files = glob.glob("./*_chain.txt")

for file in files:

    basename = file.split("_")[0]
    df = pd.read_csv(file, delimiter = ",")
    if "_chain.txt" in file:
        sindex = 7 # start column index.
    elif "_sample.txt" in file:
        sindex = 1 # start column index.
    else:
        sys.exit("Unrecognized simulation output file: " + file)

    # traceplot

    #print(df.values)
    fig = plt.figure(figsize = (8, 6))
    ax = plt.subplot(1,1,1)
    ax.plot ( range(len(df.values[:, 0]))
            , df.values[:, sindex:]
            , zorder = 1000
            )
    plt.minorticks_on()
    ax.set_xscale("log")
    ax.set_xlabel("MCMC Count", fontsize = 17)
    ax.set_ylabel("MCMC State", fontsize = 17)
    ax.tick_params(axis = "x", which = "major", labelsize = fontsize)
    ax.tick_params(axis = "y", which = "major", labelsize = fontsize)
    plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
    #ax.legend(df.columns.values, fontsize = fontsize)
    #plt.title(basename)
    plt.tight_layout()
    plt.savefig(basename + ".traceplot.png")

    # scatterplot

    if len(df.values[1, sindex:]) == 2:
        #print(df.values)
        fig = plt.figure(figsize = (8, 6))
        ax = plt.subplot(1,1,1)
        ax.scatter  ( df.values[:, sindex]
                    , df.values[:, sindex + 1]
                    , zorder = 1000
                    , s = 1
                    )
        plt.minorticks_on()
        ax.set_xscale("log")
        ax.set_xlabel(df.columns.values[sindex], fontsize = 17)
        ax.set_ylabel(df.columns.values[sindex + 1], fontsize = 17)
        ax.tick_params(axis = "x", which = "major", labelsize = fontsize)
        ax.tick_params(axis = "y", which = "major", labelsize = fontsize)
        plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
        #ax.legend(df.columns.values, fontsize = fontsize)
        #plt.title(basename)
        plt.tight_layout()
        plt.savefig(basename + ".scatterplot.png")

    # adaptation measure.

    #print(df.values)
    if "proposalAdaptation" in df.columns.values:
        if any(df["proposalAdaptation"].values != 0):
            fig = plt.figure(figsize = (8, 6))
            ax = plt.subplot(1,1,1)
            ax.scatter  ( range(len(df.values[:, 0]))
                        , df["proposalAdaptation"].values
                        , zorder = 1000
                        , c = "red"
                        , s = 1
                        )
            plt.minorticks_on()
            #ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("MCMC Count", fontsize = 17)
            ax.tick_params(axis = "x", which = "major", labelsize = fontsize)
            ax.tick_params(axis = "y", which = "major", labelsize = fontsize)
            ax.set_ylabel("proposalAdaptation", fontsize = 17)
            plt.grid(visible = True, which = "both", axis = "both", color = "0.85", linestyle = "-")
            #ax.legend(df.columns.values, fontsize = fontsize)
            #plt.title(basename)
            plt.tight_layout()
            plt.savefig(basename + ".proposalAdaptation.png")

    #sim = pm.ParaDRAM()
    #sample = sim.readSample("./out/", renabled = True)
    #sample[0].plot.grid()

######################
#### Visualize the fit
######################

# from Functions import mergeBins, findPlateauLocation

def mergeBins(counts, bins, threshold = 5):

    if len(bins)%2 == 0:
        halfway = int(np.ceil(len(bins)/2))
    else:
        halfway = int(np.floor(len(bins)/2))

    threshold = threshold
    leftCounts = []
    leftBins = [bins[0]]
    rightCounts = []
    rightBins = [bins[-1]]
    leftSumCounts = 0
    rightSumCounts = 0
    leftBinstart = 0
    rightBinstart = 0
    for i in range(0 , halfway): #goes to the middle of the distribution
       # print(counts[i],sumCounts)

        leftSumCounts += counts[i]



        if leftSumCounts >= threshold:
           # print("added value to list")
            leftCounts += [leftSumCounts]
            leftSumCounts = 0
            leftBins += [bins[i+1]]

        #if counts[-i-1]==3:print(rightSumCounts,-i-1)


        rightSumCounts += counts[-i-1]
        if rightSumCounts >= threshold:
            rightCounts += [rightSumCounts]
            rightSumCounts = 0
            rightBins += [bins[-i-2]]

    if len(bins)%2==0:
        if len(rightCounts) != 0:
            rightCounts.pop()
            rightBins.pop()
        leftBins.pop()
    else:
        rightCounts += [counts[halfway]]
        #rightBins += [bins[halfway]]

    newCounts = leftCounts + rightCounts[-1::-1]
    newBins = leftBins + rightBins[-1::-1]
    #print(len(newCounts),len(newBins))
    origionaldNdT = np.array(counts) / np.array([bins[i + 1] - bins[i] for i in range(0, len(bins) - 1) ])
    newdNdT = np.array(newCounts) / np.array([newBins[i + 1] - newBins[i] for i in range(0, len(newBins) - 1) ])

    return np.array(newCounts), np.array(newBins), np.array(origionaldNdT), np.array(newdNdT)

# Read the duration data.

data = pd.read_csv("data.csv", delimiter = ",")
durs = data["T90"].values

nbins = 50
binning = 10**np.linspace(np.log10(min(durs)), np.log10(max(durs)), nbins)
countsAndBins = np.histogram(durs, bins = binning)
counts = countsAndBins[0]
bins = countsAndBins[1]
bins = bins.tolist()
newCounts, newBins, origionaldNdT, dNdT = mergeBins(counts, binning)

# get the middle of the bin edges for location of uncertainties.
#midBins = [newBins[i] + (newBins[i + 1] - newBins[i]) / 2 for i in range(0, len(newBins) - 1)]
# use algorithm to find the longest plateau within one std of right most point
#plateau, binEdge, plateauErrors = findPlateauLocation(dNdT, dNdTerrors, newBins)
# find the best value for the plateau
#bestParams = opt.minimize(chisqfunc, [plateau[-1]], args = (plateau, plateauErrors))
#x = np.logspace(np.log10(min(binEdge)),np.log10(max(binEdge)),len(midBins))
#plt.plot(x,(bestParams["x"]).tolist()*len(x),color = 'k',lw = 2)
#x = np.logspace(np.log10(min(binEdge)), np.log10(max(binEdge)), len(midBins))

histFiles = glob.glob("./*.hist")

for histFile in histFiles:

    fig = plt.figure()
    plt.stairs  ( dNdT
                , newBins
                , linewidth = 2
                , baseline = None
                , color = "green"
                , label = 'BATSE'
                )
    plt.xscale('log')
    plt.yscale('log')

    # Read the fitting data to plot.

    model = pd.read_csv(histFile, delimiter = ",")
    x = np.exp(model["logx"].values)
    y = len(durs) * np.exp(model["logPDF"].values) / x
    plt.grid(which = "both")
    plt.plot(x, y, c = "magenta", label = histFile)
    plt.xlabel(r"$T_{90} ~ [s]$", fontsize = fontsize)
    plt.ylabel(r"$dN ~ / ~ dT_{90} ~ [1 / s]$", fontsize = fontsize)
    #plt.legend(fontsize = 13)
    plt.xticks(size = fontsize)
    plt.yticks(size = fontsize)
    plt.ylim(7*10**-3, 2*10**3)
    plt.tight_layout()
    plt.savefig(histFile + ".png")
    #plt.show()