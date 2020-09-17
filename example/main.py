#!/usr/bin/env python
#!C:\ProgramData\Anaconda3\python.exe
####################################################################################################################################
##  
##  Description:
##      +   Run the Monte Carlo sampler of the ParaMonte library given the input log-target density function `getLogFunc()`.
##  Output:
##      +   The simulation output files.
##  Author:
##      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
##  Visit:
##      +   https://www.cdslab.org/paramonte
##
####################################################################################################################################

import os
import numpy as np
import paramonte as pm
from logfunc import getLogFunc, NDIM

pmpd = pm.ParaDRAM() # define a ParaDRAM sampler instance

pmpd.runSampler ( ndim = NDIM               # number of dimensions of the objective function
                , getLogFunc = getLogFunc   # the objective function: multivariate normal distribution
                # NOTE: inputFile is optional: all simulation specifications can be set as attributes of pmpd.spec
                # NOTE: However, if specified, then all simulation specification will be read exclusively from the inputFile.
                , inputFile = os.path.dirname(os.path.abspath(__file__)) + "/paramonte.in"
                )
