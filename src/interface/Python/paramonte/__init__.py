#!/usr/bin/python
# Author:  Amir Shahmoradi
# Contact: shahmoradi@utexas.edu
#**********************************************************************************************************************************
#**********************************************************************************************************************************
#
#  ParaMonte: plain powerful parallel Monte Carlo library.
#
#  Copyright (C) 2012-present, The Computational Data Science Lab
#
#  This file is part of ParaMonte library. 
#
#  ParaMonte is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, version 3 of the License.
#
#  ParaMonte is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
#
#**********************************************************************************************************************************
#**********************************************************************************************************************************
"""
This is the Python interface to ParaMonte: Plain Powerful Parallel Monte Carlo library.

What is ParaMonte?
==================

ParaMonte is a serial/parallel library of Monte Carlo routines for sampling mathematical 
objective functions of arbitrary-dimensions, in particular, the posterior distributions 
of Bayesian models in data science, Machine Learning, and scientific inference, with the 
design goal of unifying the 

    **automation** (of Monte Carlo simulations), 
    **user-friendliness** (of the library), 
    **accessibility** (from multiple programming environments), 
    **high-performance** (at runtime), and 
    **scalability** (across many parallel processors).  

For more information on the installation, usage, and examples, visit: 

    https://www.cdslab.org/paramonte  

The routines currently supported by the Python interface of ParaMonte include:

    ParaDRAM
    ========
        
        Parallel Delayed-Rejection Adaptive Metropolis-Hastings Markov Chain Monte Carlo Sampler.

        EXAMPLE USAGE:

            import paramonte as pm
            import numpy as np
            def getLogFunc(Point):
                # return the log of the standard multivariate 
                # Normal density function with ndim dimensions
                return -0.5 * np.sum( np.double( Point )**2 )
            pmpd = pm.ParaDRAM()
            pmpd.runSampler ( ndim = 2
                            , getLogFunc = getLogFunc
                            )

"""

import os as _os
import sys as _sys
_sys.path.append(_os.path.dirname(__file__))

import _paramonte as _pm
from paradram import ParaDRAM
from _pmreqs import verify, build

__author__  = "Computational Data Science Lab @ The University of Texas"
__credits__ = "Peter O'Donnell Fellowship"
__version__ = _pm.version

verify(reset=False)