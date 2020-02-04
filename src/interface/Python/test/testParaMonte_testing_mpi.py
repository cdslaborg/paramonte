#!/usr/bin/env python
#!C:\ProgramData\Anaconda3\python.exe

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

import os
import numpy as np
import paramonte as pm
import importlib
importlib.reload(pm)
from paramonte.mvn import NDIM, getLogFunc

pd = pm.ParaDRAM()
pd.runSampler( ndim = NDIM
             , getLogFunc = getLogFunc_pntr
             , inputFilePath = os.path.dirname(os.path.abspath(__file__)) + "/input/paramonte.nml"
             , mpiEnabled = True
             , buildMode = "testing"
             )
