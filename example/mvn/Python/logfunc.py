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
from scipy.stats import multivariate_normal

NDIM = 4 # number of dimensions of the distribution

mvn = multivariate_normal   ( mean =  [0.0,0.0,0.0,0.0]
                            , cov = [ [1.0,0.5,0.5,0.5]
                                    , [0.5,1.0,0.5,0.5]
                                    , [0.5,0.5,1.0,0.5]
                                    , [0.5,0.5,0.5,1.0]
                                    ]
                            )

def getLogFunc(point): return np.log(mvn.pdf(point))
