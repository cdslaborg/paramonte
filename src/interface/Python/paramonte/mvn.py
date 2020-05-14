####################################################################################################################################
####################################################################################################################################
####
####   ParaMonte: plain powerful parallel Monte Carlo library.
####
####   Copyright (C) 2012-present, The Computational Data Science Lab
####
####   This file is part of the ParaMonte library.
####
####   ParaMonte is free software: you can redistribute it and/or modify it
####   under the terms of the GNU Lesser General Public License as published
####   by the Free Software Foundation, version 3 of the License.
####
####   ParaMonte is distributed in the hope that it will be useful,
####   but WITHOUT ANY WARRANTY; without even the implied warranty of
####   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
####   GNU Lesser General Public License for more details.
####
####   You should have received a copy of the GNU Lesser General Public License
####   along with the ParaMonte library. If not, see,
####
####       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
####
####   ACKNOWLEDGMENT
####
####   As per the ParaMonte library license agreement terms,
####   if you use any parts of this library for any purposes,
####   we ask you to acknowledge the ParaMonte library's usage
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import numpy as np

# number of dimensions of the distribution

NDIM = 1

# the coefficient of the Standard Multivariate Normal Distribution: log(1/sqrt(2*Pi)^ndim)

LOG_SMVN_COEF = NDIM * np.log( 1.0 / np.sqrt(2.0*np.pi) )

# define Python function

#def getLogFunc_pntr(ndim,Point): return LOG_SMVN_COEF - 0.5 * np.sum( np.double( Point[0:ndim[0]] )**2 )

def getLogFunc_pntr(ndim,Point): return LOG_SMVN_COEF - 0.5 * np.sum( np.double( Point[0:ndim] )**2 )

def getLogFunc(Point): return LOG_SMVN_COEF - 0.5 * np.sum( np.double( Point )**2 )
