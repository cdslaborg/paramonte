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
####   we ask you to acknowledge the use of the ParaMonte library
####   in your work (education/research/industry/development/...)
####   by citing the ParaMonte library as described on this page:
####
####       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
####
####################################################################################################################################
####################################################################################################################################

import numpy as np

# The number of dimensions of the domain of the objective function.

NDIM = 4

# This is the mean of the MVN distribution.

MEAN =  [0.0,0.0,0.0,0.0]

# This is the covariance matrix of the MVN distribution.

COVMAT =    [ [1.0,0.5,0.5,0.5]
            , [0.5,1.0,0.5,0.5]
            , [0.5,0.5,1.0,0.5]
            , [0.5,0.5,0.5,1.0]
            ]

# This is the inverse of the covariance matrix of the MVN distribution.

INVCOV = np.linalg.inv(COVMAT)

# This is the log of the coefficient used in the definition of the MVN.

MVN_COEF = NDIM * np.log( 1. / np.sqrt(2.*np.pi) ) + np.log( np.sqrt(np.linalg.det(INVCOV)) )

def getLogFunc(point):
    """
    Return the natural logarithm of an NDIM-dimensional Multivariate Normal distribution
    with the mean and covariance matrix as given in the above.
    Reference: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    """
    normedPoint = MEAN - point
    return MVN_COEF - 0.5 * ( np.dot(normedPoint,np.matmul(INVCOV,normedPoint)) )
