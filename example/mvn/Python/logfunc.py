####################################################################################################################################
##  
##  Description:
##      +   Return the natural logarithm of an ndim-dimensional Multivariate Normal (MVN) 
##          probability density function (PDF) with the Mean and Covariance Matrix as defined below.
##          Reference: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
##  Input:
##      +   point:      The input 64-bit real-valued vector of length ndim, 
##                      at which the natural logarithm of objective function is computed.
##  Output:
##      +   logFunc:    A 64-bit real scalar number representing the natural logarithm of the objective function.
##  Author:
##      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
##  Visit:
##      +   https://www.cdslab.org/paramonte
##
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
