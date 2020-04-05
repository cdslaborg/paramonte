//**********************************************************************************************************************************
//**********************************************************************************************************************************
//
//  ParaMonte: plain powerful parallel Monte Carlo library.
//
//  Copyright (C) 2012-present, The Computational Data Science Lab
//
//  This file is part of ParaMonte library. 
//
//  ParaMonte is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, version 3 of the License.
//
//  ParaMonte is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
//
//**********************************************************************************************************************************
//**********************************************************************************************************************************

//**********************************************************************************************************************************
//**********************************************************************************************************************************
//  
//  Description:
//        - returns the log(PDF) (in Neper base) of an ndim-dimensional Multivariate Normal (MVN) probability density function (PDF)
//  Input:
//        - ndim: number of dimensions of the domain of the MVN function (length of the vector X)
//        - Point: the input real-valued vector of length ndim, at which log(PDF) is computed.
//  Output:
//        - logFunc: real scalar number representing the log MVN PDF
//  Author:
//        - Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
//
//**********************************************************************************************************************************
//**********************************************************************************************************************************

#include "logfunc.h"

double getLogFunc   (
                    int32_t ndim, 
                    double Point[]
                    )
{

    // multivariate normal (MVN) distribution specifications: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

    const double LOG_INVERSE_SQRT_TWO_PI = log(0.398942280401432);  // log(1/sqrt(2*Pi))
    const double MEAN[NDIM] = {0., 0., 0., 0.};                     // mean vector of the MVN
    const double COVMAT[NDIM][NDIM] =       {
                                            {1.0,0.5,0.5,0.5},
                                            {0.5,1.0,0.5,0.5},
                                            {0.5,0.5,1.0,0.5},
                                            {0.5,0.5,0.5,1.0}
                                            };                      // covariance matrix of the MVN
    const double INVCOVMAT[NDIM][NDIM] =    {
                                            {+1.6, -0.4, -0.4, -0.4},
                                            {-0.4, +1.6, -0.4, -0.4},
                                            {-0.4, -0.4, +1.6, -0.4},
                                            {-0.4, -0.4, -0.4, +1.6}
                                            };                      // inverse covariance matrix of the MVN
    const double LOG_SQRT_DET_INV_COV = 0.581575404902840;          // logarithm of square root of the determinant of the inverse covariance matrix

    // subtract mean vector from the input point

    double NormedPoint[NDIM];
	for(int i = 0; i < NDIM; i++){
        NormedPoint[i] = Point[i] - MEAN[i];
    }

    // compute the log probability density function of the MVN

    double exponentTerm = 0.;
	for(int i = 0; i < NDIM; i++){

        double MatMulResult[NDIM];
        MatMulResult[i] = 0.;
        for(int j = 0; j < NDIM; j++){
            MatMulResult[i] += INVCOVMAT[i][j] * NormedPoint[j];
        }

        exponentTerm += NormedPoint[i] * MatMulResult[i];

    }

    double logFunc = NDIM * LOG_INVERSE_SQRT_TWO_PI + LOG_SQRT_DET_INV_COV - 0.5*exponentTerm;

    return logFunc;

}

//**********************************************************************************************************************************
//**********************************************************************************************************************************