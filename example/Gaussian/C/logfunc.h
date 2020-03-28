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
//        - returns the log(probability) (in Neper base) of an ndim-dimensional standard Gaussian probability density (pdf),
//        - that is, with mean zero, and an identity covariance matrix.
//  Input:
//        - ndim: number of dimensions of the domain of the Gaussian function (length of the vector X)
//        - Point: the input real-valued vector of length ndim, at which log(pdf) is computed.
//  Output:
//        - logProb: real scalar number representing the log Gaussian pdf
//  Author:
//        - Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
//
//**********************************************************************************************************************************
//**********************************************************************************************************************************

#ifndef PM_LOG_FUNC
#define PM_LOG_FUNC

#include <math.h>
#include <stdint.h>

#define NDIM 4  // dimension of the domain of the Gaussian distribution function (bivariate Gaussian)

double getLogFunc   (
                    int32_t ,         // ndim
                    double []         // Point
                    );

#endif
