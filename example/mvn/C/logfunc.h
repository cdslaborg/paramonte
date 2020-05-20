////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   ParaMonte: plain powerful parallel Monte Carlo library.
//
//   Copyright (C) 2012-present, The Computational Data Science Lab
//
//   This file is part of the ParaMonte library.
//
//   ParaMonte is free software: you can redistribute it and/or modify it
//   under the terms of the GNU Lesser General Public License as published
//   by the Free Software Foundation, version 3 of the License.
//
//   ParaMonte is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//   GNU Lesser General Public License for more details.
//
//   You should have received a copy of the GNU Lesser General Public License
//   along with the ParaMonte library. If not, see,
//
//       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
//
//   ACKNOWLEDGMENT
//
//   As per the ParaMonte library license agreement terms,
//   if you use any parts of this library for any purposes,
//   we ask you to acknowledge the use of the ParaMonte library
//   in your work (education/research/industry/development/...)
//   by citing the ParaMonte library as described on this page:
//
//       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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