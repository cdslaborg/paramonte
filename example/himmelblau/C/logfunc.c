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

#include "logfunc.h"

double getLogFunc   (
                    int32_t ndim, 
                    double Point[]
                    )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  
    //  Description:
    //      -   Return the negative natural logarithm of Himmelblau's function.
    //          Reference: https://en.wikipedia.org/wiki/Himmelblau%27s_function
    //  Input:
    //      -   ndim:   The number of dimensions of the domain of the objective function.
    //      -   point:  The input 64-bit real-valued vector of length ndim, 
    //                  at which the natural logarithm of objective function is computed.
    //  Output:
    //      - logFunc:  A 64-bit real scalar number representing the natural logarithm of the objective function.
    //  Author:
    //      - Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double logFunc = -log   ( pow( pow(point[0], 2) + point[1] - 11, 2 ) 
                            + pow( point[0] + pow(point[1], 2) - 7 , 2 )
                            )

    return logFunc;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
