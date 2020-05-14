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
//   we ask you to acknowledge the ParaMonte library's usage
//   in your work (education/research/industry/development/...)
//   by citing the ParaMonte library as described on this page:
//
//       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  
//  Description:
//        - The C header file containing the prototypes of the ParaMonte library routines to be called from C/C++.
//  Prototypes:
//        - runParaDRAM()
//  Author:
//        - Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ParaMonte
#define ParaMonte

#include <stdint.h>

#if defined DLL_ENABLED
extern __declspec(dllimport) 
#endif
void runParaDRAM(
                // ndim: dimension of the domain of the LogFunc
                int32_t , 
                // getLogFunc(ndim, Point(ndim)): procedure pointer to the LogFunc
                double (*)  (
                            int32_t , 
                            double []
                            ), 
                // inputFile: ParaDRAM input file, containing the path the file containing a list of
                // all optional input variables and values
                char [], 
                // inputFileLen: the length of the inputFile char vector: int32_t inputFileLen = strlen(inputFile);
                int32_t
                );

#endif
