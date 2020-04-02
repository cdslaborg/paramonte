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

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "logfunc.h"    // getLogFunc, NDIM
#include "paramonte.h"  // runParaDRAM

int main(int argc, char *argv[])
{
    char inputFile[] = "./paramonte.in";
    int32_t inputFileLen = strlen(inputFile);
    int32_t ndim = NDIM;

    // C rules for argument passing apply here
    runParaDRAM ( ndim
                , &getLogFunc
                , inputFile
                , inputFileLen
                );
    return 0;
}
