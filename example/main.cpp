////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  
//  Description:
//      +   Run the Monte Carlo sampler of the ParaMonte library given the input log-target density function `getLogFunc()`.
//  Output:
//      +   The simulation output files.
//  Author:
//      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
//  Visit:
//      +   https://www.cdslab.org/paramonte
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cstring>
#include "logfunc.hpp"      // getLogFunc, NDIM
#include "paramonte.hpp"    // runParaDRAM

int main()
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////