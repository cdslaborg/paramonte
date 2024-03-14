//
//  Description
//  -----------
//
//      Run the Monte Carlo samplers of the ParaMonte library
//      given the input log-target density function `getLogFunc()`.
//
//  Author
//  ------
//
//      Computational Data Science Lab, Monday 9:03 AM, May 16, 2016,
//      Institute for Computational Engineering and Sciences (ICES),
//      The University of Texas at Austin
//
//  Documentation
//  -------------
//
//      https://www.cdslab.org/paramonte
//
//  LICENSE
//  -------
//
//      https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
//
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "pm_sampling.h"
// The `runParaDRAM` and `REAL` macros can be set to other possibilities.
#define runParaDRAM runParaDRAMD
#define REAL double
#define NDIM 2

REAL getLogFunc(REAL state[], int32_t ndim){
    //
    //  Return the negative natural logarithm of Himmelblau function at the specified state.
    //  See also: https://en.wikipedia.org/wiki/Himmelblau%27s_function
    //
    REAL logFunc = -log ( pow(pow(state[0], 2) + state[1] - 11, 2)
                        + pow(state[0] + pow(state[1], 2) -  7, 2)
                        + 0.1
                        );
    return logFunc;
}

int main()
{
    int32_t failed;
    const char* input = "./input.nml"; // null-terminated string. It can be empty or NULL.

    // Call the ParaDRAM Adaptive MCMC sampler with the requested floating point precision.

    failed = runParaDRAM(&getLogFunc, NDIM, input);
    if (failed != 0) exit(failed);

    // Call the ParaDRAM Adaptive MCMC sampler with the requested floating point precision with an internal namelist input specifications.

    failed = runParaDRAM(&getLogFunc, NDIM, "&paradram parallelismNumThread = 16, outputChainSize = 30000, parallelismMpiFinalizeEnabled = false /");
    if (failed != 0) exit(failed);

    // Call the ParaDRAM Adaptive MCMC sampler with the requested floating point precision without input specifications.

    failed = runParaDRAM(&getLogFunc, NDIM, NULL);
    if (failed != 0) exit(failed);
}