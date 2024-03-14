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
#define NDIM 4

REAL getLogFunc(REAL state[], int32_t ndim){
    //
    //  Return the natural logarithm of an ndim-dimensional Multivariate Normal (MVN)
    //  probability density function (PDF) with the Mean and Covariance Matrix as defined below.
    //  See also: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    //
    const REAL LOG_INVERSE_SQRT_TWO_PI = log(0.398942280401432);  // log(1/sqrt(2*Pi))
    const REAL MEAN[NDIM] = {-6., -2., +2., +6.}; // mean vector of the MVN.
    const REAL COVMAT[NDIM][NDIM] = {
                                    {+1.0, +0.5, +0.5, +0.5},
                                    {+0.5, +1.0, +0.5, +0.5},
                                    {+0.5, +0.5, +1.0, +0.5},
                                    {+0.5, +0.5, +0.5, +1.0}
                                    }; // covariance matrix of the MVN.
    const REAL INVCOV[NDIM][NDIM] = {
                                    {+1.6, -0.4, -0.4, -0.4},
                                    {-0.4, +1.6, -0.4, -0.4},
                                    {-0.4, -0.4, +1.6, -0.4},
                                    {-0.4, -0.4, -0.4, +1.6}
                                    }; // inverse covariance matrix of the MVN.
    const REAL LOG_SQRT_DET_INVCOV = 0.581575404902840; // logarithm of square root of the determinant of the inverse covariance matrix.

    // subtract mean vector from the input point.

    int idim;
    REAL stateNormed[NDIM];
    for(idim = 0; idim < NDIM; idim++){
        stateNormed[idim] = state[idim] - MEAN[idim];
    }

    // compute the log probability density function of the MVN

    REAL exponentTerm = 0.;
    for(idim = 0; idim < NDIM; idim++){
        REAL matMulResult[NDIM];
        matMulResult[idim] = 0.;
        int jdim;
        for(jdim = 0; jdim < NDIM; jdim++){
            matMulResult[idim] += INVCOV[idim][jdim] * stateNormed[jdim];
        }
        exponentTerm += stateNormed[idim] * matMulResult[idim];
    }
    REAL logFunc = NDIM * LOG_INVERSE_SQRT_TWO_PI + LOG_SQRT_DET_INVCOV - 0.5 * exponentTerm;
    ///for(int i = 0; i < NDIM; i++) {
    ///    printf("%Le", state[i]);
    ///}
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