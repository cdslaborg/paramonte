////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   This file is part of ParaMonte: Parallel Monte Carlo and Machine Learning library.
//
//   Copyright (C) 2012-present, The Computational Data Science Lab
//
//   https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdint.h>
#include <string.h>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mex.h"
#include "matrix.h"

int32_t getImageCountMPI(void);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // prhs: method, getLogFunc, ndim, input.
    if(nlhs != 1) mexErrMsgIdAndTxt("Mex:ParaMonte:maxlhs", "Internal ParaMonte MATLAB library error occurred: There must be exactly one output arguments.");
    if (nrhs != 1) {mexErrMsgIdAndTxt("Mex:ParaMonte:invalidNumInputs", "Internal ParaMonte MATLAB library error occurred: There must be exactly one input arguments.");}
    mwSize iarg;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    iarg = 0;
    if (mxGetM(prhs[iarg]) != 1) mexErrMsgIdAndTxt("Mex:ParaMonte:inputNotVector", "Input must be a row vector of characters.");
    if (mxIsChar(prhs[iarg]) != 1) mexErrMsgIdAndTxt("Mex:ParaMonte:inputNotString", "Internal ParaMonte MATLAB library error occurred: Input #0 (method) must be a string.");

    char *method;
    method = mxArrayToString(prhs[iarg]);
    if(method == NULL) mexErrMsgIdAndTxt("Mex:ParaMonte:conversionFailed", "Internal ParaMonte MATLAB library error occurred: Could not convert input #0 (method) to string.");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (strncmp(method, "getImageCountMPI", 16) == 0) {
        // Create matrix for the return argument.
        plhs[0] = mxCreateDoubleScalar((double) getImageCountMPI());
    } else {
        mexErrMsgIdAndTxt("Mex:ParaMonte:invalidMethod", "Internal ParaMonte MATLAB library error occurred: Invalid input function.");
    }
    mxFree(method);
    //return;
    exit(0);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}