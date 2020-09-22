////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  
//  Description:
//      +   Return the natural logarithm of an ndim-dimensional Multivariate Normal (MVN) 
//          probability density function (PDF) with the Mean and Covariance Matrix as defined below.
//          Reference: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
//  Input:
//      +   ndim:       The number of dimensions of the domain of the objective function.
//      +   point:      The input 64-bit real-valued vector of length ndim, 
//                      at which the natural logarithm of objective function is computed.
//  Output:
//      +   logFunc:    A 64-bit real scalar number representing the natural logarithm of the objective function.
//  Author:
//      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
//  Visit:
//      +   https://www.cdslab.org/paramonte
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "logfunc.hpp"

double getLogFunc   (
                    int32_t ndim, 
                    double Point[]
                    )
{

    // multivariate normal (MVN) distribution specifications: 

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

    int i;
    double NormedPoint[NDIM];
	for(i = 0; i < NDIM; i++){
        NormedPoint[i] = Point[i] - MEAN[i];
    }

    // compute the log probability density function of the MVN

    double exponentTerm = 0.;
	for(i = 0; i < NDIM; i++){

        double MatMulResult[NDIM];
        MatMulResult[i] = 0.;
        int j;
        for(j = 0; j < NDIM; j++){
            MatMulResult[i] += INVCOVMAT[i][j] * NormedPoint[j];
        }

        exponentTerm += NormedPoint[i] * MatMulResult[i];

    }

    double logFunc = NDIM * LOG_INVERSE_SQRT_TWO_PI + LOG_SQRT_DET_INV_COV - 0.5*exponentTerm;

    return logFunc;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////