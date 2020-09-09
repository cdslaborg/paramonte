!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   MIT License
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   Permission is hereby granted, free of charge, to any person obtaining a 
!   copy of this software and associated documentation files (the "Software"), 
!   to deal in the Software without restriction, including without limitation 
!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!   and/or sell copies of the Software, and to permit persons to whom the 
!   Software is furnished to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be 
!   included in all copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!   ACKNOWLEDGMENT
!
!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!   As per the ParaMonte library license agreement terms, if you use any parts of 
!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!   work (education/research/industry/development/...) by citing the ParaMonte 
!   library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module LogFunc_mod

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%  
    !%  Description:
    !%      -   Return the natural logarithm of an ndim-dimensional Multivariate Normal (MVN) 
    !%          probability density function (PDF) with the Mean and Covariance Matrix as defined below.
    !%          Reference: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    !%  Input:
    !%      -   ndim:   The number of dimensions of the domain of the objective function.
    !%      -   point:  The input 64-bit real-valued vector of length ndim, 
    !%                  at which the natural logarithm of objective function is computed.
    !%  Output:
    !%      - logFunc:  A 64-bit real scalar number representing the natural logarithm of the objective function.
    !%  Author:
    !%      - Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
    !%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    use paramonte, only: IK, RK

    implicit none

    ! The number of dimensions of the domain of the objective function

    integer(IK), parameter  :: NDIM = 4_IK

    ! The mean vector of the MVN

    real(RK), parameter     :: MEAN(NDIM) = [0._RK, 0._RK, 0._RK, 0._RK]

    ! The coefficient of the Multivariate Normal Distribution: log(1/sqrt(2*Pi)^ndim)

    real(RK), parameter     :: LOG_INVERSE_SQRT_TWO_PI = log(1._RK/sqrt(2._RK*acos(-1._RK)))

    ! The covariance matrix of the MVN

    real(RK), parameter     :: COVMAT(NDIM,NDIM) =  reshape([ 1.0_RK,0.5_RK,0.5_RK,0.5_RK &
                                                            , 0.5_RK,1.0_RK,0.5_RK,0.5_RK &
                                                            , 0.5_RK,0.5_RK,1.0_RK,0.5_RK &
                                                            , 0.5_RK,0.5_RK,0.5_RK,1.0_RK ], shape=shape(COVMAT))

    ! The inverse covariance matrix of the MVN

    real(RK), parameter     :: INVCOVMAT(NDIM,NDIM) = reshape(  [ +1.6_RK, -0.4_RK, -0.4_RK, -0.4_RK &
                                                                , -0.4_RK, +1.6_RK, -0.4_RK, -0.4_RK &
                                                                , -0.4_RK, -0.4_RK, +1.6_RK, -0.4_RK &
                                                                , -0.4_RK, -0.4_RK, -0.4_RK, +1.6_RK ], shape=shape(INVCOVMAT))

    ! The logarithm of square root of the determinant of the inverse covariance matrix of the Multivariate Normal Distribution 

    real(RK), parameter     :: LOG_SQRT_DET_INV_COV = 0.581575404902840_RK

    ! MVN_COEF

    real(RK), parameter     :: MVN_COEF = NDIM * LOG_INVERSE_SQRT_TWO_PI + LOG_SQRT_DET_INV_COV

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogFunc(ndim,Point) result(logFunc)
        ! Return the negative natural logarithm of MVN distribution evaluated at the input vector point.
        implicit none
        integer(IK), intent(in) :: ndim
        real(RK), intent(in)    :: Point(ndim)
        real(RK)                :: NormedPoint(ndim)
        real(RK)                :: logFunc
        NormedPoint = Point - MEAN
        logFunc = MVN_COEF - 0.5_RK * ( dot_product(NormedPoint,matmul(INVCOVMAT,NormedPoint)) )
    end function getLogFunc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module LogFunc_mod