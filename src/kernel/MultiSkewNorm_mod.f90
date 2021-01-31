!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief This module contains the classes and procedures for the Multivariate Skew Normal distribution of
!> Azzalini, 2009, Statistical applications of the multivariate skew-normal distribution.
!> \author Amir Shahmoradi

module MultiSkewNorm_mod

    use Constants_mod, only: RK, IK
    use Err_mod, only: Err_type

    implicit none

    character(len=*), parameter :: MODULE_NAME = "@MultiSkewNorm_mod"

    interface getLogProbMSN
        module procedure :: getLogProbMSN_RK
    end interface getLogProbMSN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Return the log-probability-density-function value of obtaining `Point` from a Multivariate-Skew-Normal distribution
    !> characterize by the input `InvCovMat` and shape vector `Alpha`.
    !>
    !> \param[in]   nd                  : The number of dimensions of the Multivariate-Skew-Normal distribution.
    !> \param[in]   MeanVec             : The Mean of the Multivariate-Skew-Normal distribution.
    !> \param[in]   InvCovMat           : The inverse covariance matrix of the Multivariate-Skew-Normal distribution.
    !> \param[in]   logSqrtDetInvCovMat : The natural logarithm of the square-root of the determinant of the `InvCovMat`.
    !> \param[in]   Alpha               : The shape parameter of the Multivariate-Skew-Normal distribution.
    !> \param[in]   Point               : The point at which the log-probability density must be computed.
    !>
    !> \return
    !> `logProbSkewNorm`    :   the log-probability-density-function value of obtaining `Point`.
    !>
    !> \remark
    !> A positive value for any element of `Alpha` causes a negative skeness along the corresponding dimension.
    !> This looks as if the distribution is truncated or cut toward negative values.
    !> Amir Shahmoradi, January 30, 2021, 7:49 AM, Dallas, TX
    pure function getLogProbMSN_RK(nd,MeanVec,InvCovMat,logSqrtDetInvCovMat,Alpha,Point) result(logProbSkewNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMSN_RK
#endif
        use Constants_mod, only: LOGINVSQRT2PI, INVSQRT2
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)
        real(RK)   , intent(in) :: logSqrtDetInvCovMat
        real(RK)   , intent(in) :: Alpha(nd)
        real(RK)   , intent(in) :: Point(nd)
        real(RK)                :: NormedPoint(nd)
        real(RK)                :: logProbSkewNorm
        real(RK)    , parameter :: LOG2 = log(2._RK)
        real(RK)    , parameter :: LOGHALF = log(0.5_RK)
        NormedPoint = Point - MeanVec
        logProbSkewNorm = LOG2 + nd*LOGINVSQRT2PI + logSqrtDetInvCovMat & ! LCOV_EXCL_LINE
                        - 0.5_RK * dot_product(NormedPoint, matmul(InvCovMat,NormedPoint)) & ! LCOV_EXCL_LINE
                        + LOGHALF + log( 1._RK + erf(INVSQRT2 * dot_product(NormedPoint,Alpha)) )
    end function getLogProbMSN_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a MultiVariate Skew Normal (MSN) random vector with the given mean and
    !> covariance matrix represented by the the input Cholesky factorization.
    !>
    !> \param[in]   nd              :   The number of dimensions of the MSN distribution.
    !> \param[in]   MeanVec         :   The mean vector of the MSN distribution.
    !> \param[in]   AugChoLow       :   The Cholesky lower triangle of the augmented covariance matrix of the MSN distribution.
    !> \param[in]   AugChoDia       :   The Diagonal elements of the Cholesky lower triangle of the augmented covariance matrix of the MSN distribution.
    !>
    !> \return
    !> `RandMSN` : The randomly generated MSN vector.
    !>
    !> \author
    !> Amir Shahmoradi, January 30, 2021, 5:36 AM, Dallas, TX
    function getRandMSN(nd,MeanVec,AugChoLow,AugChoDia) result(RandMSN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandMSN
#endif
        use Statistics_mod, only: getRandGaus, getRandMVN
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: AugChoLow(0:nd,0:nd), AugChoDia(0:nd)
        real(RK)                :: RandMSN(nd), normrnd
        logical                 :: switchDisabled
        integer(IK)             :: j, i
        normrnd = getRandGaus()
        switchDisabled = normrnd > 0._RK
        RandMSN(1:nd) = AugChoLow(1:nd,0) * normrnd
        do j = 1, nd
            normrnd = getRandGaus()
            RandMSN(j) = RandMSN(j) + AugChoDia(j) * normrnd
            do i = j+1, nd
                RandMSN(i) = RandMSN(i) + AugChoLow(i,j) * normrnd
            end do
        end do
        if (switchDisabled) then
            RandMSN = MeanVec(1:nd) + RandMSN
        else
            RandMSN = MeanVec(1:nd) - RandMSN
        end if
    end function getRandMSN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the augmented Cholesky Factorization of the input Covariance Matrix with the shape vector `Alpha`.
    !> The output of this procedure is to be passed to [getRandMSN()](@ref getrandmsn) for MultiVariate Skew Normal
    !> random number generation.
    !>
    !> \param[in]   nd              :   The number of dimensions of the MSN distribution.
    !> \param[in]   Alpha           :   The input shape vector of the target MSN distribution.
    !> \param[in]   CovMat          :   The FULL input covariance matrix.
    !> \param[out]  AugChoLow       :   The Cholesky Lower triangle of the augmented covariance matrix of the MSN distribution.
    !> \param[out]  AugChoDia       :   The Diagonal elements of the Cholesky lower triangle of the augmented covariance matrix of the MSN distribution.
    !>
    !> \remark
    !> If the input `CovMat` is only upper, symmetrize `CovMat` via [symmetrizeUpperSquareMatrix()](@ref matrix_mod::symmetrizeuppersquarematrix).
    !>
    !> \author
    !> Amir Shahmoradi, Jan 30, 2021, 2:19 AM, Dallas, TX
    subroutine getAugChoFac(nd,Alpha,CovMat,AugChoLow,AugChoDia)
        use Matrix_mod, only: getCholeskyFactor
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: nd
        real(RK)    , intent(in)    :: Alpha(nd)
        real(RK)    , intent(in)    :: CovMat(nd,nd)
        real(RK)    , intent(out)   :: AugChoLow(0:nd,0:nd), AugChoDia(0:nd)
        real(RK)                    :: Delta(nd), coef
        integer(IK)                 :: j
        Delta(1:nd) = matmul(CovMat, Alpha)
        coef = 1._RK / sqrt(1._RK + dot_product(Alpha,Delta))
        Delta = Delta * coef
        AugChoLow(0,0) = 1._RK
        do j = 1, nd
            AugChoLow(0,j) = Delta(j)
            AugChoLow(1:j,j) = CovMat(1:j,j)
        end do
        call getCholeskyFactor(nd+1,AugChoLow,AugChoDia)
    end subroutine getAugChoFac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module MultiSkewNorm_mod ! LCOV_EXCL_LINE