!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%  
!%  Description:
!%      +   Return the natural logarithm of an ndim-dimensional Multivariate Normal (MVN) 
!%          probability density function (PDF) with the Mean and Covariance Matrix as defined below.
!%          Reference: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
!%  Input:
!%      +   ndim:       The number of dimensions of the domain of the objective function.
!%      +   point:      The input 64-bit real-valued vector of length ndim, 
!%                      at which the natural logarithm of objective function is computed.
!%  Output:
!%      +   logFunc:    A 64-bit real scalar number representing the natural logarithm of the objective function.
!%  Author:
!%      +   Computational Data Science Lab, Monday 9:03 AM, May 16 2016, ICES, UT Austin
!%  Visit:
!%      +   https://www.cdslab.org/paramonte
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module LogFunc_mod

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

contains

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

end module LogFunc_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
