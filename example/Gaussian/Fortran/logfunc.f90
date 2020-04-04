!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module LogFunc_mod

    ! MultiVariate Normal (MVN) specifications: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

    !use, intrinsic :: iso_fortran_env, only: RK => real64, IK => int32
    use paramonte, only: IK, RK

    implicit none

    ! number of dimensions of the distribution

    integer(IK), parameter  :: NDIM = 4_IK

    ! mean vector of the MVN

    real(RK), parameter     :: MEAN(NDIM) = [0._RK, 0._RK, 0._RK, 0._RK]

    ! Multivariate Normal Distribution coefficient: log(1/sqrt(2*Pi)^ndim)

    real(RK), parameter     :: LOG_INVERSE_SQRT_TWO_PI = log(1._RK/sqrt(2._RK*acos(-1._RK)))

    ! covariance matrix of the MVN

    real(RK), parameter     :: COVMAT(NDIM,NDIM) =  reshape([ 1.0_RK,0.5_RK,0.5_RK,0.5_RK &
                                                            , 0.5_RK,1.0_RK,0.5_RK,0.5_RK &
                                                            , 0.5_RK,0.5_RK,1.0_RK,0.5_RK &
                                                            , 0.5_RK,0.5_RK,0.5_RK,1.0_RK ], shape=shape(COVMAT))

    ! inverse covariance matrix of the MVN

    real(RK), parameter     :: INVCOVMAT(NDIM,NDIM) = reshape(  [ +1.6_RK, -0.4_RK, -0.4_RK, -0.4_RK &
                                                                , -0.4_RK, +1.6_RK, -0.4_RK, -0.4_RK &
                                                                , -0.4_RK, -0.4_RK, +1.6_RK, -0.4_RK &
                                                                , -0.4_RK, -0.4_RK, -0.4_RK, +1.6_RK ], shape=shape(INVCOVMAT))

    ! logarithm of square root of the determinant of the inverse covariance matrix of the Multivariate Normal Distribution 

    real(RK), parameter     :: LOG_SQRT_DET_INV_COV = 0.581575404902840_RK

    ! MVN_COEF

    real(RK), parameter     :: MVN_COEF = NDIM * LOG_INVERSE_SQRT_TWO_PI + LOG_SQRT_DET_INV_COV

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function getLogFunc(ndim,Point) result(logFunc)
        ! This function returns the probability density function of the standard multivariate normal distribution of ndim dimensions.
        implicit none
        integer(IK), intent(in) :: ndim
        real(RK), intent(in)    :: Point(ndim)
        real(RK)                :: NormedPoint(ndim)
        real(RK)                :: logFunc
        NormedPoint = Point - MEAN
        logFunc = MVN_COEF - 0.5_RK * ( dot_product(NormedPoint,matmul(INVCOVMAT,NormedPoint)) )
    end function getLogFunc

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module LogFunc_mod