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

module Test_ParaMonte_mod

    !use Constants_mod, only: IK, RK
    use ParaMonte, only: MODULE_NAME, ParaDRAM, IK, RK
    use Test_mod, only: Test_type

    implicit none

    ! Standard MultiVariate Normal (SMVN) specifications: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    integer(IK) , parameter :: NDIM = 1_IK                                              ! number of dimensions of the distribution
    real(RK)    , parameter :: LOG_SMVN_COEF = NDIM*log(1._RK/sqrt(2._RK*acos(-1._RK))) ! log(1/sqrt(2*Pi)^ndim)

    private
    public :: test_ParaDRAM

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_ParaDRAM()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_runParaDRAM()
        call Test%finalize()

    end subroutine test_ParaDRAM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_runParaDRAM()

        implicit none
        character(:), allocatable   :: internalFile
        type(ParaDRAM)              :: PD

        internalFile = "&ParaDRAM randomSeed = 1111 chainSize = 30000 /"

        call Test%testing("ParaDRAM class")

        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFunc &
                            , inputFile = Test%inDir//"paramonte.nml" &
                            !, inputFile = internalFile &
                            !, inputFile = " " &
                            )

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_runParaDRAM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

#if defined CFI_ENABLED
    function getLogFunc(ndim,Point) result(logFunc) bind(C)
#else
    function getLogFunc(ndim,Point) result(logFunc)
#endif
        ! This function returns the probability density function of the standard multivariate normal distribution of ndim dimensions.
        use Statistics_mod, only: getLogProbGausMix
        use Constants_mod, only : IK, RK
        implicit none
#if defined CFI_ENABLED
        integer(IK), intent(in), value  :: ndim
#else
        integer(IK), intent(in)         :: ndim
#endif
        real(RK), intent(in)            :: Point(ndim)
        real(RK)                        :: logFunc
        !block
        !    use System_mod, only: sleep
        !    use Err_mod, only: Err_type
        !    type(Err_type) :: Err
        !    call sleep(seconds=5000.e-6_RK,Err=Err)
        !end block
        block
            real(RK), allocatable :: unifrnd(:,:)
            allocate(unifrnd(200,20))
            call random_number(unifrnd)
            logFunc = sum(unifrnd) + LOG_SMVN_COEF - 0.5_RK * sum(Point**2) - sum(unifrnd)
            deallocate(unifrnd)
        end block
        logFunc = LOG_SMVN_COEF - 0.5_RK * sum(Point**2)
       !block
       !    integer(IK), parameter :: nmode = 2_IK
       !    real(RK) :: LogAmplitude(nmode), MeanVec(nmode), InvCovMat(nmode), LogSqrtDetInvCovMat(nmode)
       !    LogAmplitude        = [1._RK, 1._RK]
       !    MeanVec             = [0._RK, 7._RK]
       !    InvCovMat           = [1._RK,1._RK]
       !    LogSqrtDetInvCovMat = [1._RK,1._RK]
       !    logFunc = getLogProbGausMix ( nmode = 2_IK &
       !                                , nd = 1_IK &
       !                                , np = 1_IK &
       !                                , LogAmplitude = LogAmplitude &
       !                                , MeanVec = MeanVec &
       !                                , InvCovMat = InvCovMat &
       !                                , LogSqrtDetInvCovMat = LogSqrtDetInvCovMat &
       !                                , Point = Point(1) &
       !                                )
       !end block
    end function getLogFunc

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_ParaMonte_mod