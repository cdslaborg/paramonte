!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module Test_ParaMonte_mod

    use Constants_mod, only: IK, RK
    use Test_mod, only: Test_type

    implicit none

    ! Standard MultiVariate Normal (SMVN) specifications: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    integer(IK) , parameter :: NDIM = 1_IK                                              ! number of dimensions of the distribution
    real(RK)    , parameter :: LOG_SMVN_COEF = NDIM*log(1._RK/sqrt(2._RK*acos(-1._RK))) ! log(1/sqrt(2*Pi)^ndim)

    private
    public :: test_ParaMonte

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_ParaMonte()
        implicit none
        call test_ParaDRAM()
        call test_ParaDISE()
    end subroutine test_ParaMonte

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_ParaDRAM()

        use ParaDRAM_mod, only: MODULE_NAME
        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_runParaDRAM()
        call Test%finalize()

    end subroutine test_ParaDRAM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_ParaDISE()

        use ParaDISE_mod, only: MODULE_NAME
        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_runParaDISE()
        call Test%finalize()

    end subroutine test_ParaDISE

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_runParaDRAM()

!use Statistics_mod, only: paradramPrintEnabled
        use ParaDRAM_mod, only: ParaDRAM_type
        implicit none
        character(:), allocatable   :: internalFile
        type(ParaDRAM_type)         :: PD

        internalFile = "&ParaDRAM randomSeed = 1111 chainSize = 30000 /"

        call Test%testing("ParaDRAM class")

!paradramPrintEnabled = .true.
        call PD%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFunc &
                            , inputFile = Test%inDir//"paramonte.nml" &
                            , mpiFinalizeRequested = .false. &
                            !, inputFile = internalFile &
                            !, inputFile = " " &
                            )

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_runParaDRAM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_runParaDISE()

!use Statistics_mod, only: paradisePrintEnabled
        use ParaDISE_mod, only: ParaDISE_type
        implicit none
        character(:), allocatable   :: internalFile
        type(ParaDISE_type)         :: PS

        internalFile = "&ParaDISE randomSeed = 1111 chainSize = 30000 /"

        call Test%testing("ParaDISE class")

!paradisePrintEnabled = .true.
        call PS%runSampler  ( ndim = NDIM &
                            , getLogFunc = getLogFunc &
                            , inputFile = Test%inDir//"paramonte.nml" &
                            , mpiFinalizeRequested = .false. &
                            !, inputFile = internalFile &
                            !, inputFile = " " &
                            )

        !Test%assertion = .true.
        !call Test%verify()
        call Test%skipping()

    end subroutine test_runParaDISE

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