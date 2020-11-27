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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [StarFormation_mod](@ref starformation_mod).
!>  @author Amir Shahmoradi

module Test_StarFormation_mod

    use StarFormation_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_StarFormation

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_StarFormation()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getLogRateH06, "test_getLogRateH06")
        call Test%run(test_getLogRateL08, "test_getLogRateL08")
        call Test%run(test_getLogRateB10, "test_getLogRateB10")
        call Test%run(test_getLogRateM14, "test_getLogRateM14")
        call Test%run(test_getLogRateP15, "test_getLogRateP15")
        call Test%run(test_getLogRateM17, "test_getLogRateM17")
        call Test%run(test_getLogRateF18, "test_getLogRateF18")
#if defined CODECOV_ENABLED
        call Test%run(test_getBinaryMergerRate_1, "test_getBinaryMergerRate_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getBinaryMergerRate_2, "test_getBinaryMergerRate_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
#endif
        call Test%run(test_getBinaryMergerRateS15_1, "test_getBinaryMergerRateS15_1")
        call Test%run(test_getBinaryMergerRateS15_2, "test_getBinaryMergerRateS15_2")
        call Test%run(test_getBinaryMergerRateS15_3, "test_getBinaryMergerRateS15_3")
        call Test%run(test_getBinaryMergerRateS15_4, "test_getBinaryMergerRateS15_4")
        call Test%run(test_getLogBinaryMergerRateLognormH06_1, "test_getLogBinaryMergerRateLognormH06_1")
        call Test%run(test_getLogBinaryMergerRateLognormL08_1, "test_getLogBinaryMergerRateLognormL08_1")
        call Test%run(test_getLogBinaryMergerRateLognormB10_1, "test_getLogBinaryMergerRateLognormB10_1")
        call Test%run(test_getLogBinaryMergerRateLognormM14_1, "test_getLogBinaryMergerRateLognormM14_1")
        call Test%run(test_getLogBinaryMergerRateLognormM17_1, "test_getLogBinaryMergerRateLognormM17_1")
        call Test%run(test_getLogBinaryMergerRateLognormF18_1, "test_getLogBinaryMergerRateLognormF18_1")
        call Test%finalize()
    end subroutine test_StarFormation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogRateH06() result(assertion)

        use Cosmology_mod, only: getLogLumDisWicMpc
        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: logzplus1 = log(zplus1)
        real(RK), parameter :: logRate_ref = 9.87062241049176_RK
        real(RK), parameter :: tolerance = 1.e-12_RK
        real(RK)            :: difference
        real(RK)            :: logRate

        logRate = getLogRateH06( zplus1 = zplus1, logzplus1 = logzplus1, twiceLogLumDisMpc = getLogLumDisWicMpc(zplus1) )
        difference = abs( (logRate - logRate_ref) / logRate_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate_ref   = ", logRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate       = ", logRate
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogRateH06

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogRateL08() result(assertion)

        use Cosmology_mod, only: getLogLumDisWicMpc
        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: logzplus1 = log(zplus1)
        real(RK), parameter :: logRate_ref = 11.5870199718159_RK
        real(RK), parameter :: tolerance = 1.e-12_RK
        real(RK)            :: difference
        real(RK)            :: logRate

        logRate = getLogRateL08( zplus1 = zplus1, logzplus1 = logzplus1, twiceLogLumDisMpc = getLogLumDisWicMpc(zplus1) )
        difference = abs( (logRate - logRate_ref) / logRate_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate_ref   = ", logRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate       = ", logRate
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogRateL08

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogRateB10() result(assertion)

        use Cosmology_mod, only: getLogLumDisWicMpc
        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: logzplus1 = log(zplus1)
        real(RK), parameter :: logRate_ref = 13.9081968356527_RK
        real(RK), parameter :: tolerance = 1.e-12_RK
        real(RK)            :: difference
        real(RK)            :: logRate

        logRate = getLogRateB10( zplus1 = zplus1, logzplus1 = logzplus1, twiceLogLumDisMpc = getLogLumDisWicMpc(zplus1) )
        difference = abs( (logRate - logRate_ref) / logRate_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate_ref   = ", logRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate       = ", logRate
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogRateB10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogRateM14() result(assertion)

        use Cosmology_mod, only: getLogLumDisWicMpc
        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: logzplus1 = log(zplus1)
        real(RK), parameter :: logRate_ref = 7.62065413866848_RK
        real(RK), parameter :: tolerance = 1.e-12_RK
        real(RK)            :: difference
        real(RK)            :: logRate

        logRate = getLogRateM14( zplus1 = zplus1, logzplus1 = logzplus1, twiceLogLumDisMpc = getLogLumDisWicMpc(zplus1) )
        difference = abs( (logRate - logRate_ref) / logRate_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate_ref   = ", logRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate       = ", logRate
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogRateM14

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogRateP15() result(assertion)

        use Cosmology_mod, only: getLogLumDisWicMpc
        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: logzplus1 = log(zplus1)
        real(RK), parameter :: logRate_ref = 7.87332272998867_RK
        real(RK), parameter :: tolerance = 1.e-12_RK
        real(RK)            :: difference
        real(RK)            :: logRate

        logRate = getLogRateP15( zplus1 = zplus1, logzplus1 = logzplus1, twiceLogLumDisMpc = getLogLumDisWicMpc(zplus1) )
        difference = abs( (logRate - logRate_ref) / logRate_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate_ref   = ", logRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate       = ", logRate
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogRateP15

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogRateM17() result(assertion)

        use Cosmology_mod, only: getLogLumDisWicMpc
        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: logzplus1 = log(zplus1)
        real(RK), parameter :: logRate_ref = 6.85265527251550_RK
        real(RK), parameter :: tolerance = 1.e-12_RK
        real(RK)            :: difference
        real(RK)            :: logRate

        logRate = getLogRateM17( zplus1 = zplus1, logzplus1 = logzplus1, twiceLogLumDisMpc = getLogLumDisWicMpc(zplus1) )
        difference = abs( (logRate - logRate_ref) / logRate_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate_ref   = ", logRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate       = ", logRate
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogRateM17

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogRateF18() result(assertion)

        use Cosmology_mod, only: getLogLumDisWicMpc
        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: logzplus1 = log(zplus1)
        real(RK), parameter :: logRate_ref = 6.81074639841109_RK
        real(RK), parameter :: tolerance = 1.e-12_RK
        real(RK)            :: difference
        real(RK)            :: logRate

        logRate = getLogRateF18( zplus1 = zplus1, logzplus1 = logzplus1, twiceLogLumDisMpc = getLogLumDisWicMpc(zplus1) )
        difference = abs( (logRate - logRate_ref) / logRate_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate_ref   = ", logRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "logRate       = ", logRate
            write(Test%outputUnit,"(*(g0,:,' '))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogRateF18

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBinaryMergerRateS15_1() result(assertion)

        use Constants_mod, only: RK

        implicit none

        logical :: assertion
        real(RK)    , parameter :: z = 1.e-2_RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: binaryMergerRateS15_ref = .000000000000000_RK
        real(RK)                :: binaryMergerRateS15
        real(RK)                :: difference

        binaryMergerRateS15 = getBinaryMergerRateS15(z)

        difference = abs(binaryMergerRateS15 - binaryMergerRateS15_ref)
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15_ref    =", binaryMergerRateS15_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15        =", binaryMergerRateS15
            write(Test%outputUnit,"(*(g0,:,' '))") "difference                 =", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getBinaryMergerRateS15_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBinaryMergerRateS15_2() result(assertion)

        use Constants_mod, only: RK

        implicit none

        logical :: assertion
        real(RK)    , parameter :: z = 2._RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: binaryMergerRateS15_ref = .1242102965042102E-01_RK
        real(RK)                :: binaryMergerRateS15
        real(RK)                :: difference

        binaryMergerRateS15 = getBinaryMergerRateS15(z)

        difference = abs( (binaryMergerRateS15 - binaryMergerRateS15_ref) / binaryMergerRateS15_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15_ref    =", binaryMergerRateS15_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15        =", binaryMergerRateS15
            write(Test%outputUnit,"(*(g0,:,' '))") "difference                 =", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getBinaryMergerRateS15_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBinaryMergerRateS15_3() result(assertion)

        use Constants_mod, only: RK

        implicit none

        logical :: assertion
        real(RK)    , parameter :: z = 4._RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: binaryMergerRateS15_ref = .2658922023902166E-02_RK
        real(RK)                :: binaryMergerRateS15
        real(RK)                :: difference

        binaryMergerRateS15 = getBinaryMergerRateS15(z)

        difference = abs( (binaryMergerRateS15 - binaryMergerRateS15_ref) / binaryMergerRateS15_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15_ref    =", binaryMergerRateS15_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15        =", binaryMergerRateS15
            write(Test%outputUnit,"(*(g0,:,' '))") "difference                 =", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getBinaryMergerRateS15_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBinaryMergerRateS15_4() result(assertion)

        use Constants_mod, only: RK

        implicit none

        logical :: assertion
        real(RK)    , parameter :: z = 10._RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: binaryMergerRateS15_ref = 0._RK
        real(RK)                :: binaryMergerRateS15
        real(RK)                :: difference

        binaryMergerRateS15 = getBinaryMergerRateS15(z)

        difference = abs(binaryMergerRateS15 - binaryMergerRateS15_ref)
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15_ref    =", binaryMergerRateS15_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "binaryMergerRateS15        =", binaryMergerRateS15
            write(Test%outputUnit,"(*(g0,:,' '))") "difference                 =", difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getBinaryMergerRateS15_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogBinaryMergerRateLognormH06_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: Logzplus1(*) = [0.01_RK, 0.1_RK, 0.5_RK, 1._RK, 1.6_RK, 2._RK, 4._RK]
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ .000000000000000_RK &
                                                                    , -9.522643068294865_RK &
                                                                    , -5.516405042254531_RK &
                                                                    , -4.774988408220006_RK &
                                                                    , -5.693972784881680_RK &
                                                                    , -9.316295406531074_RK &
                                                                    , .000000000000000_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Logzplus1))
        real(RK)                    :: Difference(size(Logzplus1))

        do i = 1, size(Logzplus1)
            LogBinaryMergerRate(i) = getLogBinaryMergerRateLognormH06(logzplus1 = Logzplus1(i))
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref)
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogBinaryMergerRateLognormH06_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogBinaryMergerRateLognormL08_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: Logzplus1(*) = [0.01_RK, 0.1_RK, 0.5_RK, 1._RK, 1.5_RK, 2._RK, 4._RK]
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ .000000000000000_RK &
                                                                    , -9.652689811083633_RK &
                                                                    , -5.676407092361923_RK &
                                                                    , -4.791335056013564_RK &
                                                                    , -5.330328192868819_RK &
                                                                    , -8.471291899297817_RK &
                                                                    , .000000000000000_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Logzplus1))
        real(RK)                    :: Difference(size(Logzplus1))

        do i = 1, size(Logzplus1)
            LogBinaryMergerRate(i) = getLogBinaryMergerRateLognormL08(logzplus1 = Logzplus1(i))
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref)
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogBinaryMergerRateLognormL08_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogBinaryMergerRateLognormB10_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: Logzplus1(*) = [0.01_RK, 0.1_RK, 0.5_RK, 1._RK, 1.5_RK, 2._RK, 4._RK]
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ .000000000000000_RK &
                                                                    , -10.40986473440980_RK &
                                                                    , -6.500787335196973_RK &
                                                                    , -5.213677851636544_RK &
                                                                    , -5.042789011502691_RK &
                                                                    , -7.221407286252206_RK &
                                                                    , .000000000000000_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Logzplus1))
        real(RK)                    :: Difference(size(Logzplus1))

        do i = 1, size(Logzplus1)
            LogBinaryMergerRate(i) = getLogBinaryMergerRateLognormB10(logzplus1 = Logzplus1(i))
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref)
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogBinaryMergerRateLognormB10_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogBinaryMergerRateLognormM14_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: Logzplus1(*) = [0.01_RK, 0.1_RK, 0.5_RK, 1._RK, 1.5_RK, 2._RK, 4._RK]
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ .000000000000000_RK &
                                                                    , -9.329653089161376_RK &
                                                                    , -5.659712865486435_RK &
                                                                    , -4.534167866418038_RK &
                                                                    , -5.833262100335860_RK &
                                                                    , -8.197338596152411_RK &
                                                                    , .000000000000000_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Logzplus1))
        real(RK)                    :: Difference(size(Logzplus1))

        do i = 1, size(Logzplus1)
            LogBinaryMergerRate(i) = getLogBinaryMergerRateLognormM14(logzplus1 = Logzplus1(i))
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref)
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogBinaryMergerRateLognormM14_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogBinaryMergerRateLognormM17_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: Logzplus1(*) = [0.01_RK, 0.1_RK, 0.5_RK, 1._RK, 1.5_RK, 2._RK, 4._RK]
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ .000000000000000_RK &
                                                                    , -9.446407015893174_RK &
                                                                    , -5.792990136584658_RK &
                                                                    , -4.509045304093034_RK &
                                                                    , -5.807989007937440_RK &
                                                                    , -8.517966032962981_RK &
                                                                    , .000000000000000_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Logzplus1))
        real(RK)                    :: Difference(size(Logzplus1))

        do i = 1, size(Logzplus1)
            LogBinaryMergerRate(i) = getLogBinaryMergerRateLognormM17(logzplus1 = Logzplus1(i))
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref)
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogBinaryMergerRateLognormM17_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogBinaryMergerRateLognormF18_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        real(RK)    , parameter     :: Logzplus1(*) = [0.01_RK, 0.1_RK, 0.5_RK, 1._RK, 1.5_RK, 2._RK, 4._RK]
        real(RK)    , parameter     :: tolerance = 1.e-12_RK
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ .000000000000000_RK &
                                                                    , -9.190714850567522_RK &
                                                                    , -5.417287541056011_RK &
                                                                    , -4.435837002244674_RK &
                                                                    , -6.128328470549945_RK &
                                                                    , -8.694199479455726_RK &
                                                                    , .000000000000000_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Logzplus1))
        real(RK)                    :: Difference(size(Logzplus1))

        do i = 1, size(Logzplus1)
            LogBinaryMergerRate(i) = getLogBinaryMergerRateLognormF18(logzplus1 = Logzplus1(i))
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref)
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    end function test_getLogBinaryMergerRateLognormF18_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBinaryMergerRate_1() result(assertion)

        use Statistics_mod, only: getLogProbLogNorm
        use Constants_mod, only: RK, IK, HUGE_RK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: nRefinement = 5_IK
        real(RK)    , parameter     :: Zplus1(*) = exp([0.01_RK, 0.1_RK, 0.5_RK, 1._RK, 1.5_RK, 2._RK, 4._RK])
        real(RK)    , parameter     :: zplus1Max = HUGE_RK
        real(RK)    , parameter     :: maxRelativeError = 1.e-3_RK
        real(RK)    , parameter     :: tolerance = 10 * maxRelativeError
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ 47135672.01905259_RK &
                                                                    , 4039350253.898525_RK &
                                                                    , 28704208268.10552_RK &
                                                                    , 11140236299.90150_RK &
                                                                    , 2721517859.383554_RK &
                                                                    , 713346416.5410618_RK &
                                                                    , 936072.0838233130_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Zplus1))
        real(RK)                    :: Difference(size(Zplus1))

        do i = 1, size(Zplus1)
            LogBinaryMergerRate(i) = getBinaryMergerRate( zplus1  = Zplus1(i) &
                                                        , zplus1Max = zplus1Max &
                                                        , nRefinement = nRefinement &
                                                        , maxRelativeError = maxRelativeError &
                                                        , getMergerDelayTimePDF = getMergerDelayTimePDF &
                                                        , getStarFormationRateDensity = getLogRateDensityB10 &
                                                        )
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref) / LogBinaryMergerRate_ref
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    contains

        function getMergerDelayTimePDF(mergerDelayTime) result(mergerDelayTimePDF)
            implicit none
            real(RK), intent(in)    :: mergerDelayTime
            real(RK)                :: mergerDelayTimePDF
            real(RK), parameter     :: inverseVariance = 1._RK
            real(RK), parameter     :: logSqrtInverseVariance = log(sqrt(inverseVariance))
            mergerDelayTimePDF = getLogProbLogNorm  ( logMean = 1._RK &
                                                    , inverseVariance = inverseVariance &
                                                    , logSqrtInverseVariance = logSqrtInverseVariance &
                                                    , logPoint = mergerDelayTime &
                                                    )
        end function 

    end function test_getBinaryMergerRate_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBinaryMergerRate_2() result(assertion)

        use Statistics_mod, only: getLogProbLogNorm
        use Constants_mod, only: RK, IK, HUGE_RK

        implicit none

        logical                     :: assertion
        integer(IK)                 :: i
        integer(IK) , parameter     :: nRefinement = 7_IK
        real(RK)    , parameter     :: zplus1Max = HUGE_RK
        real(RK)    , parameter     :: maxRelativeError = 1.e-5_RK
        real(RK)    , parameter     :: tolerance = 50 * maxRelativeError
        real(RK)    , parameter     :: Zplus1(*) = exp( [ 0.01_RK &
                                                        !, 0.1_RK
                                                        !, 0.5_RK
                                                        !, 1._RK
                                                        !, 1.5_RK
                                                        !, 2._RK
                                                        !, 4._RK
                                                        ] )
        real(RK)    , parameter     :: LogBinaryMergerRate_ref(*) = [ 80262853.72557181_RK &
                                                                    !, 6496053990.625437_RK &
                                                                    !, 40843786839.47808_RK &
                                                                    !, 15009263609.39050_RK &
                                                                    !, 3431922966.674707_RK &
                                                                    !, 873502037.5052428_RK &
                                                                    !, 1218944.513809759_RK &
                                                                    ]
        real(RK)                    :: LogBinaryMergerRate(size(Zplus1))
        real(RK)                    :: Difference(size(Zplus1))

        do i = 1, size(Zplus1)
            LogBinaryMergerRate(i) = getBinaryMergerRate( zplus1  = Zplus1(i) &
                                                        !, zplus1Max = zplus1Max &
                                                        !, nRefinement = nRefinement &
                                                        !, maxRelativeError = maxRelativeError &
                                                        , getMergerDelayTimePDF = getMergerDelayTimePDF &
                                                        , getStarFormationRateDensity = getLogRateDensityB10 &
                                                        )
        end do

        Difference = abs(LogBinaryMergerRate - LogBinaryMergerRate_ref) / LogBinaryMergerRate_ref
        assertion = all(Difference < tolerance)

        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,' '))")
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate_ref    =", LogBinaryMergerRate_ref
            write(Test%outputUnit,"(*(g0,:,' '))") "LogBinaryMergerRate        =", LogBinaryMergerRate
            write(Test%outputUnit,"(*(g0,:,' '))") "Difference                 =", Difference
            write(Test%outputUnit,"(*(g0,:,' '))")
        end if

    contains

        function getMergerDelayTimePDF(mergerDelayTime) result(mergerDelayTimePDF)
            implicit none
            real(RK), intent(in)    :: mergerDelayTime
            real(RK)                :: mergerDelayTimePDF
            real(RK), parameter     :: inverseVariance = 1._RK
            real(RK), parameter     :: logSqrtInverseVariance = log(sqrt(inverseVariance))
            mergerDelayTimePDF = getLogProbLogNorm  ( logMean = 1._RK &
                                                    , inverseVariance = inverseVariance &
                                                    , logSqrtInverseVariance = logSqrtInverseVariance &
                                                    , logPoint = mergerDelayTime &
                                                    )
        end function 

    end function test_getBinaryMergerRate_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_StarFormation_mod