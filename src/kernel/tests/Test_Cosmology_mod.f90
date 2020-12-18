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

!>  \brief This module contains tests of the module [Cosmology_mod](@ref cosmology_mod).
!>  \author Amir Shahmoradi

module Test_Cosmology_mod

    use Cosmology_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Cosmology

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_Cosmology()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getlogdvdz, "test_getlogdvdz")
        call Test%run(test_ldiswickram, "test_ldiswickram")
        call Test%run(test_getLogLumDisWicMpc, "test_getLogLumDisWicMpc")
        call Test%run(test_getLookBackTime_1, "test_getLookBackTime_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getLookBackTime_2, "test_getLookBackTime_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getUniverseAgeDerivative, "test_getUniverseAgeDerivative")
        call Test%finalize()
    end subroutine test_Cosmology

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getlogdvdz() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: LOGzplus1 = log(zplus1)
        real(RK), parameter :: twiceLogLumDisMpc = 0._RK
        real(RK), parameter :: logdvdz_ref = 3.421655392580980_RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK)            :: difference
        real(RK)            :: logdvdz

        logdvdz = getlogdvdz(zplus1,LOGzplus1,twiceLogLumDisMpc)
        difference = abs( (logdvdz - logdvdz_ref) / logdvdz_ref )
        assertion = difference < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15))")
            write(Test%outputUnit,"(*(g0.15))") "logdvdz_ref   = ", logdvdz_ref
            write(Test%outputUnit,"(*(g0.15))") "logdvdz       = ", logdvdz
            write(Test%outputUnit,"(*(g0.15))") "difference    = ", difference
            write(Test%outputUnit,"(*(g0.15))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getlogdvdz

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_ldiswickram() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: tolerance = 1.e-8_RK
        real(RK), parameter :: lumDisWicMpc_ref = 90887.79805125466_RK
        real(RK)            :: lumDisWicMpc
        real(RK)            :: difference

        lumDisWicMpc = ldiswickram(zplus1)
        difference = abs( (lumDisWicMpc - lumDisWicMpc_ref) / lumDisWicMpc_ref )
        assertion = difference < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15))")
            write(Test%outputUnit,"(*(g0.15))") "lumDisWicMpc_ref  = ", lumDisWicMpc_ref
            write(Test%outputUnit,"(*(g0.15))") "lumDisWicMpc      = ", lumDisWicMpc
            write(Test%outputUnit,"(*(g0.15))") "difference        = ", difference
            write(Test%outputUnit,"(*(g0.15))")
        end if
        ! LCOV_EXCL_STOP

    end function test_ldiswickram

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLogLumDisWicMpc() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical :: assertion
        real(RK), parameter :: zplus1 = 1.e1_RK
        real(RK), parameter :: tolerance = 1.e-8_RK
        real(RK), parameter :: logLumDisWicMpc_ref = 11.417381027586867_RK
        real(RK)            :: logLumDisWicMpc
        real(RK)            :: difference

        logLumDisWicMpc = getLogLumDisWicMpc(zplus1)
        difference = abs( (logLumDisWicMpc - logLumDisWicMpc_ref) / logLumDisWicMpc_ref )
        assertion = difference < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15))")
            write(Test%outputUnit,"(*(g0.15))") "logLumDisWicMpc_ref   = ", logLumDisWicMpc_ref
            write(Test%outputUnit,"(*(g0.15))") "logLumDisWicMpc       = ", logLumDisWicMpc
            write(Test%outputUnit,"(*(g0.15))") "difference            = ", difference
            write(Test%outputUnit,"(*(g0.15))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLogLumDisWicMpc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLookBackTime_1() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: zplus1 = 1.e1_RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: lookBackTime_ref = 12.7740854156781_RK
        real(RK)    , parameter :: maxRelativeError = 1.e-5_RK
        integer(IK) , parameter :: nRefinement = 6
        real(RK)                :: lookBackTime
        real(RK)                :: difference

        lookBackTime = getLookBackTime(zplus1 = zplus1, maxRelativeError = maxRelativeError, nRefinement = nRefinement)
        difference = abs( (lookBackTime - lookBackTime_ref) / lookBackTime_ref )
        assertion = difference < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15))")
            write(Test%outputUnit,"(*(g0.15))") "lookBackTime_ref  = ", lookBackTime_ref
            write(Test%outputUnit,"(*(g0.15))") "lookBackTime      = ", lookBackTime
            write(Test%outputUnit,"(*(g0.15))") "difference        = ", difference
            write(Test%outputUnit,"(*(g0.15))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLookBackTime_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLookBackTime_2() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: zplus1 = 1.e1_RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: lookBackTime_ref = 12.7736258974511_RK
        real(RK)                :: lookBackTime
        real(RK)                :: difference

        lookBackTime = getLookBackTime(zplus1 = zplus1)
        difference = abs( (lookBackTime - lookBackTime_ref) / lookBackTime_ref )
        assertion = difference < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15))")
            write(Test%outputUnit,"(*(g0.15))") "lookBackTime_ref  = ", lookBackTime_ref
            write(Test%outputUnit,"(*(g0.15))") "lookBackTime      = ", lookBackTime
            write(Test%outputUnit,"(*(g0.15))") "difference        = ", difference
            write(Test%outputUnit,"(*(g0.15))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getLookBackTime_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getUniverseAgeDerivative() result(assertion)

        use Constants_mod, only: RK, IK

        implicit none

        logical                 :: assertion
        real(RK)    , parameter :: zplus1 = 1.e1_RK
        real(RK)    , parameter :: tolerance = 1.e-12_RK
        real(RK)    , parameter :: universeAgeDerivative_ref = 8.122223525986105e-05_RK
        real(RK)                :: universeAgeDerivative
        real(RK)                :: difference

        universeAgeDerivative = getUniverseAgeDerivative(zplus1 = zplus1)
        difference = abs( (universeAgeDerivative - universeAgeDerivative_ref) / universeAgeDerivative_ref )
        assertion = difference < tolerance

        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0.15))")
            write(Test%outputUnit,"(*(g0.15))") "universeAgeDerivative_ref = ", universeAgeDerivative_ref
            write(Test%outputUnit,"(*(g0.15))") "universeAgeDerivative     = ", universeAgeDerivative
            write(Test%outputUnit,"(*(g0.15))") "difference                = ", difference
            write(Test%outputUnit,"(*(g0.15))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getUniverseAgeDerivative

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_Cosmology_mod ! LCOV_EXCL_LINE