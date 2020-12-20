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

!>  \brief This module contains tests of the module [BandSpectrum_mod](@ref bandspectrum_mod).
!>  \author Amir Shahmoradi

module Test_BandSpectrum_mod

    use BandSpectrum_mod
    use Test_mod, only: Test_type
    !use Constants_mod, only: RK, IK
    implicit none

    private
    public :: test_BandSpectrum

    type(Test_type) :: Test

    type :: BandSpec_type
        real(RK) :: epk, alpha, beta, ebrk, coef
        real(RK) :: Limit(2)        ! energy window for fluence computation in keV units
        real(RK) :: photonFluence   ! integral of the Band spectrum in the given energy window
        real(RK) :: energyFluence   ! integral of the Band spectrum in the given energy window, in units of keV
        real(RK) :: tolerance       ! acceptable tolerance (accuracy) in the numerical computation of fluence
    end type BandSpec_type

    ! involves integration of both upper and lower tails

    type(BandSpec_type), parameter :: BAND_SPEC1 = BandSpec_type( epk = 700._RK &
                                                                , alpha = -0.5_RK &
                                                                , beta = -2.5_RK &
                                                                , ebrk = 9.333333333333334e2_RK &
                                                                , coef = 1.178920689527826e5_RK &
                                                                , Limit = [1._RK,10000._RK] &
                                                                , photonFluence = 37.226409565894123_RK &
                                                                , energyFluence = 1.195755906912896e4_RK &
                                                                , tolerance = 1.e-7_RK &
                                                                )

    ! involves only the lower tail integration

    type(BandSpec_type), parameter :: BAND_SPEC2 = BandSpec_type( epk = 700._RK &
                                                                , alpha = -0.5_RK &
                                                                , beta = -2.5_RK &
                                                                , ebrk = 9.333333333333334e2_RK &
                                                                , coef = 1.178920689527826e5_RK &
                                                                , Limit = [1._RK,500._RK] &
                                                                , photonFluence = 30.806431300618090_RK &
                                                                , energyFluence = 4.079656304178462e3_RK &
                                                                , tolerance = 1.e-7_RK &
                                                                )

    ! involves only the upper tail integration

    type(BandSpec_type), parameter :: BAND_SPEC3 = BandSpec_type( epk = 700._RK &
                                                                , alpha = -0.5_RK &
                                                                , beta = -2.5_RK &
                                                                , ebrk = 9.333333333333334e2_RK &
                                                                , coef = 1.178920689527826e5_RK &
                                                                , Limit = [1000._RK,10000._RK] &
                                                                , photonFluence = 2.406788327100909_RK &
                                                                , energyFluence = 5.098307740152641e3_RK &
                                                                , tolerance = 1.e-7_RK &
                                                                )

    ! involves integration of both upper and lower tails

    type(BandSpec_type), parameter :: BAND_SPEC4 = BandSpec_type( epk = 700._RK &
                                                                , alpha = -1.9_RK &
                                                                , beta = -3.5_RK &
                                                                , ebrk = 1.119999999999999e4_RK &
                                                                , coef = 6.079637221508616e5_RK &
                                                                , Limit = [1._RK,10000._RK] &
                                                                , photonFluence = 1.108858577351433_RK &
                                                                , energyFluence = 12.769650780448401_RK &
                                                                , tolerance = 1.e-6_RK &
                                                                )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_BandSpectrum()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getEbreak, "test_getEbreak")
        call Test%run(test_getBandParam_1, "test_getBandParam_1")
        call Test%run(test_getPhotonFlux_1, "test_getPhotonFlux_1")
        call Test%run(test_getPhotonFlux_2, "test_getPhotonFlux_2")
        call Test%run(test_getPhotonFlux_3, "test_getPhotonFlux_3")
        call Test%run(test_getPhotonFlux_4, "test_getPhotonFlux_4")
        call Test%run(test_getPhotonFluence_3, "test_getPhotonFluence_3")
        call Test%run(test_getPhotonFluence_5, "test_getPhotonFluence_5")
        call Test%run(test_getPhotonFluence_6, "test_getPhotonFluence_6")
        call Test%run(test_getPhotonFluence_7, "test_getPhotonFluence_7")
        call Test%run(test_getPhotonFluence_8, "test_getPhotonFluence_8")
        call Test%run(test_getEnergyFluence_3, "test_getEnergyFluence_3")
        call Test%run(test_getEnergyFluence_6, "test_getEnergyFluence_6")
        call Test%run(test_getEnergyFluence_7, "test_getEnergyFluence_7")
        call Test%run(test_getEnergyFluence_8, "test_getEnergyFluence_8")
        call Test%run(test_getPhotonFluxLower_1, "test_getPhotonFluxLower_1")
        call Test%run(test_getPhotonFluenceFromEnergyFluence_6, "test_getPhotonFluenceFromEnergyFluence_6")
        call Test%run(test_getPhotonFluenceFromEnergyFluence_7, "test_getPhotonFluenceFromEnergyFluence_7")
        call Test%run(test_getPhotonFluenceFromEnergyFluence_8, "test_getPhotonFluenceFromEnergyFluence_8")
        call Test%run(test_getPhotonFluence_1, "test_getPhotonFluence_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getPhotonFluence_2, "test_getPhotonFluence_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getPhotonFluence_4, "test_getPhotonFluence_4") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getEnergyFluence_1, "test_getEnergyFluence_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getEnergyFluence_2, "test_getEnergyFluence_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getEnergyFluence_4, "test_getEnergyFluence_4") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getPhotonFluenceFromEnergyFluence_1, "test_getPhotonFluenceFromEnergyFluence_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getPhotonFluenceFromEnergyFluence_2, "test_getPhotonFluenceFromEnergyFluence_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getPhotonFluenceFromEnergyFluence_3, "test_getPhotonFluenceFromEnergyFluence_3") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getPhotonFluenceFromEnergyFluence_4, "test_getPhotonFluenceFromEnergyFluence_4") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%run(test_getPhotonFluenceFromEnergyFluence_5, "test_getPhotonFluenceFromEnergyFluence_5") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call Test%finalize()
    end subroutine test_BandSpectrum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getEbreak() result(assertion)
        use Constants_mod, only: RK, IK
        implicit none
        logical  :: assertion
        real(RK) :: ebrk, difference
        ebrk = getEbreak(BAND_SPEC1%epk,BAND_SPEC1%alpha,BAND_SPEC1%beta)
        difference = 2._RK * abs(ebrk - BAND_SPEC1%ebrk) / (ebrk + BAND_SPEC1%ebrk)
        assertion = difference < 1.e-7_RK
        ! LCOV_EXCL_START
        if (Test%isDebugMode .and. .not. assertion) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Ebreak, Reference Ebreak, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") ebrk, BAND_SPEC1%ebrk, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEbreak

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandParam_1() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        real(RK), parameter :: alphaPlusTwo_ref = 0.5_RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK), parameter :: ebrk_ref = 6.e2_RK
        real(RK), parameter :: coef_ref = 220.72766470286541_RK
        real(RK), parameter :: alpha = -1.5_RK
        real(RK), parameter :: beta = -2.5_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical  :: assertion
        real(RK) :: difference
        real(RK) :: alphaPlusTwo
        real(RK) :: ebrk
        real(RK) :: coef

        call getBandParam(epk,alpha,beta,ebrk,coef,alphaPlusTwo)

        difference = abs( (ebrk - ebrk_ref) / ebrk_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "ebrk_ref   ", ebrk_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "ebrk       ", ebrk
            write(Test%outputUnit,"(*(g0,:,', '))") "difference ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        difference = abs( (coef - coef_ref) / coef_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "coef_ref   ", coef_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "coef       ", coef
            write(Test%outputUnit,"(*(g0,:,', '))") "difference ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        difference = abs( (alphaPlusTwo - alphaPlusTwo_ref) / alphaPlusTwo_ref )
        assertion = difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "alphaPlusTwo_ref   ", alphaPlusTwo_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "alphaPlusTwo       ", alphaPlusTwo
            write(Test%outputUnit,"(*(g0,:,', '))") "difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getBandParam_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFluxLower_1() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        real(RK), parameter :: photonFluxLower_ref = 0.59727709940405714E-5_RK
        real(RK), parameter :: alphaPlusTwo = 0.5_RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK), parameter :: energy = 1.e3_RK
        real(RK), parameter :: alpha = -1.5_RK
        real(RK), parameter :: beta = -2.5_RK
        real(RK), parameter :: coef = 220.72766470286541_RK
        real(RK), parameter :: ebrk = 6.e2_RK
        real(RK), parameter :: epk = 3.e2_RK
        real(RK), parameter :: alphaPlusTwoOverEpk = alphaPlusTwo / epk

        logical  :: assertion
        real(RK) :: difference
        real(RK) :: photonFluxLower

        photonFluxLower = getPhotonFluxLower(energy,alpha,alphaPlusTwoOverEpk)
        assertion = photonFluxLower > 0._RK

        difference = abs( (photonFluxLower - photonFluxLower_ref) / photonFluxLower_ref )
        assertion = assertion .and. difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFluxLower_ref", photonFluxLower_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFluxLower    ", photonFluxLower
            write(Test%outputUnit,"(*(g0,:,', '))") "difference         ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFluxLower_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFlux_1() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        real(RK), parameter :: photonFlux_ref = 0.69800216307100779E-5_RK
        real(RK), parameter :: alphaPlusTwo = 0.5_RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK), parameter :: energy = 1.e3_RK
        real(RK), parameter :: alpha = -1.5_RK
        real(RK), parameter :: beta = -2.5_RK
        real(RK), parameter :: coef = 220.72766470286541_RK
        real(RK), parameter :: ebrk = 6.e2_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical  :: assertion
        real(RK) :: difference
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux > 0._RK

        difference = abs( (photonFlux - photonFlux_ref) / photonFlux_ref )
        assertion = assertion .and. difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFlux_ref ", photonFlux_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFlux_2() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        real(RK), parameter :: photonFlux_ref = 0.31100098078334186E-1_RK
        real(RK), parameter :: alphaPlusTwo = 0.5_RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK), parameter :: energy = 1.e1_RK
        real(RK), parameter :: alpha = -1.5_RK
        real(RK), parameter :: beta = -2.5_RK
        real(RK), parameter :: coef = 220.72766470286541_RK
        real(RK), parameter :: ebrk = 6.e2_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical  :: assertion
        real(RK) :: difference
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux > 0._RK

        difference = abs( (photonFlux - photonFlux_ref) / photonFlux_ref )
        assertion = assertion .and. difference < tolerance

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFlux_ref ", photonFlux_ref
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(Test%outputUnit,"(*(g0,:,', '))") "difference     ", difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! alpha < -2
    function test_getPhotonFlux_3() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        real(RK), parameter :: alphaPlusTwo = 0.5_RK
        real(RK), parameter :: energy = 1.e1_RK
        real(RK), parameter :: alpha = -3.5_RK
        real(RK), parameter :: beta = -2.5_RK
        real(RK), parameter :: coef = 220.72766470286541_RK
        real(RK), parameter :: ebrk = 6.e2_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical  :: assertion
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux < 0._RK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! alpha < beta
    function test_getPhotonFlux_4() result(assertion)

        use Constants_mod, only: RK, IK
        implicit none

        real(RK), parameter :: alphaPlusTwo = 0.5_RK
        real(RK), parameter :: energy = 1.e1_RK
        real(RK), parameter :: alpha = -1.5_RK
        real(RK), parameter :: beta = -0.5_RK
        real(RK), parameter :: coef = 220.72766470286541_RK
        real(RK), parameter :: ebrk = 6.e2_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical  :: assertion
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux < 0._RK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of both upper and lower tails.
    function test_getPhotonFluence_3() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC3%photonFluence ) / ( photonFluence + BAND_SPEC3%photonFluence )
        assertion = difference < BAND_SPEC3%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluence_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getPhotonFluence_5() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: tolerance = 1.e-10_RK
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluence   ( lowerLim      = BAND_SPEC3%Limit(2)   &
                                , upperLim      = BAND_SPEC3%Limit(1)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = abs( photonFluence - 0._RK )
        assertion = difference < tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluence_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getPhotonFluence_6() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence
        type(Err_type)  :: Err
        call getPhotonFluence   ( lowerLim      = BAND_SPEC3%Limit(2)   &
                                , upperLim      = BAND_SPEC3%Limit(1)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        assertion = abs(photonFluence - 0._RK) < 1.e-12_RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence", photonFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluence_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test with conflicting alpha photon index.
    function test_getPhotonFluence_7() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence
        type(Err_type)  :: Err
        call getPhotonFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = -1.e1_RK              &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        assertion = Err%occurred .and. photonFluence < 0._RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence", photonFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluence_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test with conflicting alpha < beta photon indices.
    function test_getPhotonFluence_8() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%beta       &
                                , beta          = BAND_SPEC3%alpha      &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        assertion = Err%occurred .and. photonFluence < 0._RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluence_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of both upper and lower tails.
    function test_getEnergyFluence_3() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: energyFluence, difference
        type(Err_type)  :: Err
        call getEnergyFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( energyFluence - BAND_SPEC3%energyFluence ) / ( energyFluence + BAND_SPEC3%energyFluence )
        assertion = difference < BAND_SPEC3%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC3%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEnergyFluence_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getEnergyFluence_6() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: EnergyFluence
        type(Err_type)  :: Err
        call getEnergyFluence   ( lowerLim      = BAND_SPEC3%Limit(2)   &
                                , upperLim      = BAND_SPEC3%Limit(1)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , EnergyFluence = EnergyFluence         &
                                , Err           = Err                   &
                                )
        assertion = abs(EnergyFluence - 0._RK) < 1.e-12_RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "EnergyFluence", EnergyFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getEnergyFluence_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test with conflicting alpha photon index.
    function test_getEnergyFluence_7() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: EnergyFluence
        type(Err_type)  :: Err
        call getEnergyFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = -1.e1_RK              &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , EnergyFluence = EnergyFluence         &
                                , Err           = Err                   &
                                )
        assertion = Err%occurred .and. EnergyFluence < 0._RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "EnergyFluence", EnergyFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getEnergyFluence_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test with conflicting alpha < beta photon indices.
    function test_getEnergyFluence_8() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: EnergyFluence
        type(Err_type)  :: Err
        call getEnergyFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%beta       &
                                , beta          = BAND_SPEC3%alpha      &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , EnergyFluence = EnergyFluence         &
                                , Err           = Err                   &
                                )
        assertion = Err%occurred .and. EnergyFluence < 0._RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "EnergyFluence", EnergyFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getEnergyFluence_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getPhotonFluenceFromEnergyFluence_6() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
                                                , lowerLim      = BAND_SPEC2%Limit(2)       &
                                                , upperLim      = BAND_SPEC2%Limit(1)       &
                                                , epk           = BAND_SPEC2%epk            &
                                                , alpha         = BAND_SPEC2%alpha          &
                                                , beta          = BAND_SPEC2%alpha          &
                                                , tolerance     = BAND_SPEC2%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        assertion = abs(photonFluence - 0._RK) < 1.e-12_RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFluence", photonFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluenceFromEnergyFluence_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test with conflicting alpha photon index alpha < -2.
    function test_getPhotonFluenceFromEnergyFluence_7() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
                                                , lowerLim      = BAND_SPEC2%Limit(1)       &
                                                , upperLim      = BAND_SPEC2%Limit(2)       &
                                                , epk           = BAND_SPEC2%epk            &
                                                , alpha         = -1.e1_RK                  &
                                                , beta          = BAND_SPEC2%beta           &
                                                , tolerance     = BAND_SPEC2%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        assertion = Err%occurred .and. photonFluence < 0._RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFluence", photonFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluenceFromEnergyFluence_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test with conflicting alpha < beta photon indices.
    function test_getPhotonFluenceFromEnergyFluence_8() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
                                                , lowerLim      = BAND_SPEC2%Limit(1)       &
                                                , upperLim      = BAND_SPEC2%Limit(2)       &
                                                , epk           = BAND_SPEC2%epk            &
                                                , alpha         = BAND_SPEC2%beta           &
                                                , beta          = BAND_SPEC2%alpha          &
                                                , tolerance     = BAND_SPEC2%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        assertion = Err%occurred .and. photonFluence < 0._RK
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photonFluence", photonFluence
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getPhotonFluenceFromEnergyFluence_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if !defined OS_IS_WSL || !defined CODECOV_ENABLED || defined DLL_ENABLED

    !> \brief
    !> Test the integration of both the upper and lower tails.
    function test_getPhotonFluence_1() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluence   ( lowerLim      = BAND_SPEC1%Limit(1)   &
                                , upperLim      = BAND_SPEC1%Limit(2)   &
                                , epk           = BAND_SPEC1%epk        &
                                , alpha         = BAND_SPEC1%alpha      &
                                , beta          = BAND_SPEC1%beta       &
                                , tolerance     = BAND_SPEC1%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC1%photonFluence ) / ( photonFluence + BAND_SPEC1%photonFluence )
        assertion = difference < BAND_SPEC1%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC1%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluence_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of only the upper tail.
    function test_getPhotonFluence_2() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluence   ( lowerLim      = BAND_SPEC2%Limit(1)   &
                                , upperLim      = BAND_SPEC2%Limit(2)   &
                                , epk           = BAND_SPEC2%epk        &
                                , alpha         = BAND_SPEC2%alpha      &
                                , beta          = BAND_SPEC2%beta       &
                                , tolerance     = BAND_SPEC2%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC2%photonFluence ) / ( photonFluence + BAND_SPEC2%photonFluence )
        assertion = difference < BAND_SPEC2%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluence_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of both upper and upper tails with steep slopes.
    function test_getPhotonFluence_4() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err

        call getPhotonFluence   ( lowerLim      = BAND_SPEC4%Limit(1)   &
                                , upperLim      = BAND_SPEC4%Limit(2)   &
                                , epk           = BAND_SPEC4%epk        &
                                , alpha         = BAND_SPEC4%alpha      &
                                , beta          = BAND_SPEC4%beta       &
                                , tolerance     = BAND_SPEC4%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC4%photonFluence ) / ( photonFluence + BAND_SPEC4%photonFluence )
        assertion = difference < BAND_SPEC4%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC4%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluence_4

    !> \brief
    !> Test the integration of both upper and upper tails.
    function test_getEnergyFluence_1() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: energyFluence, difference
        type(Err_type)  :: Err
        call getEnergyFluence   ( lowerLim      = BAND_SPEC1%Limit(1)   &
                                , upperLim      = BAND_SPEC1%Limit(2)   &
                                , epk           = BAND_SPEC1%epk        &
                                , alpha         = BAND_SPEC1%alpha      &
                                , beta          = BAND_SPEC1%beta       &
                                , tolerance     = BAND_SPEC1%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( energyFluence - BAND_SPEC1%energyFluence ) / ( energyFluence + BAND_SPEC1%energyFluence )
        assertion = difference < BAND_SPEC1%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC1%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEnergyFluence_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of only the upper tail.
    function test_getEnergyFluence_2() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: energyFluence, difference
        type(Err_type)  :: Err
        call getEnergyFluence   ( lowerLim      = BAND_SPEC2%Limit(1)   &
                                , upperLim      = BAND_SPEC2%Limit(2)   &
                                , epk           = BAND_SPEC2%epk        &
                                , alpha         = BAND_SPEC2%alpha      &
                                , beta          = BAND_SPEC2%beta       &
                                , tolerance     = BAND_SPEC2%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( energyFluence - BAND_SPEC2%energyFluence ) / ( energyFluence + BAND_SPEC2%energyFluence )
        assertion = difference < BAND_SPEC2%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC2%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEnergyFluence_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the integration of both upper and upper tails with steep slopes.
    function test_getEnergyFluence_4() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: energyFluence, difference
        type(Err_type)  :: Err
        call getEnergyFluence   ( lowerLim      = BAND_SPEC4%Limit(1)   &
                                , upperLim      = BAND_SPEC4%Limit(2)   &
                                , epk           = BAND_SPEC4%epk        &
                                , alpha         = BAND_SPEC4%alpha      &
                                , beta          = BAND_SPEC4%beta       &
                                , tolerance     = BAND_SPEC4%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2._RK * abs( energyFluence - BAND_SPEC4%energyFluence ) / ( energyFluence + BAND_SPEC4%energyFluence )
        assertion = difference < BAND_SPEC4%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC4%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getEnergyFluence_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFluenceFromEnergyFluence_1() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC1%energyFluence  &
                                                , lowerLim      = BAND_SPEC1%Limit(1)       &
                                                , upperLim      = BAND_SPEC1%Limit(2)       &
                                                , epk           = BAND_SPEC1%epk            &
                                                , alpha         = BAND_SPEC1%alpha          &
                                                , beta          = BAND_SPEC1%beta           &
                                                , tolerance     = BAND_SPEC1%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC1%photonFluence ) / ( photonFluence + BAND_SPEC1%photonFluence )
        assertion = difference < BAND_SPEC1%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC1%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluenceFromEnergyFluence_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFluenceFromEnergyFluence_2() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
                                                , lowerLim      = BAND_SPEC2%Limit(1)       &
                                                , upperLim      = BAND_SPEC2%Limit(2)       &
                                                , epk           = BAND_SPEC2%epk            &
                                                , alpha         = BAND_SPEC2%alpha          &
                                                , beta          = BAND_SPEC2%beta           &
                                                , tolerance     = BAND_SPEC2%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC2%photonFluence ) / ( photonFluence + BAND_SPEC2%photonFluence )
        assertion = difference < BAND_SPEC2%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluenceFromEnergyFluence_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFluenceFromEnergyFluence_3() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC3%energyFluence  &
                                                , lowerLim      = BAND_SPEC3%Limit(1)       &
                                                , upperLim      = BAND_SPEC3%Limit(2)       &
                                                , epk           = BAND_SPEC3%epk            &
                                                , alpha         = BAND_SPEC3%alpha          &
                                                , beta          = BAND_SPEC3%beta           &
                                                , tolerance     = BAND_SPEC3%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC3%photonFluence ) / ( photonFluence + BAND_SPEC3%photonFluence )
        assertion = difference < BAND_SPEC3%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluenceFromEnergyFluence_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFluenceFromEnergyFluence_4() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC4%energyFluence  &
                                                , lowerLim      = BAND_SPEC4%Limit(1)       &
                                                , upperLim      = BAND_SPEC4%Limit(2)       &
                                                , epk           = BAND_SPEC4%epk            &
                                                , alpha         = BAND_SPEC4%alpha          &
                                                , beta          = BAND_SPEC4%beta           &
                                                , tolerance     = BAND_SPEC4%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        difference = 2 * abs( photonFluence - BAND_SPEC4%photonFluence ) / ( photonFluence + BAND_SPEC4%photonFluence )
        assertion = difference < BAND_SPEC4%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC4%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluenceFromEnergyFluence_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFluenceFromEnergyFluence_5() result(assertion)
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        logical         :: assertion
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err
        call getPhotonFluenceFromEnergyFluence  ( energyFluence = BAND_SPEC1%energyFluence  &
                                                , lowerLim      = BAND_SPEC1%Limit(1)       &
                                                , upperLim      = BAND_SPEC1%Limit(2)       &
                                                , epk           = BAND_SPEC1%epk            &
                                                , alpha         = BAND_SPEC1%alpha          &
                                                , beta          = BAND_SPEC1%beta           &
                                                , tolerance     = BAND_SPEC1%tolerance      &
                                                , lowerLimNew   = BAND_SPEC2%Limit(1)       &
                                                , upperLimNew   = BAND_SPEC2%Limit(2)       &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        difference = 2._RK * abs( photonFluence - BAND_SPEC2%photonFluence ) / ( photonFluence + BAND_SPEC2%photonFluence )
        assertion = difference < BAND_SPEC4%tolerance
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getPhotonFluenceFromEnergyFluence_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#endif

end module Test_BandSpectrum_mod ! LCOV_EXCL_LINE