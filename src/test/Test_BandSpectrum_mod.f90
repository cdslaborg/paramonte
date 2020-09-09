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
        real(RK) :: Limit(2)    ! energy window for fluence computation in keV units
        real(RK) :: photonFluence     ! integral of the Band spectrum in the given energy window
        real(RK) :: energyFluence     ! integral of the Band spectrum in the given energy window, in units of keV
        real(RK) :: tolerance   ! acceptable tolerance (accuracy) in the numerical computation of fluence
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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_BandSpectrum()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call test_getEbreak()
        call test_getPhotonFluence()
        call test_getEnergyFluence()
        call test_getPhotonFluenceFromEnergyFluence()
        call Test%finalize()
    end subroutine test_BandSpectrum

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getEbreak()
        use Constants_mod, only: RK, IK
        implicit none
        real(RK)                :: ebrk, difference
        if (Test%Image%isFirst) call Test%testing("getEbreak()")
        ebrk = getEbreak(BAND_SPEC1%epk,BAND_SPEC1%alpha,BAND_SPEC1%beta)
        difference = 2*abs(ebrk-BAND_SPEC1%ebrk)/(ebrk+BAND_SPEC1%ebrk)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "Ebreak, Reference Ebreak, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") ebrk, BAND_SPEC1%ebrk, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < 1.e-7_RK
        call Test%verify()
    end subroutine test_getEbreak

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getPhotonFluence()
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err

        if (Test%Image%isFirst) call Test%testing("getPhotonFluence()")

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of both upper and upper tails
        call getPhotonFluence   ( lowerLim      = BAND_SPEC1%Limit(1)   &
                                , upperLim      = BAND_SPEC1%Limit(2)   &
                                , epk           = BAND_SPEC1%epk        &
                                , alpha         = BAND_SPEC1%alpha      &
                                , beta          = BAND_SPEC1%beta       &
                                , tolerance     = BAND_SPEC1%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( photonFluence - BAND_SPEC1%photonFluence ) / ( photonFluence + BAND_SPEC1%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC1%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC1%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of only the upper tail
        call getPhotonFluence   ( lowerLim      = BAND_SPEC2%Limit(1)   &
                                , upperLim      = BAND_SPEC2%Limit(2)   &
                                , epk           = BAND_SPEC2%epk        &
                                , alpha         = BAND_SPEC2%alpha      &
                                , beta          = BAND_SPEC2%beta       &
                                , tolerance     = BAND_SPEC2%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( photonFluence - BAND_SPEC2%photonFluence ) / ( photonFluence + BAND_SPEC2%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC2%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of both upper and lower tails
        call getPhotonFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( photonFluence - BAND_SPEC3%photonFluence ) / ( photonFluence + BAND_SPEC3%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC3%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of both upper and upper tails with steep slopes
        call getPhotonFluence   ( lowerLim      = BAND_SPEC4%Limit(1)   &
                                , upperLim      = BAND_SPEC4%Limit(2)   &
                                , epk           = BAND_SPEC4%epk        &
                                , alpha         = BAND_SPEC4%alpha      &
                                , beta          = BAND_SPEC4%beta       &
                                , tolerance     = BAND_SPEC4%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( photonFluence - BAND_SPEC4%photonFluence ) / ( photonFluence + BAND_SPEC4%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC4%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC4%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    end subroutine test_getPhotonFluence

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getEnergyFluence()
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        real(RK)        :: energyFluence, difference
        type(Err_type)  :: Err

        if (Test%Image%isFirst) call Test%testing("getEnergyFluence()")

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of both upper and upper tails
        call getEnergyFluence   ( lowerLim      = BAND_SPEC1%Limit(1)   &
                                , upperLim      = BAND_SPEC1%Limit(2)   &
                                , epk           = BAND_SPEC1%epk        &
                                , alpha         = BAND_SPEC1%alpha      &
                                , beta          = BAND_SPEC1%beta       &
                                , tolerance     = BAND_SPEC1%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( energyFluence - BAND_SPEC1%energyFluence ) / ( energyFluence + BAND_SPEC1%energyFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC1%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC1%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of only the upper tail
        call getEnergyFluence   ( lowerLim      = BAND_SPEC2%Limit(1)   &
                                , upperLim      = BAND_SPEC2%Limit(2)   &
                                , epk           = BAND_SPEC2%epk        &
                                , alpha         = BAND_SPEC2%alpha      &
                                , beta          = BAND_SPEC2%beta       &
                                , tolerance     = BAND_SPEC2%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( energyFluence - BAND_SPEC2%energyFluence ) / ( energyFluence + BAND_SPEC2%energyFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC2%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC2%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of both upper and lower tails
        call getEnergyFluence   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( energyFluence - BAND_SPEC3%energyFluence ) / ( energyFluence + BAND_SPEC3%energyFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC3%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC3%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! test integration of both upper and upper tails with steep slopes
        call getEnergyFluence   ( lowerLim      = BAND_SPEC4%Limit(1)   &
                                , upperLim      = BAND_SPEC4%Limit(2)   &
                                , epk           = BAND_SPEC4%epk        &
                                , alpha         = BAND_SPEC4%alpha      &
                                , beta          = BAND_SPEC4%beta       &
                                , tolerance     = BAND_SPEC4%tolerance  &
                                , energyFluence = energyFluence         &
                                , Err           = Err                   &
                                )
        difference = 2 * abs( energyFluence - BAND_SPEC4%energyFluence ) / ( energyFluence + BAND_SPEC4%energyFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") energyFluence, BAND_SPEC4%energyFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC4%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    end subroutine test_getEnergyFluence

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getPhotonFluenceFromEnergyFluence()
        use Constants_mod, only: RK, IK
        use Err_mod, only: Err_type
        implicit none
        real(RK)        :: photonFluence, difference
        type(Err_type)  :: Err

        if (Test%Image%isFirst) call Test%testing("getPhotonFluenceFromEnergyFluence()")

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
        difference = 2 * abs( photonFluence - BAND_SPEC1%photonFluence ) / ( photonFluence + BAND_SPEC1%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC1%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC1%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
        difference = 2 * abs( photonFluence - BAND_SPEC2%photonFluence ) / ( photonFluence + BAND_SPEC2%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC2%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
        difference = 2 * abs( photonFluence - BAND_SPEC3%photonFluence ) / ( photonFluence + BAND_SPEC3%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC3%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC4%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC4%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
        difference = 2 * abs( photonFluence - BAND_SPEC2%photonFluence ) / ( photonFluence + BAND_SPEC2%photonFluence )
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0,:,', '))")
            write(Test%outputUnit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(Test%outputUnit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(Test%outputUnit,"(*(g0,:,', '))")
        end if
        Test%assertion = difference < BAND_SPEC4%tolerance
        call Test%verify()

        !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    end subroutine test_getPhotonFluenceFromEnergyFluence

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_BandSpectrum_mod