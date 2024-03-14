!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [pm_distBand](@ref pm_distBand).
!>  \author Amir Shahmoradi

module test_pm_distBand

    use pm_distBand
    use pm_test, only: test_type, LK
    !use pm_kind, only: RK, IK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

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

    subroutine setTest()

        test = test_type(MODULE_NAME)
        call test%run(test_getBandEbreak, SK_"test_getBandEbreak")
        call test%run(test_getBandParam_1, SK_"test_getBandParam_1")
        call test%run(test_getPhotonFlux_1, SK_"test_getPhotonFlux_1")
        call test%run(test_getPhotonFlux_2, SK_"test_getPhotonFlux_2")
        call test%run(test_getPhotonFlux_3, SK_"test_getPhotonFlux_3")
        call test%run(test_getPhotonFlux_4, SK_"test_getPhotonFlux_4")
        call test%run(test_getBandPhoton_1, SK_"test_getBandPhoton_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getBandPhoton_2, SK_"test_getBandPhoton_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getBandPhoton_3, SK_"test_getBandPhoton_3")
        call test%run(test_getBandPhoton_4, SK_"test_getBandPhoton_4") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getBandPhoton_5, SK_"test_getBandPhoton_5")
        call test%run(test_getBandPhoton_6, SK_"test_getBandPhoton_6")
        call test%run(test_getBandPhoton_7, SK_"test_getBandPhoton_7")
        call test%run(test_getBandPhoton_8, SK_"test_getBandPhoton_8")
        call test%run(test_getFluenceEnergy_1, SK_"test_getFluenceEnergy_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getFluenceEnergy_2, SK_"test_getFluenceEnergy_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getFluenceEnergy_3, SK_"test_getFluenceEnergy_3")
        call test%run(test_getFluenceEnergy_4, SK_"test_getFluenceEnergy_4") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getFluenceEnergy_6, SK_"test_getFluenceEnergy_6")
        call test%run(test_getFluenceEnergy_7, SK_"test_getFluenceEnergy_7")
        call test%run(test_getFluenceEnergy_8, SK_"test_getFluenceEnergy_8")
        call test%run(test_getPhotonFluxLower_1, SK_"test_getPhotonFluxLower_1")
        call test%run(test_getBandPhotonFromEnergyFluence_6, SK_"test_getBandPhotonFromEnergyFluence_6")
        call test%run(test_getBandPhotonFromEnergyFluence_7, SK_"test_getBandPhotonFromEnergyFluence_7")
        call test%run(test_getBandPhotonFromEnergyFluence_8, SK_"test_getBandPhotonFromEnergyFluence_8")
        call test%run(test_getBandPhotonFromEnergyFluence_1, SK_"test_getBandPhotonFromEnergyFluence_1") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getBandPhotonFromEnergyFluence_2, SK_"test_getBandPhotonFromEnergyFluence_2") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getBandPhotonFromEnergyFluence_3, SK_"test_getBandPhotonFromEnergyFluence_3") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getBandPhotonFromEnergyFluence_4, SK_"test_getBandPhotonFromEnergyFluence_4") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%run(test_getBandPhotonFromEnergyFluence_5, SK_"test_getBandPhotonFromEnergyFluence_5") ! The internal function passing as actual argument causes segfault with Gfortran (any version) on Windows subsystem for Linux.
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandEbreak() result(assertion)
        use pm_kind, only: RK, IK
        implicit none
        logical(LK)  :: assertion
        real(RK) :: ebrk, difference
        ebrk = getBandEbreak(BAND_SPEC1%epk,BAND_SPEC1%alpha,BAND_SPEC1%beta)
        difference = 2._RK * abs(ebrk - BAND_SPEC1%ebrk) / (ebrk + BAND_SPEC1%ebrk)
        assertion = difference < 1.e-7_RK
        ! LCOV_EXCL_START
        if (test%traceable .and. .not. assertion) then
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "Ebreak, Reference Ebreak, difference"
            write(test%disp%unit,"(*(g0,:,', '))") ebrk, BAND_SPEC1%ebrk, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandEbreak

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandParam_1() result(assertion)

        use pm_kind, only: RK, IK
        implicit none

        real(RK), parameter :: alphaPlusTwo_ref = 0.5_RK
        real(RK), parameter :: tolerance = 1.e-10_RK
        real(RK), parameter :: ebrk_ref = 6.e2_RK
        real(RK), parameter :: coef_ref = 220.72766470286541_RK
        real(RK), parameter :: alpha = -1.5_RK
        real(RK), parameter :: beta = -2.5_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical(LK)  :: assertion
        real(RK) :: difference
        real(RK) :: alphaPlusTwo
        real(RK) :: ebrk
        real(RK) :: coef

        call getBandParam(epk,alpha,beta,ebrk,coef,alphaPlusTwo)

        difference = abs( (ebrk - ebrk_ref) / ebrk_ref )
        assertion = difference < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "ebrk_ref   ", ebrk_ref
            write(test%disp%unit,"(*(g0,:,', '))") "ebrk       ", ebrk
            write(test%disp%unit,"(*(g0,:,', '))") "difference ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        difference = abs( (coef - coef_ref) / coef_ref )
        assertion = difference < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "coef_ref   ", coef_ref
            write(test%disp%unit,"(*(g0,:,', '))") "coef       ", coef
            write(test%disp%unit,"(*(g0,:,', '))") "difference ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

        difference = abs( (alphaPlusTwo - alphaPlusTwo_ref) / alphaPlusTwo_ref )
        assertion = difference < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "alphaPlusTwo_ref   ", alphaPlusTwo_ref
            write(test%disp%unit,"(*(g0,:,', '))") "alphaPlusTwo       ", alphaPlusTwo
            write(test%disp%unit,"(*(g0,:,', '))") "difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getBandParam_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFluxLower_1() result(assertion)

        use pm_kind, only: RK, IK
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

        logical(LK)  :: assertion
        real(RK) :: difference
        real(RK) :: photonFluxLower

        photonFluxLower = getPhotonFluxLower(energy,alpha,alphaPlusTwoOverEpk)
        assertion = photonFluxLower > 0._RK

        difference = abs( (photonFluxLower - photonFluxLower_ref) / photonFluxLower_ref )
        assertion = assertion .and. difference < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFluxLower_ref", photonFluxLower_ref
            write(test%disp%unit,"(*(g0,:,', '))") "photonFluxLower    ", photonFluxLower
            write(test%disp%unit,"(*(g0,:,', '))") "difference         ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFluxLower_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFlux_1() result(assertion)

        use pm_kind, only: RK, IK
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

        logical(LK)  :: assertion
        real(RK) :: difference
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux > 0._RK

        difference = abs( (photonFlux - photonFlux_ref) / photonFlux_ref )
        assertion = assertion .and. difference < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFlux_ref ", photonFlux_ref
            write(test%disp%unit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getPhotonFlux_2() result(assertion)

        use pm_kind, only: RK, IK
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

        logical(LK)  :: assertion
        real(RK) :: difference
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux > 0._RK

        difference = abs( (photonFlux - photonFlux_ref) / photonFlux_ref )
        assertion = assertion .and. difference < tolerance

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFlux_ref ", photonFlux_ref
            write(test%disp%unit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(test%disp%unit,"(*(g0,:,', '))") "difference     ", difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! alpha < -2
    function test_getPhotonFlux_3() result(assertion)

        use pm_kind, only: RK, IK
        implicit none

        real(RK), parameter :: alphaPlusTwo = 0.5_RK
        real(RK), parameter :: energy = 1.e1_RK
        real(RK), parameter :: alpha = -3.5_RK
        real(RK), parameter :: beta = -2.5_RK
        real(RK), parameter :: coef = 220.72766470286541_RK
        real(RK), parameter :: ebrk = 6.e2_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical(LK)  :: assertion
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux < 0._RK

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! alpha < beta
    function test_getPhotonFlux_4() result(assertion)

        use pm_kind, only: RK, IK
        implicit none

        real(RK), parameter :: alphaPlusTwo = 0.5_RK
        real(RK), parameter :: energy = 1.e1_RK
        real(RK), parameter :: alpha = -1.5_RK
        real(RK), parameter :: beta = -0.5_RK
        real(RK), parameter :: coef = 220.72766470286541_RK
        real(RK), parameter :: ebrk = 6.e2_RK
        real(RK), parameter :: epk = 3.e2_RK

        logical(LK)  :: assertion
        real(RK) :: photonFlux

        photonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
        assertion = photonFlux < 0._RK

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFlux     ", photonFlux
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP

    end function test_getPhotonFlux_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of both upper and lower tails.
    function test_getBandPhoton_3() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhoton   ( lowerLim      = BAND_SPEC3%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhoton_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getBandPhoton_5() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: tolerance = 1.e-10_RK
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhoton   ( lowerLim      = BAND_SPEC3%Limit(2)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhoton_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getBandPhoton_6() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence
        type(err_type)  :: Err
        call getBandPhoton   ( lowerLim      = BAND_SPEC3%Limit(2)   &
                                , upperLim      = BAND_SPEC3%Limit(1)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        assertion = abs(photonFluence - 0._RK) < 1.e-12_RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence", photonFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhoton_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test with conflicting alpha photon index.
    function test_getBandPhoton_7() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence
        type(err_type)  :: Err
        call getBandPhoton   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = -1.e1_RK              &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        assertion = err%occurred .and. photonFluence < 0._RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence", photonFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhoton_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test with conflicting alpha < beta photon indices.
    function test_getBandPhoton_8() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhoton   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%beta       &
                                , beta          = BAND_SPEC3%alpha      &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , photonFluence = photonFluence         &
                                , Err           = Err                   &
                                )
        assertion = err%occurred .and. photonFluence < 0._RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhoton_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of both upper and lower tails.
    function test_getFluenceEnergy_3() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: energyFluence, difference
        type(err_type)  :: Err
        call getFluenceEnergy   ( lowerLim      = BAND_SPEC3%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") energyFluence, BAND_SPEC3%energyFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getFluenceEnergy_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getFluenceEnergy_6() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: EnergyFluence
        type(err_type)  :: Err
        call getFluenceEnergy   ( lowerLim      = BAND_SPEC3%Limit(2)   &
                                , upperLim      = BAND_SPEC3%Limit(1)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%alpha      &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , EnergyFluence = EnergyFluence         &
                                , Err           = Err                   &
                                )
        assertion = abs(EnergyFluence - 0._RK) < 1.e-12_RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "EnergyFluence", EnergyFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getFluenceEnergy_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test with conflicting alpha photon index.
    function test_getFluenceEnergy_7() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: EnergyFluence
        type(err_type)  :: Err
        call getFluenceEnergy   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = -1.e1_RK              &
                                , beta          = BAND_SPEC3%beta       &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , EnergyFluence = EnergyFluence         &
                                , Err           = Err                   &
                                )
        assertion = err%occurred .and. EnergyFluence < 0._RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "EnergyFluence", EnergyFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getFluenceEnergy_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test with conflicting alpha < beta photon indices.
    function test_getFluenceEnergy_8() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: EnergyFluence
        type(err_type)  :: Err
        call getFluenceEnergy   ( lowerLim      = BAND_SPEC3%Limit(1)   &
                                , upperLim      = BAND_SPEC3%Limit(2)   &
                                , epk           = BAND_SPEC3%epk        &
                                , alpha         = BAND_SPEC3%beta       &
                                , beta          = BAND_SPEC3%alpha      &
                                , tolerance     = BAND_SPEC3%tolerance  &
                                , EnergyFluence = EnergyFluence         &
                                , Err           = Err                   &
                                )
        assertion = err%occurred .and. EnergyFluence < 0._RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "EnergyFluence", EnergyFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getFluenceEnergy_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of when lower limit is larger than upper limit.
    function test_getBandPhotonFromEnergyFluence_6() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFluence", photonFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhotonFromEnergyFluence_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test with conflicting alpha photon index alpha < -2.
    function test_getBandPhotonFromEnergyFluence_7() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
                                                , lowerLim      = BAND_SPEC2%Limit(1)       &
                                                , upperLim      = BAND_SPEC2%Limit(2)       &
                                                , epk           = BAND_SPEC2%epk            &
                                                , alpha         = -1.e1_RK                  &
                                                , beta          = BAND_SPEC2%beta           &
                                                , tolerance     = BAND_SPEC2%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        assertion = err%occurred .and. photonFluence < 0._RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFluence", photonFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhotonFromEnergyFluence_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test with conflicting alpha < beta photon indices.
    function test_getBandPhotonFromEnergyFluence_8() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
                                                , lowerLim      = BAND_SPEC2%Limit(1)       &
                                                , upperLim      = BAND_SPEC2%Limit(2)       &
                                                , epk           = BAND_SPEC2%epk            &
                                                , alpha         = BAND_SPEC2%beta           &
                                                , beta          = BAND_SPEC2%alpha          &
                                                , tolerance     = BAND_SPEC2%tolerance      &
                                                , photonFluence = photonFluence             &
                                                , Err           = Err                       &
                                                )
        assertion = err%occurred .and. photonFluence < 0._RK
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photonFluence", photonFluence
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
     end function test_getBandPhotonFromEnergyFluence_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if !WSL_ENABLED || !CODECOV_ENABLED || DLL_ENABLED

    !>  \brief
    !> Test the integration of both the upper and lower tails.
    function test_getBandPhoton_1() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhoton   ( lowerLim      = BAND_SPEC1%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC1%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhoton_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of only the upper tail.
    function test_getBandPhoton_2() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhoton   ( lowerLim      = BAND_SPEC2%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhoton_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of both upper and upper tails with steep slopes.
    function test_getBandPhoton_4() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err

        call getBandPhoton   ( lowerLim      = BAND_SPEC4%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC4%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhoton_4

    !>  \brief
    !> Test the integration of both upper and upper tails.
    function test_getFluenceEnergy_1() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: energyFluence, difference
        type(err_type)  :: Err
        call getFluenceEnergy   ( lowerLim      = BAND_SPEC1%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") energyFluence, BAND_SPEC1%energyFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getFluenceEnergy_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of only the upper tail.
    function test_getFluenceEnergy_2() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: energyFluence, difference
        type(err_type)  :: Err
        call getFluenceEnergy   ( lowerLim      = BAND_SPEC2%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") energyFluence, BAND_SPEC2%energyFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getFluenceEnergy_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the integration of both upper and upper tails with steep slopes.
    function test_getFluenceEnergy_4() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: energyFluence, difference
        type(err_type)  :: Err
        call getFluenceEnergy   ( lowerLim      = BAND_SPEC4%Limit(1)   &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") energyFluence, BAND_SPEC4%energyFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getFluenceEnergy_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandPhotonFromEnergyFluence_1() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC1%energyFluence  &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC1%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhotonFromEnergyFluence_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandPhotonFromEnergyFluence_2() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC2%energyFluence  &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhotonFromEnergyFluence_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandPhotonFromEnergyFluence_3() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC3%energyFluence  &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC3%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhotonFromEnergyFluence_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandPhotonFromEnergyFluence_4() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC4%energyFluence  &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC4%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhotonFromEnergyFluence_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getBandPhotonFromEnergyFluence_5() result(assertion)
        use pm_kind, only: RK, IK
        use pm_err, only: err_type
        implicit none
        logical(LK)     :: assertion
        real(RK)        :: photonFluence, difference
        type(err_type)  :: Err
        call getBandPhotonFromEnergyFluence  ( energyFluence = BAND_SPEC1%energyFluence  &
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
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0,:,', '))")
            write(test%disp%unit,"(*(g0,:,', '))") "photon fluence, Reference photon fluence, difference"
            write(test%disp%unit,"(*(g0,:,', '))") photonFluence, BAND_SPEC2%photonFluence, difference
            write(test%disp%unit,"(*(g0,:,', '))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBandPhotonFromEnergyFluence_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#endif

end module test_pm_distBand ! LCOV_EXCL_LINE