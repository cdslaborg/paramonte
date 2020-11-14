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

module BandSpectrum_mod
! This module contains methods for computing the Band Spectrum properties of GRBs.
! The Band spectral form is given by Eqn A6 of Shahmoradi and Nemiroff 2015: 
! "Short versus long gamma-ray bursts: a comprehensive study of energetics and prompt gamma-ray correlations"
! Amir Shahmoradi, Tuesday April 30, 2019, 12:58 PM, SEIR, UTA

    use Constants_mod, only: IK, RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@BandSpectrum_mod"

    real(RK), parameter     :: AMPLITUDE_DEFAULT = 1._RK    ! default amplitude of the Band function
    real(RK), parameter     :: ALPHA_DEFAULT = -1.1_RK      ! default low-energy index of the Band function
    real(RK), parameter     :: BETA_DEFAULT = -2.3_RK       ! default high-energy index of the Band function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns epk*(alpha-beta)/(alpha+2), which is the normalization factor in the exponent of the first component of the Band model
    pure function getEbreak(epk,alpha,beta) result(ebrk)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEbreak
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)            :: epk
        real(RK), intent(in), optional  :: alpha, beta
        real(RK)                        :: ebrk
        ebrk = epk * (alpha-beta) / (alpha+2._RK)
    end function getEbreak

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the break energy ebrk o fthe Band model, as well as the coefficient by which the high energy component of the
    ! Band model must be multiplied in order to get a smooth function (assuming the coefficient of the lower energy component
    ! is unity).
    pure subroutine getBandParam(epk,alpha,beta,ebrk,coef,alphaPlusTwo)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBandParam
#endif
        use Constants_mod, only: IK, RK
        implicit none
        real(RK), intent(in)    :: epk,alpha,beta
        real(RK), intent(out)   :: ebrk, coef, alphaPlusTwo
        real(RK)                :: alphaMinusBeta
        alphaPlusTwo = alpha + 2._RK
        alphaMinusBeta = alpha - beta
        ebrk = epk * alphaMinusBeta / alphaPlusTwo
        coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta)
    end subroutine getBandParam

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the value the Band differential spectrum for the given energy and input parameters.
    ! It is expected that the input energy is in units of KeV, although it does not affect the computations here.
    ! NOTE: A negative huge output value is used to signal error has occurred. Under normal conditions, output is always positive.
    pure function getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo) result(photonFlux)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getPhotonFlux
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: energy, epk, ebrk, alpha, beta,coef,alphaPlusTwo
        real(RK)                        :: photonFlux
        ! check if the photon indices are consistent with the mathematical rules
        if (alpha<beta .or. alpha<-2._RK) then
            photonFlux = -HUGE_RK
            return
        end if
        ! compute the spectrum
        if (energy<=ebrk) then
            photonFlux = energy**alpha * exp(-energy*alphaPlusTwo/epk)
        else
            photonFlux = coef * energy**beta
        end if
    end function getPhotonFlux

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the differential photon count as a function of energy in keV, given the lower component of the Band model.
    ! NOTE: energy values beyond ebrk should not be passed to this function.
    pure function getPhotonFluxLower(energy,alpha,alphaPlusTwoOverEpk) result(photonFluxLower)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getPhotonFluxLower
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: energy, alpha, alphaPlusTwoOverEpk
        real(RK)                        :: photonFluxLower
        photonFluxLower = energy**alpha * exp(-energy*alphaPlusTwoOverEpk)  ! the lower-energy spectrum
    end function getPhotonFluxLower

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the integral of the Band differential spectrum over the input energy range.
    ! It is expected that the input energy is in units of KeV, although it does not affect the computations here.
    ! NOTE: A negative huge output value is used to signal error has occurred. Under normal conditions, output is always positive.
    subroutine getPhotonFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,photonFluence,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getPhotonFluence
#endif

       !use Integration_mod, only: doQuadRombClosed
        use Constants_mod, only: IK, RK, HUGE_RK
        use QuadPackSPR_mod, only: qag
        use Err_mod, only: Err_type
        implicit none
        real(RK), intent(in)            :: lowerLim, upperLim, epk, alpha, beta, tolerance
        real(RK), intent(out)           :: photonFluence
        type(Err_type), intent(out)     :: Err
        character(*), parameter         :: PROCEDURE_NAME = "@getPhotonFluence()"
        real(RK)                        :: ebrk, alphaPlusTwo
        real(RK)                        :: thisUpperLim, alphaPlusTwoOverEpk, betaPlusOne
        real(RK)                        :: alphaMinusBeta, coef
        real(RK)                        :: abserr
        integer(IK)                     :: neval
        integer(IK)                     :: ierr
       !real(RK)                        :: alphaPlusOne, logGammaAlphaPlusOne

        Err%occurred = .false.

        if (lowerLim>=upperLim) then
            photonFluence = 0._RK
            return
        end if

        ! check if the photon indices are consistent with the mathematical rules
        if (alpha<beta .or. alpha<-2._RK) then
            photonFluence = -HUGE_RK
            Err%occurred = .true.
            Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred: alpha<beta .or. alpha<-2._RK"
            return
        end if

        ! integrate the spectrum
        alphaPlusTwo = alpha + 2._RK
        alphaMinusBeta = alpha - beta
        ebrk = epk*alphaMinusBeta/alphaPlusTwo

        if (lowerLim>ebrk) then

            ! there is only the high energy component in the photonFluence
            betaPlusOne = beta + 1._RK
            coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
            photonFluence = coef * ( upperLim**betaPlusOne - lowerLim**betaPlusOne ) / betaPlusOne
            return

        elseif (lowerLim<ebrk) then

            alphaPlusTwoOverEpk = alphaPlusTwo / epk
            thisUpperLim = min(upperLim,ebrk)
            !alphaPlusOne = alpha + 1._RK
            !if (alpha>-1._RK) then
            !    logGammaAlphaPlusOne = log_gamma( alphaPlusOne )
            !    ! use the analytical approach to compute the photonFluence:
            !    ! https://www.wolframalpha.com/input/?i=integrate+x%5Ea+*+exp(-b*x)
            !    photonFluence = getUpperGamma( exponent = alphaPlusOne &
            !                            , logGammaExponent = logGammaAlphaPlusOne &
            !                            , lowerLim = alphaPlusTwoOverEpk * lowerLim &
            !                            , tolerance = tolerance &
            !                            ) &
            !             - getUpperGamma( exponent = alphaPlusOne &
            !                            , logGammaExponent = logGammaAlphaPlusOne &
            !                            , lowerLim = alphaPlusTwoOverEpk * thisUpperLim &
            !                            , tolerance = tolerance &
            !                            )
            !    photonFluence = photonFluence / alphaPlusTwoOverEpk**alphaPlusOne
            !else
                ! use brute-force integration
                call qag( f             = getBandCompLowPhoton  &
                        , a             = lowerLim              &
                        , b             = thisUpperLim          &
                        , epsabs        = 0._RK                 &
                        , epsrel        = tolerance             &
                        , key           = 1_IK                  &
                        , result        = photonFluence         &
                        , abserr        = abserr                &
                        , neval         = neval                 &
                        , ier           = ierr                  &
                        )
                !write(*,*) neval
                !call doQuadRombClosed   ( getFunc       = getBandCompLowPhoton &
                !                        , xmin          = lowerLim             &
                !                        , xmax          = thisUpperLim         &
                !                        , tolerance     = 1.e-7_RK             &
                !                        , nRefinement   = 10_IK                &
                !                        , photonFluence = photonFluence        &
                !                        , ierr          = ierr                 &
                !                        )
                if (ierr/=0_IK) then
                    photonFluence = -HUGE_RK
                    Err%occurred = .true.
                    Err%stat = ierr
                    Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred at QuadPack routine. Check the error code to identify the root cause."
                    return
                end if
            !end if

            if (upperLim>ebrk) then
                ! add the remaining part of the photonFluence from the high-energy component
                betaPlusOne = beta + 1._RK
                alphaMinusBeta = alpha - beta
                coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
                photonFluence = photonFluence + coef * ( upperLim**betaPlusOne - ebrk**betaPlusOne ) / betaPlusOne
                return
            end if

        end if

    contains

        pure function getBandCompLowPhoton(energy) result(bandCompLow)
            implicit none
            real(RK), intent(in)    :: energy
            real(RK)                :: bandCompLow
            bandCompLow = energy**alpha * exp(-alphaPlusTwoOverEpk*energy)
        end function getBandCompLowPhoton

    end subroutine getPhotonFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the integral of the Band differential spectrum over the input energy range, in units of the input energy.
    ! It is expected that the input energy is in units of KeV, although it does not affect the computations here.
    ! NOTE: A negative huge output value is used to signal error has occurred. Under normal conditions, output is always positive.
    subroutine getEnergyFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,energyFluence,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEnergyFluence
#endif

        use Constants_mod, only: IK, RK, HUGE_RK
        use QuadPackSPR_mod, only: qag
        use Err_mod, only: Err_type
        implicit none
        real(RK), intent(in)            :: lowerLim, upperLim, epk, alpha, beta, tolerance
        type(Err_type), intent(out)     :: Err
        real(RK), intent(out)           :: energyFluence
        character(*), parameter         :: PROCEDURE_NAME = "@getEnergyFluence()"
        real(RK)                        :: ebrk, alphaPlusTwo
        real(RK)                        :: thisUpperLim, alphaPlusTwoOverEpk, betaPlusTwo
        real(RK)                        :: alphaMinusBeta, coef, alphaPlusOne
        real(RK)                        :: abserr
        integer(IK)                     :: neval
        integer(IK)                     :: ierr

        Err%occurred = .false.

        if (lowerLim>=upperLim) then
            energyFluence = 0._RK
            return
        end if

        ! check if the photon indices are consistent with the mathematical rules
        if (alpha<beta .or. alpha<-2._RK) then
            energyFluence = -HUGE_RK
            Err%occurred = .true.
            Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred: alpha<beta .or. alpha<-2._RK"
            return
        end if

        ! integrate the spectrum
        alphaPlusTwo = alpha + 2._RK
        alphaMinusBeta = alpha - beta
        ebrk = epk*alphaMinusBeta/alphaPlusTwo

        if (lowerLim>ebrk) then

            ! there is only the high energy component in the energyFluence
            betaPlusTwo = beta + 2._RK
            coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
            energyFluence = coef * ( upperLim**betaPlusTwo - lowerLim**betaPlusTwo ) / betaPlusTwo
            return

        elseif (lowerLim<ebrk) then

            alphaPlusTwoOverEpk = alphaPlusTwo / epk
            thisUpperLim = min(upperLim,ebrk)
            alphaPlusOne = alpha + 1._RK
            ! use brute-force integration
            call qag( f             = getBandCompLowEnergy  &
                    , a             = lowerLim              &
                    , b             = thisUpperLim          &
                    , epsabs        = 0._RK                 &
                    , epsrel        = tolerance             &
                    , key           = 1_IK                  &
                    , result        = energyFluence         &
                    , abserr        = abserr                &
                    , neval         = neval                 &
                    , ier           = ierr                  &
                    )

            if (ierr/=0_IK) then
                energyFluence = -HUGE_RK
                Err%occurred = .true.
                Err%stat = ierr
                Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred at QuadPack routine. Check the error code to identify the root cause."
                return
            end if

            if (upperLim>ebrk) then ! add the remaining part of the energyFluence from the high-energy component
                betaPlusTwo = beta + 2._RK
                alphaMinusBeta = alpha - beta
                coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
                energyFluence = energyFluence + coef * ( upperLim**betaPlusTwo - ebrk**betaPlusTwo ) / betaPlusTwo
                return
            end if

        end if

    contains

        function getBandCompLowEnergy(energy) result(bandCompLow)
            implicit none
            real(RK), intent(in)    :: energy
            real(RK)                :: bandCompLow
            bandCompLow = energy**alphaPlusOne * exp(-alphaPlusTwoOverEpk*energy)
        end function getBandCompLowEnergy

    end subroutine getEnergyFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! converts an input energy fluence in [lowerLim, upperLim] energy window, all in units of keV (or ALL in some other units),
    ! to photon fluence, within the same energy range, or if [lowerLimNew, upperLimNew] is provided, then in that range.
    ! tolerance is passed to integration units as a measure of the desired accuracy for the computation of the fluence.
    ! NOTE: A negative huge output value is used to signal error has occurred. Under normal conditions, output is always positive.
    subroutine getPhotonFluenceFromEnergyFluence( energyFluence, lowerLim, upperLim, epk, alpha, beta, tolerance &
                                                , photonFluence, Err, lowerLimNew, upperLimNew )

#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getPhotonFluenceFromEnergyFluence
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        use Err_mod, only: Err_type
        implicit none
        real(RK), intent(in)            :: energyFluence, lowerLim, upperLim, epk, alpha, beta, tolerance
        real(RK), intent(out)           :: photonFluence
        type(Err_type), intent(out)     :: Err
        character(*), parameter         :: PROCEDURE_NAME = "@getPhotonFluenceFromEnergyFluence()"
        real(RK), intent(in), optional  :: lowerLimNew, upperLimNew
        real(RK)                        :: lowLimNew, uppLimNew, amplitude

        Err%occurred = .false.

        ! check if the photon indices are consistent with the mathematical rules
        if (lowerLim>=upperLim .or. alpha<beta .or. alpha<-2._RK) then
            Err%occurred = .true.
            Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred: lowerLim>=upperLim .or. alpha<beta .or. alpha<-2._RK"
            photonFluence = -HUGE_RK
            return
        end if

        lowLimNew = lowerLim
        if (present(lowerLimNew)) lowLimNew = lowerLimNew
        uppLimNew = upperLim
        if (present(upperLimNew)) uppLimNew = upperLimNew

        call getEnergyFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,amplitude,Err)
        if (Err%occurred) then
            photonFluence = -HUGE_RK
            Err%msg = MODULE_NAME // PROCEDURE_NAME // Err%msg
            return
        end if
        amplitude = energyFluence / amplitude

        call getPhotonFluence(lowLimNew,uppLimNew,epk,alpha,beta,tolerance,photonFluence,Err)
        if (Err%occurred) then
            photonFluence = -HUGE_RK
            Err%msg = MODULE_NAME // PROCEDURE_NAME // Err%msg
            return
        end if
        photonFluence = amplitude * photonFluence

    end subroutine getPhotonFluenceFromEnergyFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module BandSpectrum_mod