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

!> \brief
!> This module contains the procedures for computing the various Band Spectrum properties of GRBs.
!> See Eqn. A6 of [Shahmoradi and Nemiroff, 2015, MNRAS, 451:4645-4662](https://www.cdslab.org/pubs/2015_Shahmoradi_I.pdf).
!>
!> \author
!> Amir Shahmoradi, Tuesday April 30, 2019, 12:58 PM, SEIR, UTA

module BandSpectrum_mod

    use Constants_mod, only: IK, RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@BandSpectrum_mod"

    real(RK), parameter     :: AMPLITUDE_DEFAULT = 1._RK    ! default amplitude of the Band function
    real(RK), parameter     :: ALPHA_DEFAULT = -1.1_RK      ! default low-energy index of the Band function
    real(RK), parameter     :: BETA_DEFAULT = -2.3_RK       ! default high-energy index of the Band function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module subroutine getEnergyFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,energyFluence,Err)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEnergyFluence
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        use QuadPackSPR_mod, only: qag
        use Err_mod, only: Err_type
        implicit none
        real(RK), intent(in)            :: lowerLim, upperLim, epk, alpha, beta, tolerance
        real(RK), intent(out)           :: energyFluence
        type(Err_type), intent(out)     :: Err
    end subroutine getEnergyFluence
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module subroutine getPhotonFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,photonFluence,Err)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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
    end subroutine getPhotonFluence
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return `epk*(alpha-beta)/(alpha+2)`, which is the normalization
    !> factor in the exponent of the first component of the Band model.
    !> \param[in]   epk     :   The spectral peak energy in units of [keV].
    !> \param[in]   alpha   :   The lower spectral exponent of the Band model.
    !> \param[in]   beta    :   The upper spectral exponent of the Band model.
    !>
    !> \return
    !> `ebrk` : The spectral break energy in units of [keV].
    pure function getEbreak(epk,alpha,beta) result(ebrk)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> Compute the break energy `ebrk` of the Band model, as well as the coefficient by which the high energy component of the
    !> Band model must be multiplied in order to get a smooth function (assuming the normalization coefficient of the lower
    !> energy component is unity).
    !>
    !> \param[in]   epk             :   The spectral peak energy in units of [keV].
    !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
    !> \param[in]   beta            :   The upper spectral exponent of the Band model.
    !> \param[out]  ebrk            :   The spectral break energy in units of [keV].
    !> \param[out]  coef            :   The spectral continuity coefficient.
    !> \param[out]  alphaPlusTwo    :   The lower spectral exponent of the Band model plus two.
    pure subroutine getBandParam(epk,alpha,beta,ebrk,coef,alphaPlusTwo)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> Compute the differential photon flux according to the Band differential spectrum at the given input `energy` value.
    !>
    !> \param[in]   energy          :   The energy (in units of [keV]) at which the differential photon flux must be computed.
    !> \param[in]   epk             :   The spectral peak energy in units of [keV].
    !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
    !> \param[in]   beta            :   The upper spectral exponent of the Band model.
    !> \param[in]   ebrk            :   The spectral break energy in units of [keV].
    !> \param[in]   coef            :   The spectral continuity coefficient.
    !> \param[in]   alphaPlusTwo    :   The lower spectral exponent of the Band model plus two.
    !>
    !> \return
    !> `photonFlux` : The energy flux in units of photon counts.
    !>
    !> \warning
    !> A negative huge output value is used to signal error has occurred. Under normal conditions, the output is always positive.
    !>
    !> \warning
    !> The input energy values `energy`, `epk`, `ebrk`, must all have the same units.
    !> It is expected that the input energy is in units of keV, although it does not affect the accuracy of the results.
    pure function getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo) result(photonFlux)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> Compute the differential photon flux according to the lower component of the Band spectrum at the given input `energy` value.
    !>
    !> \param[in]   energy              :   The energy (in units of [keV]) at which the differential photon flux must be computed.
    !> \param[in]   alpha               :   The lower spectral exponent of the Band model.
    !> \param[in]   alphaPlusTwoOverEpk :   The lower spectral exponent of the Band model plus two.
    !>
    !> \return
    !> `photonFluxLower` : The energy flux in units of photon counts.
    !>
    !> \warning
    !> The input `energy` values must be less than `ebrk`.
    pure function getPhotonFluxLower(energy,alpha,alphaPlusTwoOverEpk) result(photonFluxLower)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getPhotonFluxLower
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: energy, alpha, alphaPlusTwoOverEpk
        real(RK)                        :: photonFluxLower
        photonFluxLower = energy**alpha * exp(-energy*alphaPlusTwoOverEpk)  ! the lower-energy spectrum
    end function getPhotonFluxLower

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Bizarrely and frustratingly, Microsoft Windows Subsystem for Linux Ubuntu with GFortran yields Segmentation faults with internal procedure calls.
! This is so unfortunate. To bypass this issue for now, the following subroutine is implemented as separate submodule 
! so that the internal shared parameters can be safely passed as submodule parameters.

!WSL_GFORTRAN_BUG     !> \brief
!WSL_GFORTRAN_BUG     !> Integrate the Band differential spectrum over the input energy range.
!WSL_GFORTRAN_BUG     !>
!WSL_GFORTRAN_BUG     !> \param[in]   lowerLim        :   The lower limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   upperLim        :   The upper limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   epk             :   The spectral peak energy in units of [keV].
!WSL_GFORTRAN_BUG     !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   beta            :   The upper spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   tolerance       :   The relative accuracy tolerance of the integration.
!WSL_GFORTRAN_BUG     !> \param[out]  photonFluence   :   The fluence in units of photon counts within the input energy range.
!WSL_GFORTRAN_BUG     !> \param[out]  Err             :   An object of class [Err_type](@ref err_mod::err_type) containing error-handling information.
!WSL_GFORTRAN_BUG     subroutine getPhotonFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,photonFluence,Err)
!WSL_GFORTRAN_BUG #if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
!WSL_GFORTRAN_BUG         !DEC$ ATTRIBUTES DLLEXPORT :: getPhotonFluence
!WSL_GFORTRAN_BUG #endif
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG        !use Integration_mod, only: doQuadRombClosed
!WSL_GFORTRAN_BUG         use Constants_mod, only: IK, RK, HUGE_RK
!WSL_GFORTRAN_BUG         use QuadPackSPR_mod, only: qag
!WSL_GFORTRAN_BUG         use Err_mod, only: Err_type
!WSL_GFORTRAN_BUG         implicit none
!WSL_GFORTRAN_BUG         real(RK), intent(in)            :: lowerLim, upperLim, epk, alpha, beta, tolerance
!WSL_GFORTRAN_BUG         real(RK), intent(out)           :: photonFluence
!WSL_GFORTRAN_BUG         type(Err_type), intent(out)     :: Err
!WSL_GFORTRAN_BUG         character(*), parameter         :: PROCEDURE_NAME = "@getPhotonFluence()"
!WSL_GFORTRAN_BUG         real(RK)                        :: ebrk, alphaPlusTwo
!WSL_GFORTRAN_BUG         real(RK)                        :: thisUpperLim, alphaPlusTwoOverEpk, betaPlusOne
!WSL_GFORTRAN_BUG         real(RK)                        :: alphaMinusBeta, coef
!WSL_GFORTRAN_BUG         real(RK)                        :: abserr
!WSL_GFORTRAN_BUG         integer(IK)                     :: neval
!WSL_GFORTRAN_BUG         integer(IK)                     :: ierr
!WSL_GFORTRAN_BUG        !real(RK)                        :: alphaPlusOne, logGammaAlphaPlusOne
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         Err%occurred = .false.
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         if (lowerLim>=upperLim) then
!WSL_GFORTRAN_BUG             photonFluence = 0._RK
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         ! check if the photon indices are consistent with the mathematical rules
!WSL_GFORTRAN_BUG         if (alpha<beta .or. alpha<-2._RK) then
!WSL_GFORTRAN_BUG             photonFluence = -HUGE_RK
!WSL_GFORTRAN_BUG             Err%occurred = .true.
!WSL_GFORTRAN_BUG             Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred: alpha<beta .or. alpha<-2._RK"
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         ! integrate the spectrum
!WSL_GFORTRAN_BUG         alphaPlusTwo = alpha + 2._RK
!WSL_GFORTRAN_BUG         alphaMinusBeta = alpha - beta
!WSL_GFORTRAN_BUG         ebrk = epk*alphaMinusBeta/alphaPlusTwo
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         if (lowerLim>ebrk) then
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG             ! there is only the high energy component in the photonFluence
!WSL_GFORTRAN_BUG             betaPlusOne = beta + 1._RK
!WSL_GFORTRAN_BUG             coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
!WSL_GFORTRAN_BUG             photonFluence = coef * ( upperLim**betaPlusOne - lowerLim**betaPlusOne ) / betaPlusOne
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG !#if defined OS_IS_WSL && defined CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_START
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         elseif (lowerLim<ebrk) then
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG             alphaPlusTwoOverEpk = alphaPlusTwo / epk
!WSL_GFORTRAN_BUG             thisUpperLim = min(upperLim,ebrk)
!WSL_GFORTRAN_BUG             !alphaPlusOne = alpha + 1._RK
!WSL_GFORTRAN_BUG             !if (alpha>-1._RK) then
!WSL_GFORTRAN_BUG             !    logGammaAlphaPlusOne = log_gamma( alphaPlusOne )
!WSL_GFORTRAN_BUG             !    ! use the analytical approach to compute the photonFluence:
!WSL_GFORTRAN_BUG             !    ! https://www.wolframalpha.com/input/?i=integrate+x%5Ea+*+exp(-b*x)
!WSL_GFORTRAN_BUG             !    photonFluence = getUpperGamma( exponent = alphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , logGammaExponent = logGammaAlphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , lowerLim = alphaPlusTwoOverEpk * lowerLim &
!WSL_GFORTRAN_BUG             !                            , tolerance = tolerance &
!WSL_GFORTRAN_BUG             !                            ) &
!WSL_GFORTRAN_BUG             !             - getUpperGamma( exponent = alphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , logGammaExponent = logGammaAlphaPlusOne &
!WSL_GFORTRAN_BUG             !                            , lowerLim = alphaPlusTwoOverEpk * thisUpperLim &
!WSL_GFORTRAN_BUG             !                            , tolerance = tolerance &
!WSL_GFORTRAN_BUG             !                            )
!WSL_GFORTRAN_BUG             !    photonFluence = photonFluence / alphaPlusTwoOverEpk**alphaPlusOne
!WSL_GFORTRAN_BUG             !else
!WSL_GFORTRAN_BUG                 ! use brute-force integration
!WSL_GFORTRAN_BUG                 call qag( f             = getBandCompLowPhoton  &
!WSL_GFORTRAN_BUG                         , a             = lowerLim              &
!WSL_GFORTRAN_BUG                         , b             = thisUpperLim          &
!WSL_GFORTRAN_BUG                         , epsabs        = 0._RK                 &
!WSL_GFORTRAN_BUG                         , epsrel        = tolerance             &
!WSL_GFORTRAN_BUG                         , key           = 1_IK                  &
!WSL_GFORTRAN_BUG                         , result        = photonFluence         &
!WSL_GFORTRAN_BUG                         , abserr        = abserr                &
!WSL_GFORTRAN_BUG                         , neval         = neval                 &
!WSL_GFORTRAN_BUG                         , ier           = ierr                  &
!WSL_GFORTRAN_BUG                         )
!WSL_GFORTRAN_BUG                 !write(*,*) neval
!WSL_GFORTRAN_BUG                 !call doQuadRombClosed   ( getFunc       = getBandCompLowPhoton &
!WSL_GFORTRAN_BUG                 !                        , xmin          = lowerLim             &
!WSL_GFORTRAN_BUG                 !                        , xmax          = thisUpperLim         &
!WSL_GFORTRAN_BUG                 !                        , tolerance     = 1.e-7_RK             &
!WSL_GFORTRAN_BUG                 !                        , nRefinement   = 10_IK                &
!WSL_GFORTRAN_BUG                 !                        , photonFluence = photonFluence        &
!WSL_GFORTRAN_BUG                 !                        , ierr          = ierr                 &
!WSL_GFORTRAN_BUG                 !                        )
!WSL_GFORTRAN_BUG                 if (ierr/=0_IK) then
!WSL_GFORTRAN_BUG                     photonFluence = -HUGE_RK
!WSL_GFORTRAN_BUG                     Err%occurred = .true.
!WSL_GFORTRAN_BUG                     Err%stat = ierr
!WSL_GFORTRAN_BUG                     Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred at QuadPack routine. Check the error code to identify the root cause."
!WSL_GFORTRAN_BUG                     return
!WSL_GFORTRAN_BUG                 end if
!WSL_GFORTRAN_BUG             !end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG             if (upperLim>ebrk) then
!WSL_GFORTRAN_BUG                 ! add the remaining part of the photonFluence from the high-energy component
!WSL_GFORTRAN_BUG                 betaPlusOne = beta + 1._RK
!WSL_GFORTRAN_BUG                 alphaMinusBeta = alpha - beta
!WSL_GFORTRAN_BUG                 coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
!WSL_GFORTRAN_BUG                 photonFluence = photonFluence + coef * ( upperLim**betaPlusOne - ebrk**betaPlusOne ) / betaPlusOne
!WSL_GFORTRAN_BUG                 return
!WSL_GFORTRAN_BUG             end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG     contains
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         pure function getBandCompLowPhoton(energy) result(bandCompLow)
!WSL_GFORTRAN_BUG             implicit none
!WSL_GFORTRAN_BUG             real(RK), intent(in)    :: energy
!WSL_GFORTRAN_BUG             real(RK)                :: bandCompLow
!WSL_GFORTRAN_BUG             bandCompLow = energy**alpha * exp(-alphaPlusTwoOverEpk*energy)
!WSL_GFORTRAN_BUG         end function getBandCompLowPhoton
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG !#if defined OS_IS_WSL && defined CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_STOP
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG     end subroutine getPhotonFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Bizarrely and frustratingly, Microsoft Windows Subsystem for Linux Ubuntu with GFortran yields Segmentation faults with internal procedure calls.
! This is so unfortunate. To bypass this issue for now, the following subroutine is implemented as separate submodule 
! so that the internal shared parameters can be safely passed as submodule parameters.

!WSL_GFORTRAN_BUG     !> \brief
!WSL_GFORTRAN_BUG     !> Integrate the Band differential spectrum over the input energy range in units of the input energy.
!WSL_GFORTRAN_BUG     !>
!WSL_GFORTRAN_BUG     !> \param[in]   lowerLim        :   The lower limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   upperLim        :   The upper limit energy (in units of [keV]) of the integration.
!WSL_GFORTRAN_BUG     !> \param[in]   epk             :   The spectral peak energy in units of [keV].
!WSL_GFORTRAN_BUG     !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   beta            :   The upper spectral exponent of the Band model.
!WSL_GFORTRAN_BUG     !> \param[in]   tolerance       :   The relative accuracy tolerance of the integration.
!WSL_GFORTRAN_BUG     !> \param[out]  energyFluence   :   The fluence in units of the input energy (keV) within the input energy range.
!WSL_GFORTRAN_BUG     !> \param[out]  Err             :   An object of class [Err_type](@ref err_mod::err_type) containing error-handling information.
!WSL_GFORTRAN_BUG     subroutine getEnergyFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,energyFluence,Err)
!WSL_GFORTRAN_BUG #if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
!WSL_GFORTRAN_BUG         !DEC$ ATTRIBUTES DLLEXPORT :: getEnergyFluence
!WSL_GFORTRAN_BUG #endif
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         use Constants_mod, only: IK, RK, HUGE_RK
!WSL_GFORTRAN_BUG         use QuadPackSPR_mod, only: qag
!WSL_GFORTRAN_BUG         use Err_mod, only: Err_type
!WSL_GFORTRAN_BUG         implicit none
!WSL_GFORTRAN_BUG         real(RK), intent(in)            :: lowerLim, upperLim, epk, alpha, beta, tolerance
!WSL_GFORTRAN_BUG         type(Err_type), intent(out)     :: Err
!WSL_GFORTRAN_BUG         real(RK), intent(out)           :: energyFluence
!WSL_GFORTRAN_BUG         character(*), parameter         :: PROCEDURE_NAME = "@getEnergyFluence()"
!WSL_GFORTRAN_BUG         real(RK)                        :: ebrk, alphaPlusTwo
!WSL_GFORTRAN_BUG         real(RK)                        :: thisUpperLim, alphaPlusTwoOverEpk, betaPlusTwo
!WSL_GFORTRAN_BUG         real(RK)                        :: alphaMinusBeta, coef, alphaPlusOne
!WSL_GFORTRAN_BUG         real(RK)                        :: abserr
!WSL_GFORTRAN_BUG         integer(IK)                     :: neval
!WSL_GFORTRAN_BUG         integer(IK)                     :: ierr
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         Err%occurred = .false.
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         if (lowerLim>=upperLim) then
!WSL_GFORTRAN_BUG             energyFluence = 0._RK
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         ! check if the photon indices are consistent with the mathematical rules
!WSL_GFORTRAN_BUG         if (alpha<beta .or. alpha<-2._RK) then
!WSL_GFORTRAN_BUG             energyFluence = -HUGE_RK
!WSL_GFORTRAN_BUG             Err%occurred = .true.
!WSL_GFORTRAN_BUG             Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred: alpha<beta .or. alpha<-2._RK"
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         ! integrate the spectrum
!WSL_GFORTRAN_BUG         alphaPlusTwo = alpha + 2._RK
!WSL_GFORTRAN_BUG         alphaMinusBeta = alpha - beta
!WSL_GFORTRAN_BUG         ebrk = epk*alphaMinusBeta/alphaPlusTwo
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         if (lowerLim>ebrk) then
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG             ! there is only the high energy component in the energyFluence
!WSL_GFORTRAN_BUG             betaPlusTwo = beta + 2._RK
!WSL_GFORTRAN_BUG             coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
!WSL_GFORTRAN_BUG             energyFluence = coef * ( upperLim**betaPlusTwo - lowerLim**betaPlusTwo ) / betaPlusTwo
!WSL_GFORTRAN_BUG             return
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG !#if defined OS_IS_WSL && defined CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_START
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         elseif (lowerLim<ebrk) then
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG             alphaPlusTwoOverEpk = alphaPlusTwo / epk
!WSL_GFORTRAN_BUG             thisUpperLim = min(upperLim,ebrk)
!WSL_GFORTRAN_BUG             alphaPlusOne = alpha + 1._RK
!WSL_GFORTRAN_BUG             ! use brute-force integration
!WSL_GFORTRAN_BUG             call qag( f             = getBandCompLowEnergy  &
!WSL_GFORTRAN_BUG                     , a             = lowerLim              &
!WSL_GFORTRAN_BUG                     , b             = thisUpperLim          &
!WSL_GFORTRAN_BUG                     , epsabs        = 0._RK                 &
!WSL_GFORTRAN_BUG                     , epsrel        = tolerance             &
!WSL_GFORTRAN_BUG                     , key           = 1_IK                  &
!WSL_GFORTRAN_BUG                     , result        = energyFluence         &
!WSL_GFORTRAN_BUG                     , abserr        = abserr                &
!WSL_GFORTRAN_BUG                     , neval         = neval                 &
!WSL_GFORTRAN_BUG                     , ier           = ierr                  &
!WSL_GFORTRAN_BUG                     )
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG             if (ierr/=0_IK) then
!WSL_GFORTRAN_BUG                 energyFluence = -HUGE_RK
!WSL_GFORTRAN_BUG                 Err%occurred = .true.
!WSL_GFORTRAN_BUG                 Err%stat = ierr
!WSL_GFORTRAN_BUG                 Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred at QuadPack routine. Check the error code to identify the root cause."
!WSL_GFORTRAN_BUG                 return
!WSL_GFORTRAN_BUG             end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG             if (upperLim>ebrk) then ! add the remaining part of the energyFluence from the high-energy component
!WSL_GFORTRAN_BUG                 betaPlusTwo = beta + 2._RK
!WSL_GFORTRAN_BUG                 alphaMinusBeta = alpha - beta
!WSL_GFORTRAN_BUG                 coef = ebrk**alphaMinusBeta * exp(-alphaMinusBeta);
!WSL_GFORTRAN_BUG                 energyFluence = energyFluence + coef * ( upperLim**betaPlusTwo - ebrk**betaPlusTwo ) / betaPlusTwo
!WSL_GFORTRAN_BUG                 return
!WSL_GFORTRAN_BUG             end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         end if
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG     contains
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG         pure function getBandCompLowEnergy(energy) result(bandCompLow)
!WSL_GFORTRAN_BUG             implicit none
!WSL_GFORTRAN_BUG             real(RK), intent(in)    :: energy
!WSL_GFORTRAN_BUG             real(RK)                :: bandCompLow
!WSL_GFORTRAN_BUG             bandCompLow = energy**alphaPlusOne * exp(-alphaPlusTwoOverEpk*energy)
!WSL_GFORTRAN_BUG         end function getBandCompLowEnergy
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG !#if defined OS_IS_WSL && defined CODECOV_ENABLED
!WSL_GFORTRAN_BUG !! LCOV_EXCL_STOP
!WSL_GFORTRAN_BUG !#endif
!WSL_GFORTRAN_BUG 
!WSL_GFORTRAN_BUG     end subroutine getEnergyFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Convert an **input energy fluence** in `[lowerLim, upperLim]` energy window, all in units of keV (or ALL in some other units),
    !> **to photon fluence**, within the same energy range, or if `[lowerLimNew, upperLimNew]` is provided, then in that range.
    !> The input `tolerance` is passed to integrators as a measure of the desired accuracy for the computation of `photonFluence`.
    !>
    !> \param[in]   energyFluence   :   The input energy fluence (in units of [keV]) to be converted to the output photon fluence.
    !> \param[in]   lowerLim        :   The lower limit energy (in units of [keV]) of the integration.
    !> \param[in]   upperLim        :   The upper limit energy (in units of [keV]) of the integration.
    !> \param[in]   epk             :   The spectral peak energy in units of [keV].
    !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
    !> \param[in]   beta            :   The upper spectral exponent of the Band model.
    !> \param[in]   tolerance       :   The relative accuracy tolerance of the integration.
    !> \param[out]  photonFluence   :   The output fluence in units of photon counts.
    !> \param[out]  Err             :   An object of class [Err_type](@ref err_mod::err_type) containing error-handling information.
    !> \param[in]   lowerLimNew     :   The lower limit of energy windows (in keV) to be used for the computation of `photonFluence` (optional).
    !> \param[in]   upperLimNew     :   The upper limit of energy windows (in keV) to be used for the computation of `photonFluence` (optional).
    !>
    !> \remark
    !> If the optional `[lowerLimNew, upperLimNew]` are provided, each will replace the
    !> corresponding input `[lowerLim, upperLim]` in the computation of the output `photonFluence`.
    subroutine getPhotonFluenceFromEnergyFluence( energyFluence, lowerLim, upperLim, epk, alpha, beta, tolerance, photonFluence, Err, lowerLimNew, upperLimNew )
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

!#if defined OS_IS_WSL && defined CODECOV_ENABLED
!! LCOV_EXCL_START
!#endif

        lowLimNew = lowerLim
        if (present(lowerLimNew)) lowLimNew = lowerLimNew
        uppLimNew = upperLim
        if (present(upperLimNew)) uppLimNew = upperLimNew

        call getEnergyFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,amplitude,Err)
        if (Err%occurred) then
        ! LCOV_EXCL_START
            photonFluence = -HUGE_RK
            Err%msg = MODULE_NAME // PROCEDURE_NAME // Err%msg
            return
        end if
        ! LCOV_EXCL_STOP
        amplitude = energyFluence / amplitude

        call getPhotonFluence(lowLimNew,uppLimNew,epk,alpha,beta,tolerance,photonFluence,Err)
        if (Err%occurred) then
        ! LCOV_EXCL_START
            photonFluence = -HUGE_RK
            Err%msg = MODULE_NAME // PROCEDURE_NAME // Err%msg
            return
        end if
        ! LCOV_EXCL_STOP
        photonFluence = amplitude * photonFluence

!#if defined OS_IS_WSL && defined CODECOV_ENABLED
!! LCOV_EXCL_STOP
!#endif

    end subroutine getPhotonFluenceFromEnergyFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module BandSpectrum_mod ! LCOV_EXCL_LINE