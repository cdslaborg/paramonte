submodule (BandSpectrum_mod) PhotonFluence_smod

    implicit none

    real(RK)                        :: mv_alpha
    real(RK)                        :: mv_alphaPlusTwoOverEpk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Integrate the Band differential spectrum over the input energy range.
    !>
    !> \param[in]   lowerLim        :   The lower limit energy (in units of [keV]) of the integration.
    !> \param[in]   upperLim        :   The upper limit energy (in units of [keV]) of the integration.
    !> \param[in]   epk             :   The spectral peak energy in units of [keV].
    !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
    !> \param[in]   beta            :   The upper spectral exponent of the Band model.
    !> \param[in]   tolerance       :   The relative accuracy tolerance of the integration.
    !> \param[out]  photonFluence   :   The fluence in units of photon counts within the input energy range.
    !> \param[out]  Err             :   An object of class [Err_type](@ref err_mod::err_type) containing error-handling information.
    module subroutine getPhotonFluence(lowerLim,upperLim,epk,alpha,beta,tolerance,photonFluence,Err)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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
        real(RK)                        :: thisUpperLim, betaPlusOne
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

            mv_alpha = alpha

            mv_alphaPlusTwoOverEpk = alphaPlusTwo / epk
            thisUpperLim = min(upperLim,ebrk)
            !alphaPlusOne = alpha + 1._RK
            !if (alpha>-1._RK) then
            !    logGammaAlphaPlusOne = log_gamma( alphaPlusOne )
            !    ! use the analytical approach to compute the photonFluence:
            !    ! https://www.wolframalpha.com/input/?i=integrate+x%5Ea+*+exp(-b*x)
            !    photonFluence = getUpperGamma( exponent = alphaPlusOne &
            !                            , logGammaExponent = logGammaAlphaPlusOne &
            !                            , lowerLim = mv_alphaPlusTwoOverEpk * lowerLim &
            !                            , tolerance = tolerance &
            !                            ) &
            !             - getUpperGamma( exponent = alphaPlusOne &
            !                            , logGammaExponent = logGammaAlphaPlusOne &
            !                            , lowerLim = mv_alphaPlusTwoOverEpk * thisUpperLim &
            !                            , tolerance = tolerance &
            !                            )
            !    photonFluence = photonFluence / mv_alphaPlusTwoOverEpk**alphaPlusOne
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
                ! LCOV_EXCL_START
                    photonFluence = -HUGE_RK
                    Err%occurred = .true.
                    Err%stat = ierr
                    Err%msg = MODULE_NAME // PROCEDURE_NAME // ": Error occurred at QuadPack routine. Check the error code to identify the root cause."
                    return
                ! LCOV_EXCL_STOP
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

    end subroutine getPhotonFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getBandCompLowPhoton(energy) result(bandCompLow)
        implicit none
        real(RK), intent(in)    :: energy
        real(RK)                :: bandCompLow
        bandCompLow = energy**mv_alpha * exp(-mv_alphaPlusTwoOverEpk*energy)
    end function getBandCompLowPhoton

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule PhotonFluence_smod ! LCOV_EXCL_LINE