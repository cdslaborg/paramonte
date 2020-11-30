submodule (BandSpectrum_mod) EnergyFluence_smod

    implicit none

    real(RK)                        :: mv_alphaPlusOne
    real(RK)                        :: mv_alphaPlusTwoOverEpk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Integrate the Band differential spectrum over the input energy range in units of the input energy.
    !>
    !> \param[in]   lowerLim        :   The lower limit energy (in units of [keV]) of the integration.
    !> \param[in]   upperLim        :   The upper limit energy (in units of [keV]) of the integration.
    !> \param[in]   epk             :   The spectral peak energy in units of [keV].
    !> \param[in]   alpha           :   The lower spectral exponent of the Band model.
    !> \param[in]   beta            :   The upper spectral exponent of the Band model.
    !> \param[in]   tolerance       :   The relative accuracy tolerance of the integration.
    !> \param[out]  energyFluence   :   The fluence in units of the input energy (keV) within the input energy range.
    !> \param[out]  Err             :   An object of class [Err_type](@ref err_mod::err_type) containing error-handling information.
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
        character(*), parameter         :: PROCEDURE_NAME = "@getEnergyFluence()"
        real(RK)                        :: ebrk, alphaPlusTwo
        real(RK)                        :: thisUpperLim, betaPlusTwo
        real(RK)                        :: alphaMinusBeta, coef
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

!#if defined OS_IS_WSL && defined CODECOV_ENABLED
!! LCOV_EXCL_START
!#endif

        elseif (lowerLim<ebrk) then

            mv_alphaPlusTwoOverEpk = alphaPlusTwo / epk
            thisUpperLim = min(upperLim,ebrk)
            mv_alphaPlusOne = alpha + 1._RK
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

    end subroutine getEnergyFluence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getBandCompLowEnergy(energy) result(bandCompLow)
        implicit none
        real(RK), intent(in)    :: energy
        real(RK)                :: bandCompLow
        bandCompLow = energy**mv_alphaPlusOne * exp(-mv_alphaPlusTwoOverEpk*energy)
    end function getBandCompLowEnergy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end submodule EnergyFluence_smod ! LCOV_EXCL_LINE