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

module Cosmology_mod

    use Constants_mod, only: RK, PI

    implicit none

    character(*), parameter :: MODULE_NAME = "@Cosmology_mod"

    ! Cosmological constants

    real(RK), parameter :: LIGHT_SPEED = 3.e5_RK                                    ! LIGHT_SPEED is the speed of light (Km/s).
    real(RK), parameter :: HUBBLE_CONST = 7.1e1_RK                                  ! HUBBLE_CONST is the Hubble constant in units of km/s/MPc.
    real(RK), parameter :: HUBBLE_TIME_GYRS = 13.8_RK                               ! hubble time (liddle 2003, page 57) in units of gyrs
    real(RK), parameter :: INVERSE_HUBBLE_CONST = 1._RK / HUBBLE_CONST              ! inverse of HUBBLE_CONST.
    real(RK), parameter :: LS2HC = LIGHT_SPEED / HUBBLE_CONST                       ! the speed of light in units of km/s divided by the Hubble constant.
    real(RK), parameter :: LOGLS2HC = log(LIGHT_SPEED) - log(HUBBLE_CONST)          ! log speed of light in units of km/s divided by the Hubble constant.
    real(RK), parameter :: MPC2CM = 3.09e24_RK                                      ! 1 Mega Parsec = MPC2CM centimeters.
    real(RK), parameter :: LOGMPC2CMSQ4PI   = log(4._RK*PI) + 2._RK*log(MPC2CM)     ! log(MegaParsec2centimeters).
    real(RK), parameter :: LOG10MPC2CMSQ4PI = log10(4._RK*PI) + 2._RK*log10(MPC2CM) ! log10(MegaParsec2centimeters).
    real(RK), parameter :: OMEGA_DE = 0.7_RK                                        ! Dark Energy density.
    real(RK), parameter :: OMEGA_DM = 0.3_RK                                        ! Dark Matter density.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getlogdvdz(zplus1,logzplus1,twiceLogLumDisMpc) result(logdvdz)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getlogdvdz
#endif
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logdvdz
        logdvdz = LOG_COEF + twiceLogLumDisMpc - ( 2._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+omega_DE) )
    end function getlogdvdz

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogLumDisWicMpc(zplus1) result(logLumDisWicMpc)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogLumDisWicMpc
#endif
        ! This function calculates the luminosity distance in units of Mpc. Input is redshift+1 (zplus1) and must be zplus1>=1.1.
        ! The distance is calculated according to approximate algorithm of Wickramasinghe & Okwatta (2010). 
        ! Note that for redshifts less than 0.1, the error in the calculated luminosity distance grows to more than 0.001.
        ! This algorithm should therefore not be used for zplus1<0.1.
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1   ! redshift + 1
        real(RK)                :: logLumDisWicMpc
        real(RK)                :: alpha1, x1, x1Sq, psix1
        real(RK), parameter     :: TWICE_OMEGA_DE_OVER_OMEGA_DM = 2._RK * OMEGA_DE / OMEGA_DM
        real(RK), parameter     :: ALPHA0 = 1._RK + TWICE_OMEGA_DE_OVER_OMEGA_DM 
        real(RK), parameter     :: X0 = log( ALPHA0 + sqrt(ALPHA0**2-1._RK) )
        real(RK), parameter     :: X0Sq = X0**2
        real(RK), parameter     :: PSI_COEF1 = 2._RK**(2._RK/3._RK)
        real(RK), parameter     :: PSI_COEF2 = -PSI_COEF1 / 252._RK
        real(RK), parameter     :: PSI_COEF3 = +PSI_COEF1 / 21060._RK
        real(RK), parameter     :: PSIX0 = X0**0.333333333333333_RK * ( PSI_COEF1 + X0Sq * (PSI_COEF2 + X0Sq*PSI_COEF3) )
        real(RK), parameter     :: LOG_COEF = log(LS2HC / (OMEGA_DE**0.1666666666666667_RK*OMEGA_DM**0.3333333333333333_RK))
        alpha1          = 1._RK + TWICE_OMEGA_DE_OVER_OMEGA_DM / zplus1**3
        x1              = log( alpha1 + sqrt( alpha1**2 - 1._RK ) )
        x1Sq            = x1**2
        psix1           = x1**0.333333333333333_RK * ( PSI_COEF1 + x1Sq * (PSI_COEF2 + x1Sq*PSI_COEF3) )
        logLumDisWicMpc = LOG_COEF +  log(zplus1*(PSIX0-psix1))
    end function getLogLumDisWicMpc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function ldiswickram(zplus1)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: ldiswickram
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1
        real(RK)                :: ldiswickram,alpha,alpha0,x,x0,psi,psi0
        alpha=1.d0+2.d0*OMEGA_DE/(OMEGA_DM*zplus1**3)
        alpha0=1.d0+2.d0*OMEGA_DE/OMEGA_DM
        x=log(alpha+sqrt(alpha*alpha-1.d0))
        x0=log(alpha0+sqrt(alpha0*alpha0-1.d0))
        psi=x**0.33333333*(1.5874010519682-6.2992105236833d-3*x*x+7.5375168659459d-5*x**4)
        psi0=x0**0.33333333*(1.5874010519682-6.2992105236833d-3*x0*x0+7.5375168659459d-5*x0**4)
        ldiswickram=LS2HC*zplus1*(psi0-psi)/(OMEGA_DE**0.16666666*OMEGA_DM**0.33333333)
    end function ldiswickram

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLookBackTime(zplus1,maxRelativeError,nRefinement) result(lookBackTime)
        !   This function calculates the lookback time in GYrs at given redshift for the assumed cosmological parameters.
        !   The input parameters are:
        !        -           zplus1: the redshift plus one
        !        - maxRelativeError: the maximum tolerance for error in the numerical integration  (default = 1.e-6)
        !        -      nRefinement: the number of refinements in Romberg's numerical integration  (default = 5)
        !   The rest of the required parameters are taken from the module:
        !       -  HUBBLE_TIME_GYRS: The Hubble time in units of Gyrs (Giga Years). It should be a number close to 13.8 Gyrs (Liddle, 2003, Page 57).
        !       -          OMEGA_DE: the dark energy density,
        !       -          OMEGA_DM: the dark matter density.
        !   The integration is performed using Romberg's interation method.
        !   Amir Shahmoradi, Sunday 2:31 PM, January 6, 2013, IFS, The University of Texas at Austin.

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Integration_mod, only: doQuadRombClosed, ErrorMessage
        use Constants_mod, only: RK, IK
        implicit none
        real(RK)    , intent(in)            :: zplus1
        real(RK)    , intent(in), optional  :: maxRelativeError
        integer(IK) , intent(in), optional  :: nRefinement
        real(RK)    , parameter             :: ZPLUS1_MIN = 1._RK
        real(RK)                            :: lookBackTime
        real(RK)                            :: maxRelativeErrorDefault, relerr
        integer(IK)                         :: neval, ierr, nRefinementDefault

        nRefinementDefault = 7_IK; if (present(nRefinement)) nRefinementDefault = nRefinement
        maxRelativeErrorDefault = 1.e-6_RK; if (present(maxRelativeError)) maxRelativeErrorDefault = maxRelativeError

        call doQuadRombClosed   ( getFunc           = getIntegrand              &
                                , lowerLim          = ZPLUS1_MIN                &
                                , upperLim          = zplus1                    &
                                , maxRelativeError  = maxRelativeErrorDefault   &
                                , nRefinement       = nRefinementDefault        &
                                , integral          = lookBackTime              &
                                , relativeError     = relerr                    &
                                , numFuncEval       = neval                     &
                                , ierr              = ierr                      &
                                )
        if (ierr/=0_IK) then
            write(output_unit,"(A)") ErrorMessage(ierr)
            error stop
        end if
        lookBackTime = HUBBLE_TIME_GYRS * lookBackTime

    contains

        function getIntegrand(zplus1) result(integrand)
            use Constants_mod, only: RK
            implicit none
            real(RK), intent(in)    :: zplus1
            real(RK)                :: integrand
            integrand = 1._RK / ( zplus1 * sqrt(OMEGA_DM * zplus1**3 + OMEGA_DE) )
        end function getIntegrand

    end function getLookBackTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getUniverseAgeDerivative(zplus1) result(universeAgeDerivative)
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1
        real(RK)                :: universeAgeDerivative
        universeAgeDerivative = INVERSE_HUBBLE_CONST / ( zplus1 * sqrt(OMEGA_DM * zplus1**3 + OMEGA_DE) )
    end function getUniverseAgeDerivative

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Cosmology_mod