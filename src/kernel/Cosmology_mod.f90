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

!> \brief This module contains procedures and constants for cosmological calculations.
!> \author Amir Shahmoradi

module Cosmology_mod

    use Constants_mod, only: RK, PI

    implicit none

    character(*), parameter :: MODULE_NAME = "@Cosmology_mod"

    ! Cosmological constants

    real(RK), parameter :: LIGHT_SPEED = 3.e5_RK                                    !< @public LIGHT_SPEED is the speed of light (Km/s).
    real(RK), parameter :: HUBBLE_CONST = 7.1e1_RK                                  !< @public HUBBLE_CONST is the Hubble constant in units of km/s/MPc.
    real(RK), parameter :: HUBBLE_TIME_GYRS = 13.8_RK                               !< @public hubble time (liddle 2003, page 57) in units of gyrs
    real(RK), parameter :: INVERSE_HUBBLE_CONST = 1._RK / HUBBLE_CONST              !< @public inverse of HUBBLE_CONST: 0.014084507042254
    real(RK), parameter :: LS2HC = LIGHT_SPEED / HUBBLE_CONST                       !< @public the speed of light in units of km/s divided by the Hubble constant.
    real(RK), parameter :: LOGLS2HC = log(LIGHT_SPEED) - log(HUBBLE_CONST)          !< @public log speed of light in units of km/s divided by the Hubble constant.
    real(RK), parameter :: MPC2CM = 3.09e24_RK                                      !< @public 1 Mega Parsec = MPC2CM centimeters.
    real(RK), parameter :: LOGMPC2CMSQ4PI   = log(4._RK*PI) + 2._RK*log(MPC2CM)     !< @public log(MegaParsec2centimeters).
    real(RK), parameter :: LOG10MPC2CMSQ4PI = log10(4._RK*PI) + 2._RK*log10(MPC2CM) !< @public log10(MegaParsec2centimeters).
    real(RK), parameter :: OMEGA_DE = 0.7_RK                                        !< @public Dark Energy density.
    real(RK), parameter :: OMEGA_DM = 0.3_RK                                        !< @public Dark Matter density.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the natural logarithm of the differential comoving volume of cosmos.
    !>
    !> @param[in]   zplus1              : The redshift plus 1.
    !> @param[in]   logzplus1           : The logarithm of redshift plus 1.
    !> @param[in]   twiceLogLumDisMpc   : The term representing the logarithm of the luminosity distance squared.
    !>
    !> \return
    !> `logdvdz` : The natural logarithm of the differential comoving volume of cosmos.
    pure function getlogdvdz(zplus1,logzplus1,twiceLogLumDisMpc) result(logdvdz)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getlogdvdz
#endif
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4*PI*LS2HC)
        real(RK)                :: logdvdz
        logdvdz = LOG_COEF + twiceLogLumDisMpc - ( 2*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) )
    end function getlogdvdz

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the approximate logarithm of the cosmological luminosity distance in units of MPc.
    !>
    !> @param[in]   zplus1 : The redshift plus 1.
    !>
    !> \return
    !> `logLumDisWicMpc` : The approximate logarithm of the cosmological luminosity distance.
    !>
    !> \warning
    !> The approximation is accurate with a relative error of 0.001 for any redshift above 0.1, or `zplus1 >= 1.1`.
    !! Note that for redshifts less than 0.1, the error in the calculated luminosity distance grows to more than 0.001.
    !! This algorithm should therefore not be used for zplus1<0.1.
    !>
    !> \remark
    !> The distance is calculated according to approximate algorithm of Wickramasinghe & Okwatta (2010).
    pure function getLogLumDisWicMpc(zplus1) result(logLumDisWicMpc)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogLumDisWicMpc
#endif
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

    !> \brief
    !> Return the approximation to the cosmological luminosity distance.
    !>
    !> @param[in]   zplus1 : The redshift plus 1.
    !>
    !> \return
    !> `ldiswickram` : The approximate cosmological luminosity distance.
    !>
    !> The approximation is accurate with a relative error of 0.001 for any redshift above 0.1.
    !>
    !> \remark
    !> This function is the same as [getLogLumDisWicMpc()](@ref getloglumdiswicmpc) function,
    !> except for the fact that it return the natural value, as opposed to the natural logarithm.
    !> It is kept only for legacy reasons and should not be used in new code.
    pure function ldiswickram(zplus1)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: ldiswickram
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1
        real(RK)                :: ldiswickram,alpha,alpha0,x,x0,psi,psi0
        alpha=1._RK+2*OMEGA_DE/(OMEGA_DM*zplus1**3)
        alpha0=1._RK+2*OMEGA_DE/OMEGA_DM
        x=log(alpha+sqrt(alpha*alpha-1._RK))
        x0=log(alpha0+sqrt(alpha0*alpha0-1._RK))
        psi=x**0.333333333333333_RK*(1.5874010519682_RK-6.2992105236833e-3_RK*x*x+7.5375168659459e-5_RK*x**4)
        psi0=x0**0.333333333333333_RK*(1.5874010519682_RK-6.2992105236833e-3_RK*x0*x0+7.5375168659459e-5_RK*x0**4)
        ldiswickram=LS2HC*zplus1*(psi0-psi)/(OMEGA_DE**0.1666666666666667_RK*OMEGA_DM**0.333333333333333_RK)
    end function ldiswickram

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the cosmological lookback time in GYrs at the given redshift for the assumed cosmological parameters.
    !>
    !> @param[in]   zplus1              : The redshift plus 1.
    !> @param[in]   maxRelativeError    : The maximum tolerance for error in the numerical integration (**optional**, default = 1.e-6).
    !> @param[in]   nRefinement         : The number of refinements in the Romberg numerical integration (**optional**, default = 5).
    !> \return
    !> `lookBackTime` : The cosmological lookback time in GYrs at the given redshift.
    !>
    !> \remark
    !> The rest of the required parameters are taken from the module:
    !> - `HUBBLE_TIME_GYRS` : The Hubble time in units of Gyrs (Giga Years). It should be a number close to 13.8 Gyrs (Liddle, 2003, Page 57).
    !>
    !> \remark
    !> The integrations are performed using the Romberg integration method.
    !>
    !> \author
    !> Amir Shahmoradi, Sunday 2:31 PM, January 6, 2013, IFS, The University of Texas at Austin.
    function getLookBackTime(zplus1,maxRelativeError,nRefinement) result(lookBackTime)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLookBackTime
#endif
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

        maxRelativeErrorDefault = 1.e-6_RK; if (present(maxRelativeError)) maxRelativeErrorDefault = maxRelativeError
        nRefinementDefault = 7_IK; if (present(nRefinement)) nRefinementDefault = nRefinement

        call doQuadRombClosed   ( getFunc           = getLookBackTimeDensity    &
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
            ! LCOV_EXCL_START
            write(output_unit,"(A)") ErrorMessage(ierr)
            error stop
            ! LCOV_EXCL_STOP
        end if
        lookBackTime = HUBBLE_TIME_GYRS * lookBackTime

    end function getLookBackTime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the **differential** (w.r.t. `z`) cosmological lookback time in GYrs at the given redshift for the assumed cosmological parameters.
    !>
    !> @param[in]   zplus1              : The redshift plus 1.
    !> \return
    !> `lookBackTimeDnesity` : The cosmological lookback time in GYrs at the given redshift.
    !>
    !> \remark
    !> The rest of the required parameters are taken from the module:
    !> - `OMEGA_DE` : the dark energy density,
    !> - `OMEGA_DM` : the dark matter density.
    !>
    !> \remark
    !> This function could have become an internal function within [getLookBackTime](@ref getlookbacktime).
    !> However, the library yields segmentation fault error when compiled and run on the Windows Subsystem 
    !> for Linux Ubuntu with GFortran. As such, it is implemented as an independent function.
    !>
    !> \author
    !> Amir Shahmoradi, Sunday 2:31 PM, January 6, 2013, IFS, The University of Texas at Austin.
    pure function getLookBackTimeDensity(zplus1) result(lookBackTimeDnesity)
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1
        real(RK)                :: lookBackTimeDnesity
        lookBackTimeDnesity = 1._RK / ( zplus1 * sqrt(OMEGA_DM * zplus1**3 + OMEGA_DE) )
    end function getLookBackTimeDensity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the derivative of the age of the Universe, w.r.t. redshift for a given input redshift + 1.
    !>
    !> @param[in] zplus1 : The redshift plus 1.
    !>
    !> \return
    !> `universeAgeDerivative` : The derivative of the age of the Universe, w.r.t. redshift for a given input `zplus1`.
    !>
    !> \remark
    !> The rest of the required parameters are taken from the module:
    !> - `INVERSE_HUBBLE_CONST`
    !> - `OMEGA_DE`
    !> - `OMEGA_DM`
    !>
    !> \remark
    !> The integrations are performed using the Romberg integration method.
    !>
    !> \author
    !> Amir Shahmoradi, Sunday 2:31 PM, January 6, 2013, IFS, The University of Texas at Austin.
    pure function getUniverseAgeDerivative(zplus1) result(universeAgeDerivative)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniverseAgeDerivative
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1
        real(RK)                :: universeAgeDerivative
        universeAgeDerivative = INVERSE_HUBBLE_CONST / ( zplus1 * sqrt(OMEGA_DM * zplus1**3 + OMEGA_DE) )
    end function getUniverseAgeDerivative

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Cosmology_mod