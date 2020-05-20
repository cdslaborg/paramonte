!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module Cosmology_mod

    use Constants_mod, only: RK, PI

    implicit none

    character(*), parameter :: MODULE_NAME = "@Cosmology_mod"

    ! Cosmological constants

    real(RK), parameter :: LIGHT_SPEED = 3.e5_RK                                    ! LIGHT_SPEED is the speed of light (Km/s).
    real(RK), parameter :: HUBBLE_CONST = 7.1e1_RK                                  ! HUBBLE_CONST is the Hubble constant in units of km/s/MPc.
    real(RK), parameter :: LS2HC = LIGHT_SPEED / HUBBLE_CONST                       ! the speed of light in units of km/s divided by the Hubble constant.
    real(RK), parameter :: LOGLS2HC = log(LIGHT_SPEED) - log(HUBBLE_CONST)          ! log speed of light in units of km/s divided by the Hubble constant.
    real(RK), parameter :: MPC2CM = 3.09e24_RK                                      ! 1 Mega Parsec = MPC2CM centimeters.
    real(RK), parameter :: LOGMPC2CMSQ4PI   = log(4._RK*PI) + 2._RK*log(MPC2CM)     ! log(MegaParsec2centimeters).
    real(RK), parameter :: LOG10MPC2CMSQ4PI = log10(4._RK*PI) + 2._RK*log10(MPC2CM) ! log10(MegaParsec2centimeters).
    real(RK), parameter :: OMEGA_DE = 0.7_RK                                        ! Dark Energy density.
    real(RK), parameter :: OMEGA_DM = 0.3_RK                                        ! Dark Matter density.

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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
        real(RK), parameter     :: PSIX0 = X0**0.33333333 * ( PSI_COEF1 + X0Sq * (PSI_COEF2 + X0Sq*PSI_COEF3) )
        real(RK), parameter     :: LOG_COEF = log(LS2HC / (OMEGA_DE**0.1666666666666667_RK*OMEGA_DM**0.3333333333333333_RK))
        alpha1          = 1._RK + TWICE_OMEGA_DE_OVER_OMEGA_DM / zplus1**3
        x1              = log( alpha1 + sqrt( alpha1**2 - 1._RK ) )
        x1Sq            = x1**2
        psix1           = x1**0.33333333 * ( PSI_COEF1 + x1Sq * (PSI_COEF2 + x1Sq*PSI_COEF3) )
        logLumDisWicMpc = LOG_COEF +  log(zplus1*(PSIX0-psix1))
    end function getLogLumDisWicMpc

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Cosmology_mod