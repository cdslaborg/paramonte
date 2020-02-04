!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module StarFormation_mod

    use Constants_mod, only: IK, RK, PI, NEGINF_RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@StarFormation_mod"

    real(RK)    , parameter :: ZMIN = 0.1e0_RK
    real(RK)    , parameter :: ZMAX = 1.0e2_RK

    abstract interface
    function getLogSFR_proc(zplus1,logzplus1,twiceLogLumDisMpc) result(logSFR)
        import :: RK, IK
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK)                :: logSFR
    end function getLogSFR_proc
    end interface

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogDensitySFRH06(logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogDensitySFRH06
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: logzplus1
        real(RK), parameter     :: logz0plus1 = log(1._RK+0.97_RK)
        real(RK), parameter     :: logz1plus1 = log(1._RK+4.50_RK)
        real(RK), parameter     :: g0 = +3.4_RK
        real(RK), parameter     :: g1 = -0.3_RK
        real(RK), parameter     :: g2 = -7.8_RK
        real(RK), parameter     :: logNormFac1 = logz0plus1*(g0-g1)
        real(RK), parameter     :: logNormFac2 = logz1plus1*(g1-g2) + logNormFac1
        real(RK)                :: logDensitySFR
        if (logzplus1<0._RK) then
            logDensitySFR = NEGINF_RK
        elseif (logzplus1<logz0plus1) then
            logDensitySFR = logzplus1*g0
        elseif (logzplus1<logz1plus1) then
            logDensitySFR = logzplus1*g1 + logNormFac1
        else
            logDensitySFR = logzplus1*g2 + logNormFac2
        end if
    end function getLogDensitySFRH06

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogDensitySFRL08(logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogDensitySFRL08
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: logzplus1
        real(RK), parameter     :: logz0plus1 = log(1._RK+0.993_RK)
        real(RK), parameter     :: logz1plus1 = log(1._RK+3.800_RK)
        real(RK), parameter     :: g0 = +3.3000_RK
        real(RK), parameter     :: g1 = +0.0549_RK
        real(RK), parameter     :: g2 = -4.4600_RK
        real(RK), parameter     :: logNormFac1 = logz0plus1*(g0-g1)
        real(RK), parameter     :: logNormFac2 = logz1plus1*(g1-g2) + logNormFac1
        real(RK)                :: logDensitySFR
        if (logzplus1<0._RK) then
            logDensitySFR = NEGINF_RK
        elseif (logzplus1<logz0plus1) then
            logDensitySFR = logzplus1*g0
        elseif (logzplus1<logz1plus1) then
            logDensitySFR = logzplus1*g1 + logNormFac1
        else
            logDensitySFR = logzplus1*g2 + logNormFac2
        end if
    end function getLogDensitySFRL08

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! Metalicity corrected SFR of Butler, Bloom et al 2010
    pure function getLogDensitySFRB10(logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogDensitySFRB10
#endif
        use Constants_mod, only: IK, RK
        implicit none
        real(RK), intent(in)    :: logzplus1
        real(RK), parameter     :: logz0plus1 = log(1._RK+0.97_RK)
        real(RK), parameter     :: logz1plus1 = log(1._RK+4.00_RK)
        real(RK), parameter     :: g0 = +3.14_RK
        real(RK), parameter     :: g1 = +1.36_RK
        real(RK), parameter     :: g2 = -2.92_RK
        real(RK), parameter     :: logNormFac1 = logz0plus1*(g0-g1)
        real(RK), parameter     :: logNormFac2 = logz1plus1*(g1-g2) + logNormFac1
        real(RK)                :: logDensitySFR
        if (logzplus1<0._RK) then
            logDensitySFR = NEGINF_RK
        elseif (logzplus1<logz0plus1) then
            logDensitySFR = logzplus1*g0
        elseif (logzplus1<logz1plus1) then
            logDensitySFR = logzplus1*g1 + logNormFac1
        else
            logDensitySFR = logzplus1*g2 + logNormFac2
        end if
    end function getLogDensitySFRB10

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! returns the Comoving Star Formation Rate Density according to Eqn 15 of Madau 2014: Cosmic Star-Formation History
    ! densitySFR(z) = 0.015 * (1+z)^2.7 / ( 1 + [(1+z)/2.9]^5.6 )
    pure function getLogDensitySFRM14(zplus1,logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogDensitySFRM14
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1
        real(RK)                :: logDensitySFR
        real(RK), parameter     :: logAmplitude = log(0.015_RK)
        real(RK), parameter     :: lowerExp = 2.7_RK
        real(RK), parameter     :: upperExp = 5.6_RK
        real(RK), parameter     :: zplus1Break = 2.9_RK
        real(RK), parameter     :: zplus1Coeff = 1._RK / (zplus1Break**upperExp)
        logDensitySFR = logAmplitude + lowerExp*logzplus1 - log( 1._RK + zplus1Coeff * zplus1**upperExp )
    end function getLogDensitySFRM14

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogSFRH06(zplus1,logzplus1,twiceLogLumDisMpc) result(logSFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSFRH06
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logSFR
        logSFR  = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogDensitySFRH06(logzplus1)
    end function getLogSFRH06

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogSFRL08(zplus1,logzplus1,twiceLogLumDisMpc) result(logSFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSFRL08
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logSFR
        logSFR  = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogDensitySFRL08(logzplus1)
    end function getLogSFRL08

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogSFRB10(zplus1,logzplus1,twiceLogLumDisMpc) result(logSFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSFRB10
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logSFR
        logSFR  = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogDensitySFRB10(logzplus1)
    end function getLogSFRB10

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    pure function getLogSFRM14(zplus1,logzplus1,twiceLogLumDisMpc) result(logSFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSFRM14
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logSFR
        logSFR  = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogDensitySFRM14(zplus1,logzplus1)
    end function getLogSFRM14

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module StarFormation_mod