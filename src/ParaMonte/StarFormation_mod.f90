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

module StarFormation_mod

    use Constants_mod, only: IK, RK, PI, NEGINF_RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@StarFormation_mod"

    real(RK)    , parameter :: ZMIN = 0.1e0_RK
    real(RK)    , parameter :: ZMAX = 1.0e2_RK

    abstract interface
    function getLogRate_proc(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
        import :: RK, IK
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK)                :: logRate
    end function getLogRate_proc
    end interface

    abstract interface
    function getRateDensity_proc(zplus1) result(rateDensity)
        import :: RK, IK
        real(RK), intent(in)    :: zplus1
        real(RK)                :: rateDensity
    end function getRateDensity_proc
    end interface

    abstract interface
    function getMergerDelayTimePDF_proc(mergerDelayTime) result(mergerDelayTimeProb)
        import :: RK, IK
        real(RK), intent(in)    :: mergerDelayTime
        real(RK)                :: mergerDelayTimeProb
    end function getMergerDelayTimePDF_proc
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! GRBFR based on Petrosian et al (2015)
    pure function getLogRateDensityP15(logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityP15
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: logzplus1
        real(RK), parameter     :: logz1plus1 = log(1._RK+4.50_RK)
        real(RK), parameter     :: exponentHighZ = -7.8_RK
        real(RK), parameter     :: logNormFac2 = -exponentHighZ * logz1plus1
        real(RK)                :: logDensitySFR
        if (logzplus1<0._RK) then
            logDensitySFR = NEGINF_RK
        elseif (logzplus1<logz1plus1) then
            logDensitySFR = 0._RK
        else
            logDensitySFR = logzplus1*exponentHighZ + logNormFac2
        end if
    end function getLogRateDensityP15

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateDensityH06(logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityH06
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
    end function getLogRateDensityH06

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateDensityL08(logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityL08
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
    end function getLogRateDensityL08

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Metalicity corrected SFR of Butler, Bloom et al 2010
    pure function getLogRateDensityB10(logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityB10
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
    end function getLogRateDensityB10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the Comoving Star Formation Rate Density according to Eqn 15 of Madau 2014: Cosmic Star-Formation History
    ! densitySFR(z) = 0.015 * (1+z)^2.7 / ( 1 + [(1+z)/2.9]^5.6 )
    pure function getLogRateDensityM14(zplus1,logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM14
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
    end function getLogRateDensityM14

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the Comoving Star Formation Rate Density according to Eqn 1 of Madau 2017: Cosmic Star-Formation History
    ! densitySFR(z) = 0.01 * (1+z)^2.6 / ( 1 + [(1+z)/3.2]^6.2 )
    pure function getLogRateDensityM17(zplus1,logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityM17
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1
        real(RK)                :: logDensitySFR
        real(RK), parameter     :: logAmplitude = log(0.01_RK)
        real(RK), parameter     :: lowerExp = 2.6_RK
        real(RK), parameter     :: upperExp = 6.2_RK
        real(RK), parameter     :: zplus1Break = 3.2_RK
        real(RK), parameter     :: zplus1Coeff = 1._RK / (zplus1Break**upperExp)
        logDensitySFR = logAmplitude + lowerExp*logzplus1 - log( 1._RK + zplus1Coeff * zplus1**upperExp )
    end function getLogRateDensityM17

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Mordau Comoving Star Formation Rate Density with updated parameters from Fermi 2018
    ! densitySFR(z) = 0.013 * (1+z)^2.99 / ( 1 + [(1+z)/2.63]^6.19 )
    pure function getLogRateDensityF18(zplus1,logzplus1) result(logDensitySFR)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateDensityF18
#endif
        use Constants_mod, only: RK
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1
        real(RK)                :: logDensitySFR
        real(RK), parameter     :: logAmplitude = log(0.013_RK)
        real(RK), parameter     :: lowerExp = 2.99_RK
        real(RK), parameter     :: upperExp = 6.19_RK
        real(RK), parameter     :: zplus1Break = 2.63_RK
        real(RK), parameter     :: zplus1Coeff = 1._RK / (zplus1Break**upperExp)
        logDensitySFR = logAmplitude + lowerExp*logzplus1 - log( 1._RK + zplus1Coeff * zplus1**upperExp )
    end function getLogRateDensityF18

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateP15(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateP15
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logRate
        logRate = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogRateDensityP15(logzplus1)
    end function getLogRateP15

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateH06(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateH06
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logRate
        logRate = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogRateDensityH06(logzplus1)
    end function getLogRateH06

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateL08(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateL08
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logRate
        logRate = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogRateDensityL08(logzplus1)
    end function getLogRateL08

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateB10(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateB10
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logRate
        logRate = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogRateDensityB10(logzplus1)
    end function getLogRateB10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateM14(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateM14
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logRate
        logRate = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogRateDensityM14(zplus1,logzplus1)
    end function getLogRateM14

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateM17(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateM17
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logRate
        logRate = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogRateDensityM17(zplus1,logzplus1)
    end function getLogRateM17

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogRateF18(zplus1,logzplus1,twiceLogLumDisMpc) result(logRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogRateF18
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK), intent(in)    :: zplus1, logzplus1, twiceLogLumDisMpc
        real(RK), parameter     :: LOG_COEF = log(4._RK*PI*LS2HC)
        real(RK)                :: logRate
        logRate = LOG_COEF + twiceLogLumDisMpc - ( 3._RK*logzplus1 + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) &
                + getLogRateDensityF18(zplus1,logzplus1)
    end function getLogRateF18

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getBinaryMergerRate( zplus1 &
                                , zplus1Max &
                                , nRefinement &
                                , maxRelativeError &
                                , getMergerDelayTimePDF &
                                , getStarFormationRateDensity &
                                ) result(binaryMergerRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinaryMergerRate
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE, getLogLumDisWicMpc
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)                :: zplus1
        real(RK)    , intent(in), optional      :: zplus1Max, maxRelativeError
        integer(IK) , intent(in), optional      :: nRefinement
        procedure(getRateDensity_proc)          :: getStarFormationRateDensity
        procedure(getMergerDelayTimePDF_proc)   :: getMergerDelayTimePDF
        real(RK)                                :: binaryMergerRate
        real(RK), parameter                     :: LOG_COEF = log(4._RK*PI*LS2HC)

        binaryMergerRate    = exp( LOG_COEF + 2_IK*getLogLumDisWicMpc(zplus1) - ( 3._RK*log(zplus1) + 0.5_RK*log(OMEGA_DM*zplus1**3+OMEGA_DE) ) ) &
                            * getBinaryMergerRateDensity( zplus1 &
                                                        , zplus1Max &
                                                        , nRefinement &
                                                        , maxRelativeError &
                                                        , getMergerDelayTimePDF &
                                                        , getStarFormationRateDensity &
                                                        )

    end function getBinaryMergerRate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getBinaryMergerRateDensity ( zplus1 &
                                        , zplus1Max &
                                        , nRefinement &
                                        , maxRelativeError &
                                        , getMergerDelayTimePDF &
                                        , getStarFormationRateDensity &
                                        ) result(binaryMergerRateDensity)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinaryMergerRateDensity
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: RK, HUGE_RK
        use Cosmology_mod, only: getLookBackTime
        use Integration_mod, only: doQuadRombOpen, ErrorMessage!, midinf
        use Integration_mod, only: midexp
       !use Integration_mod, only: midinf
        implicit none
        real(RK)    , intent(in)                :: zplus1
        real(RK)    , intent(in), optional      :: zplus1Max, maxRelativeError
        integer(IK) , intent(in), optional      :: nRefinement
        procedure(getRateDensity_proc)          :: getStarFormationRateDensity
        procedure(getMergerDelayTimePDF_proc)   :: getMergerDelayTimePDF
        integer(IK)                             :: neval, ierr, nRefinementDefault
        real(RK)                                :: zplus1MaxDefault, maxRelativeErrorDefault, relerr
        real(RK)                                :: binaryMergerRateDensity, lookBackTimeRef

        nRefinementDefault = 7_IK; if (present(nRefinement)) nRefinementDefault = nRefinement
        zplus1MaxDefault = HUGE_RK; if (present(zplus1Max)) zplus1MaxDefault = zplus1Max
        maxRelativeErrorDefault = 1.e-6_RK; if (present(maxRelativeError)) maxRelativeErrorDefault = maxRelativeError

        lookBackTimeRef = getLookBackTime   ( zplus1 = zplus1 &
                                            , maxRelativeError = maxRelativeErrorDefault &
                                            , nRefinement = nRefinementDefault &
                                            )

        call doQuadRombOpen ( getFunc           = getBinaryMergerRateDensityIntegrand   &
                            , integrate         = midexp                                &
                           !, integrate         = midinf                                &
                            , lowerLim          = zplus1                                &
                            , upperLim          = zplus1MaxDefault                      &
                            , maxRelativeError  = maxRelativeErrorDefault               &
                            , nRefinement       = nRefinementDefault                    &
                            , integral          = binaryMergerRateDensity               &
                            , relativeError     = relerr                                &
                            , numFuncEval       = neval                                 &
                            , ierr              = ierr                                  &
                            )
        if (ierr/=0_IK) then
            write(output_unit,"(A)") ErrorMessage(ierr)
            error stop
        end if

    contains

        function getBinaryMergerRateDensityIntegrand(zplus1) result(binaryMergerRateIntegrand)

            use Cosmology_mod, only: getUniverseAgeDerivative
            implicit none
            real(RK)    , intent(in)    :: zplus1
            real(RK)                    :: binaryMergerRateIntegrand !,lognormpdf
            real(RK)                    :: mergerDelayTime

            ! note that zp<z always, so that delay>0.
            mergerDelayTime = getLookBackTime   ( zplus1 = zplus1 &
                                                , maxRelativeError = maxRelativeErrorDefault &
                                                , nRefinement = nRefinementDefault &
                                                ) 
            mergerDelayTime = mergerDelayTime - lookBackTimeRef    
            if (mergerDelayTime<=0._RK) then
                write(output_unit,"(A)") "The mergerDelayTime is non-positive in getBinaryMergerRateDensityIntegrand(): (zplus1, mergerDelayTime) = ", zplus1, mergerDelayTime
                error stop
            end if

            binaryMergerRateIntegrand   = getMergerDelayTimePDF(mergerDelayTime) &
                                        * getStarFormationRateDensity(zplus1) &
                                        * getUniverseAgeDerivative(zplus1)

        end function getBinaryMergerRateDensityIntegrand

    end function getBinaryMergerRateDensity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Amir Shahmoradi, Wednesday 5:43 PM, December 25, 2013, IFS, UT Austin
    ! compute the binary merger rate as a function of z, according to Shahmoradi and Nemiroff (2015)
    ! equivalent to delayed_rate_Belz_Li(z) in S15
    ! returns 0, if z>6.501_RK or z<0.09_RK
    pure function getBinaryMergerRateS15(z) result(binaryMergerRateS15)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinaryMergerRateS15
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)    :: z
        real(RK)                    :: binaryMergerRateS15
        if (z>2.5_RK .and. z<=6.501_RK) then
            binaryMergerRateS15 = &
                                - 2.09118024744342000_RK &
                                + 5.15382361299299000_RK * z &
                                - 5.46442271664195000_RK * z**2 &
                                + 3.29445310883082000_RK * z**3 &
                                - 1.24547016168265000_RK * z**4 &
                                + 0.30628893690508400_RK * z**5 &
                                - 0.04904403249641820_RK * z**6 &
                                + 0.00493757380504717_RK * z**7 &
                                - 2.84061971928750e-4_RK * z**8 &
                                + 7.12674138757750e-6_RK * z**9
        elseif (z>1.0_RK .and. z<=2.5_RK) then
            binaryMergerRateS15 = &
                                - 0.86022576265904100_RK &
                                + 4.22669545558817000_RK * z &
                                - 8.86086728534670000_RK * z**2 &
                                + 10.4863792284648000_RK * z**3 &
                                - 7.64722909221129000_RK * z**4 &
                                + 3.51616699500767000_RK * z**5 &
                                - 0.99555474471022000_RK * z**6 &
                                + 0.15876893754371900_RK * z**7 &
                                - 0.01092541997736420_RK * z**8
        elseif (z<=1._RK .and. z>=0.09_RK) then
            binaryMergerRateS15 = &
                                + 1.92595299989370e-4_RK &
                                - 0.00345273599582578_RK * z &
                                + 0.03157500615320920_RK * z**2 &
                                - 0.04470545521198460_RK * z**3 &
                                + 0.06812481521281660_RK * z**4 &
                                - 0.03846033416253570_RK * z**5
        else
            binaryMergerRateS15 = 0._RK
        end if
    end function getBinaryMergerRateS15

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the natural log of the binary merger rate as a function of logzplus1.
    ! However, note that the computed rate is dN/dz, even though the input is logzplus1.
    ! the merger delay time distribution is the same as that of Shahmoradi and Nemiroff (2015), but with the logMean:log(0.1_RK) and sigma: 0.9612813_RK.
    ! returns 0, if z>19.929999999999882_RK or z<0.03_RK
    pure function getLogBinaryMergerRateLognormB10(logzplus1) result(logBinaryMergerRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBinaryMergerRateLognormB10
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)    :: logzplus1
        real(RK)                    :: logBinaryMergerRate
        if (logzplus1>0.02955880224154443_RK .and. logzplus1<=0.20701416938432557_RK) then
            logBinaryMergerRate = & 
                                - 15.27802857671202000_RK &
                                + 94.54179164991284000_RK * logzplus1 &
                                - 687.3676159275769000_RK * logzplus1**2 &
                                + 2695.420977251770600_RK * logzplus1**3 &
                                - 4077.601406650646000_RK * logzplus1**4
        elseif (logzplus1>0.20701416938432557_RK .and. logzplus1<=0.8241754429663476_RK) then
            logBinaryMergerRate = &
                                - 13.5066182170954650_RK &
                                + 40.1985219822299200_RK * logzplus1 & 
                                - 121.506350703598660_RK * logzplus1**2 & 
                                + 224.621285123736100_RK * logzplus1**3 &
                                - 210.878836655472500_RK * logzplus1**4 &
                                + 76.3335749498628400_RK * logzplus1**5
        elseif (logzplus1>0.8241754429663476_RK .and. logzplus1<=1.4243124283074096_RK) then
            logBinaryMergerRate = &
                                - 10.0515447816113700_RK &
                                + 12.6659826494097970_RK * logzplus1 & 
                                - 13.2268991886238200_RK * logzplus1**2 & 
                                + 6.84523627043807100_RK * logzplus1**3 &
                                - 1.44645280124922220_RK * logzplus1**4
        elseif (logzplus1>1.4243124283074096_RK .and. logzplus1<=1.6104374127671848_RK) then
            logBinaryMergerRate = &
                                - 1187.90539057029950_RK &
                                + 3240.19327021926350_RK * logzplus1 & 
                                - 3330.70645904271000_RK * logzplus1**2 & 
                                + 1522.87499612399850_RK * logzplus1**3 &
                                - 261.341408956542300_RK * logzplus1**4
        elseif (logzplus1>1.6104374127671848_RK .and. logzplus1<=3.0411835364579027_RK) then
            logBinaryMergerRate = &
                                - 1.43934839576471260_RK &
                                + 1.72951867017028120_RK * logzplus1 & 
                                - 4.06729555225025000_RK * logzplus1**2 & 
                                + 1.18253386764330200_RK * logzplus1**3 &
                                - 0.15201156018584210_RK * logzplus1**4
        elseif (logzplus1<=0.02955880224154443_RK .or. logzplus1>3.0411835364579027_RK) then
            logBinaryMergerRate = 0._RK
        end if
    end function getLogBinaryMergerRateLognormB10

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the natural log of the binary merger rate as a function of logzplus1.
    ! However, note that the computed rate is dN/dz, even though the input is logzplus1.
    ! the merger delay time distribution is the same as that of Shahmoradi and Nemiroff (2015), but with the logMean:log(0.1_RK) and sigma: 0.9612813_RK..
    ! returns 0, if z>19.929999999999882_RK or z<0.03_RK
    pure function getLogBinaryMergerRateLognormF18(logzplus1) result(logBinaryMergerRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBinaryMergerRateLognormF18
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)    :: logzplus1
        real(RK)                    :: logBinaryMergerRate
        if (logzplus1>0.02955880224154443_RK .and. logzplus1<=0.16551443847757297_RK) then
            logBinaryMergerRate = & 
                                - 13.80128475140318000_RK &
                                + 79.17963739241087000_RK * logzplus1 &
                                - 420.9808813943490700_RK * logzplus1**2 &
                                + 902.4149755380632000_RK * logzplus1**3 
        elseif (logzplus1>0.16551443847757297_RK .and. logzplus1<=0.9282193027394269_RK) then
            logBinaryMergerRate = &
                                - 10.89131941401880800_RK &
                                + 21.59483273763076400_RK * logzplus1 & 
                                - 33.07662054750123000_RK * logzplus1**2 & 
                                + 29.23608723915102600_RK * logzplus1**3 &
                                - 11.33984422193848700_RK * logzplus1**4 
        elseif (logzplus1>0.9282193027394269_RK .and. logzplus1<=1.3937663759585892_RK) then
            logBinaryMergerRate = &
                                - 14.02518830695731000_RK &
                                + 24.92009885816913300_RK * logzplus1 & 
                                - 20.04762612951797000_RK * logzplus1**2 & 
                                + 4.885281389900379500_RK * logzplus1**3 &
                                - 0.168402813838908260_RK * logzplus1**4
        elseif (logzplus1>1.3937663759585892_RK .and. logzplus1<=3.0411835364579027_RK) then
            logBinaryMergerRate = &
                                - 4.348081430972656000_RK &
                                + 4.815143234949144000_RK * logzplus1 & 
                                - 6.143880845780776000_RK * logzplus1**2 & 
                                + 1.738835623950871300_RK * logzplus1**3 &
                                - 0.206972882929076480_RK * logzplus1**4
        elseif (logzplus1<=0.02955880224154443_RK .or. logzplus1>3.0411835364579027_RK) then
            logBinaryMergerRate = 0._RK
        end if
    end function getLogBinaryMergerRateLognormF18
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the natural log of the binary merger rate as a function of logzplus1.
    ! However, note that the computed rate is dN/dz, even though the input is logzplus1.
    ! the merger delay time distribution is the same as that of Shahmoradi and Nemiroff (2015), but with the logMean:log(0.1_RK) and sigma: 0.9612813_RK..
    ! returns 0, if z>19.929999999999882_RK or z<0.03_RK
    pure function getLogBinaryMergerRateLognormH06(logzplus1) result(logBinaryMergerRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBinaryMergerRateLognormH06
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)    :: logzplus1
        real(RK)                    :: logBinaryMergerRate
        if (logzplus1>0.02955880224154443_RK .and. logzplus1<=0.1441003439737565_RK) then
            logBinaryMergerRate = & 
                                - 14.26464149493092000_RK &
                                + 84.73477757043948000_RK * logzplus1 &
                                - 488.5893985602366500_RK * logzplus1**2 &
                                + 1154.414655194473900_RK * logzplus1**3 
        elseif (logzplus1>0.1441003439737565_RK .and. logzplus1<=0.6575200029167926_RK) then
            logBinaryMergerRate = &
                                - 11.19700066921606300_RK &
                                + 20.46712963401572300_RK * logzplus1 & 
                                - 24.31794334813894300_RK * logzplus1**2 & 
                                + 12.21213317590724400_RK * logzplus1**3 
        elseif (logzplus1>0.6575200029167926_RK .and. logzplus1<=1.5591966959973538_RK) then
            logBinaryMergerRate = &
                                - 9.094912666461765000_RK &
                                + 15.23119806754538900_RK * logzplus1 & 
                                - 18.77526325204311800_RK * logzplus1**2 & 
                                + 9.941360355936961000_RK * logzplus1**3 &
                                - 2.077370913197473000_RK * logzplus1**4
        elseif (logzplus1>1.5591966959973538_RK .and. logzplus1<=1.7056567701746455_RK) then
            logBinaryMergerRate = &
                                - 2392.907733171019000_RK &
                                + 6210.872744126407000_RK * logzplus1 & 
                                - 6054.866136454215000_RK * logzplus1**2 & 
                                + 2622.628785434413700_RK * logzplus1**3 &
                                - 426.0273477222719000_RK * logzplus1**4
        elseif (logzplus1>1.7056567701746455_RK .and. logzplus1<=3.0411835364579027_RK) then
            logBinaryMergerRate = &
                                + 9.538876239886940000_RK &
                                - 8.753418172517534000_RK * logzplus1 & 
                                - 0.159980818030374640_RK * logzplus1**2 & 
                                - 0.088551503657680930_RK * logzplus1**3 
        elseif (logzplus1<=0.02955880224154443_RK .or. logzplus1>3.0411835364579027_RK) then
            logBinaryMergerRate = 0._RK
        end if
    end function getLogBinaryMergerRateLognormH06

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the natural log of the binary merger rate as a function of logzplus1.
    ! However, note that the computed rate is dN/dz, even though the input is logzplus1.
    ! the merger delay time distribution is the same as that of Shahmoradi and Nemiroff (2015), but with the logMean:log(0.1_RK) and sigma: 0.9612813_RK..
    ! returns 0, if z>19.929999999999882_RK or z<0.03_RK
    pure function getLogBinaryMergerRateLognormL08(logzplus1) result(logBinaryMergerRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBinaryMergerRateLognormL08
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)    :: logzplus1
        real(RK)                    :: logBinaryMergerRate
        if (logzplus1>0.02955880224154443_RK .and. logzplus1<=0.20701416938432557_RK) then
            logBinaryMergerRate = & 
                                - 14.53696144309043900_RK &
                                + 94.70274747509626000_RK * logzplus1 &
                                - 687.3663996060040000_RK * logzplus1**2 &
                                + 2695.421036673770700_RK * logzplus1**3 &
                                - 4077.601561165490000_RK * logzplus1**4
        elseif (logzplus1>0.20701416938432557_RK .and. logzplus1<=0.8241754429663476_RK) then
            logBinaryMergerRate = &
                                - 13.5104005566057670_RK &
                                + 49.6443928683743600_RK * logzplus1 & 
                                - 164.286063097338630_RK * logzplus1**2 & 
                                + 315.721394966368100_RK * logzplus1**3 &
                                - 300.345052726248640_RK * logzplus1**4 &
                                + 108.470535327547080_RK * logzplus1**5
        elseif (logzplus1>0.8241754429663476_RK .and. logzplus1<=1.4243124283074096_RK) then
            logBinaryMergerRate = &
                                - 8.77634469738400500_RK &
                                + 13.1999684738558810_RK * logzplus1 & 
                                - 15.8698236818922140_RK * logzplus1**2 & 
                                + 8.48676936452957000_RK * logzplus1**3 &
                                - 1.83190451512279620_RK * logzplus1**4
        elseif (logzplus1>1.4243124283074096_RK .and. logzplus1<=1.6154199841116488_RK) then
            logBinaryMergerRate = &
                                + 4158.29353781047900_RK &
                                - 10954.1105856433040_RK * logzplus1 & 
                                + 10789.3451136201870_RK * logzplus1**2 & 
                                - 4713.80244702217800_RK * logzplus1**3 &
                                + 770.488645040204600_RK * logzplus1**4
        elseif (logzplus1>1.6154199841116488_RK .and. logzplus1<=3.0411835364579027_RK) then
            logBinaryMergerRate = &
                                + 0.37742655174185624_RK &
                                + 0.30883738015163340_RK * logzplus1 & 
                                - 4.04937550957291800_RK * logzplus1**2 & 
                                + 1.11680537027038170_RK * logzplus1**3 &
                                - 0.13770838345089523_RK * logzplus1**4
        elseif (logzplus1<=0.02955880224154443_RK .or. logzplus1>3.0411835364579027_RK) then
            logBinaryMergerRate = 0._RK
        end if
    end function getLogBinaryMergerRateLognormL08
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the natural log of the binary merger rate as a function of logzplus1.
    ! However, note that the computed rate is dN/dz, even though the input is logzplus1.
    ! the merger delay time distribution is the same as that of Shahmoradi and Nemiroff (2015), but with the logMean:log(0.1_RK) and sigma: 0.9612813_RK..
    ! returns 0, if z>19.929999999999882_RK or z<0.03_RK
    pure function getLogBinaryMergerRateLognormM14(logzplus1) result(logBinaryMergerRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBinaryMergerRateLognormM14
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)    :: logzplus1
        real(RK)                    :: logBinaryMergerRate
        if (logzplus1>0.02955880224154443_RK .and. logzplus1<=0.16551443847757297_RK) then
            logBinaryMergerRate = & 
                                - 13.91129314580349600_RK &
                                + 78.88963489621422000_RK * logzplus1 &
                                - 420.9801740859396700_RK * logzplus1**2 &
                                + 902.4783078800951000_RK * logzplus1**3 
        elseif (logzplus1>0.16551443847757297_RK .and. logzplus1<=0.9282193027394269_RK) then
            logBinaryMergerRate = &
                                - 11.00951046136480500_RK &
                                + 21.38817515748999000_RK * logzplus1 & 
                                - 33.29451048508970000_RK * logzplus1**2 & 
                                + 29.32158835244860400_RK * logzplus1**3 &
                                - 10.97377449040440000_RK * logzplus1**4 
        elseif (logzplus1>0.9282193027394269_RK .and. logzplus1<=1.3937663759585892_RK) then
            logBinaryMergerRate = &
                                - 8.254476015464371000_RK &
                                + 3.620963332444886000_RK * logzplus1 & 
                                + 6.734585433384001000_RK * logzplus1**2 & 
                                - 9.151412394211048000_RK * logzplus1**3 &
                                + 2.516171777428496000_RK * logzplus1**4
        elseif (logzplus1>1.3937663759585892_RK .and. logzplus1<=3.0411835364579027_RK) then
            logBinaryMergerRate = &
                                - 6.539697727782377000_RK &
                                + 8.522331572601950000_RK * logzplus1 & 
                                - 8.242990979412244000_RK * logzplus1**2 & 
                                + 2.316632169715435300_RK * logzplus1**3 &
                                - 0.266462340853027450_RK * logzplus1**4
        elseif (logzplus1<=0.02955880224154443_RK .or. logzplus1>3.0411835364579027_RK) then
            logBinaryMergerRate = 0._RK
        end if
    end function getLogBinaryMergerRateLognormM14

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! computes the natural log of the binary merger rate as a function of logzplus1.
    ! However, note that the computed rate is dN/dz, even though the input is logzplus1.
    ! the merger delay time distribution is the same as that of Shahmoradi and Nemiroff (2015), but with the logMean:log(0.1_RK) and sigma: 0.9612813_RK..
    ! returns 0, if z>19.929999999999882_RK or z<0.03_RK
    pure function getLogBinaryMergerRateLognormM17(logzplus1) result(logBinaryMergerRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogBinaryMergerRateLognormM17
#endif
        use Cosmology_mod, only: LS2HC, OMEGA_DM, OMEGA_DE
        use Constants_mod, only: RK, PI
        implicit none
        real(RK)    , intent(in)    :: logzplus1
        real(RK)                    :: logBinaryMergerRate
        if (logzplus1>0.02955880224154443_RK .and. logzplus1<=0.16551443847757297_RK) then
            logBinaryMergerRate = & 
                                - 14.01939141013502300_RK &
                                + 78.80010843737509000_RK * logzplus1 &
                                - 420.9593253775164000_RK * logzplus1**2 &
                                + 902.5668042795056000_RK * logzplus1**3 
        elseif (logzplus1>0.16551443847757297_RK .and. logzplus1<=0.9282193027394269_RK) then
            logBinaryMergerRate = &
                                - 11.12915953695671500_RK &
                                + 21.43217730985805500_RK * logzplus1 & 
                                - 33.75904206577289000_RK * logzplus1**2 & 
                                + 30.03916282499634200_RK * logzplus1**3 &
                                - 11.12086545981264500_RK * logzplus1**4 
        elseif (logzplus1>0.9282193027394269_RK .and. logzplus1<=1.3937663759585892_RK) then
            logBinaryMergerRate = &
                                - 1.802362223155230800_RK &
                                - 20.58526172567768200_RK * logzplus1 & 
                                + 38.93828966743146000_RK * logzplus1**2 & 
                                - 27.19863916580484500_RK * logzplus1**3 &
                                + 6.138928143113263000_RK * logzplus1**4
        elseif (logzplus1>1.3937663759585892_RK .and. logzplus1<=3.0411835364579027_RK) then
            logBinaryMergerRate = &
                                - 7.711815956299844000_RK &
                                + 11.68891993486079700_RK * logzplus1 & 
                                - 10.62908897824095300_RK * logzplus1**2 & 
                                + 2.945685425778300700_RK * logzplus1**3 &
                                - 0.327069839977957850_RK * logzplus1**4
        elseif (logzplus1<=0.02955880224154443_RK .or. logzplus1>3.0411835364579027_RK) then
            logBinaryMergerRate = 0._RK
        end if
    end function getLogBinaryMergerRateLognormM17

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module StarFormation_mod