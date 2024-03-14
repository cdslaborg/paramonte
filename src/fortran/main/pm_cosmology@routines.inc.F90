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

!>  \brief
!>  This file contains implementations of procedures [pm_cosmology](@ref pm_cosmology).
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%
#if     getSizeUnivNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: sindex(1000), neval_def, nint_def, err_def
        real(RKC)   :: sinfo(4, 1000), abserr
        !reltol = sqrt(epsilon(0._RKC))
        if (zplus1 < huge(zplus1)) then
            err_def = getQuadErr( getFunc = getIntegrand & ! LCOV_EXCL_LINE
                                , lb = zplus1 & ! LCOV_EXCL_LINE
                                , ub = pinf & ! LCOV_EXCL_LINE
                                , abstol = 0._RKC & ! LCOV_EXCL_LINE
                                , reltol = reltol & ! LCOV_EXCL_LINE
                                , qrule = GK21 & ! LCOV_EXCL_LINE
                                , help = weps & ! LCOV_EXCL_LINE
                                , integral = sizeUnivNormed & ! LCOV_EXCL_LINE
                                , abserr = abserr & ! LCOV_EXCL_LINE
                                , sinfo = sinfo & ! LCOV_EXCL_LINE
                                , sindex = sindex & ! LCOV_EXCL_LINE
                                , neval = neval_def & ! LCOV_EXCL_LINE
                                , nint = nint_def & ! LCOV_EXCL_LINE
                                )
            if (present(err)) then
                err = err_def
            elseif (err_def /= 0_IK) then
                error stop SK_"@getQuadErr() failed with err = "//getStr(err_def) ! LCOV_EXCL_LINE
            end if
            if (present(neval)) neval = neval_def
        else
            sizeUnivNormed = 0._RKC
            if (present(err)) err = 0_IK
            if (present(neval)) neval = 0_IK
        end if

    contains

        function getIntegrand(zplus1) result(integrand)
            real(RKC)   , intent(in)    :: zplus1
            real(RKC)                   :: integrand
#if         Z_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1)))
#elif       ZML_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL)))
#elif       ZMLR_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)))
#elif       ZMLRK_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
#else
#error      "Unrecognized interface."
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisLookbackNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: sindex(1000), neval_def, nint_def, err_def
        real(RKC)   :: sinfo(4, 1000), abserr
        if (zplus1 > 1._RKC) then
            err_def = getQuadErr( getFunc = getIntegrand & ! LCOV_EXCL_LINE
                                , lb = 1._RKC & ! LCOV_EXCL_LINE
                                , ub = zplus1 & ! LCOV_EXCL_LINE
                                , abstol = 0._RKC & ! LCOV_EXCL_LINE
                                , reltol = reltol & ! LCOV_EXCL_LINE
                                , qrule = GK21 & ! LCOV_EXCL_LINE
                                , integral = disLookbackNormed & ! LCOV_EXCL_LINE
                                , abserr = abserr & ! LCOV_EXCL_LINE
                                , sinfo = sinfo & ! LCOV_EXCL_LINE
                                , sindex = sindex & ! LCOV_EXCL_LINE
                                , neval = neval_def & ! LCOV_EXCL_LINE
                                , nint = nint_def & ! LCOV_EXCL_LINE
                                )
            if (present(err)) then
                err = err_def
            elseif (err_def /= 0_IK) then
                error stop SK_"@getQuadErr() failed with err = "//getStr(err_def) ! LCOV_EXCL_LINE
            end if
            if (present(neval)) neval = neval_def
        else
            disLookbackNormed = 0._RKC
            if (present(err)) err = 0_IK
            if (present(neval)) neval = 0_IK
        end if

    contains

        function getIntegrand(zplus1) result(integrand)
            real(RKC)   , intent(in)    :: zplus1
            real(RKC)                   :: integrand
#if         Z_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1)))
#elif       ZML_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL)))
#elif       ZMLR_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR)))
#elif       ZMLRK_ENABLED
            integrand = 1._RKC / (zplus1 * sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK)))
#else
#error      "Unrecognized interface."
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisComNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: sindex(1000), neval_def, nint_def, err_def
        real(RKC)   :: sinfo(4, 1000), abserr
        !reltol = sqrt(epsilon(0._RKC))
        if (zplus1 > 1._RKC) then
            err_def = getQuadErr( getFunc = getHubbleParamNormedInv & ! LCOV_EXCL_LINE
                                , lb = 1._RKC & ! LCOV_EXCL_LINE
                                , ub = zplus1 & ! LCOV_EXCL_LINE
                                , abstol = 0._RKC & ! LCOV_EXCL_LINE
                                , reltol = reltol & ! LCOV_EXCL_LINE
                                , qrule = GK21 & ! LCOV_EXCL_LINE
                                , integral = disComNormed & ! LCOV_EXCL_LINE
                                , abserr = abserr & ! LCOV_EXCL_LINE
                                , sinfo = sinfo & ! LCOV_EXCL_LINE
                                , sindex = sindex & ! LCOV_EXCL_LINE
                                , neval = neval_def & ! LCOV_EXCL_LINE
                                , nint = nint_def & ! LCOV_EXCL_LINE
                                )
            if (present(err)) then
                err = err_def
            elseif (err_def /= 0_IK) then
                error stop SK_"@getQuadErr() failed with err = "//getStr(err_def) ! LCOV_EXCL_LINE
            end if
            if (present(neval)) neval = neval_def
        else
            disComNormed = 0._RKC
            if (present(err)) err = 0_IK
            if (present(neval)) neval = 0_IK
        end if

    contains

        function getHubbleParamNormedInv(zplus1) result(hubbleParamNormedInv)
            real(RKC)   , intent(in)    :: zplus1
            real(RKC)                   :: hubbleParamNormedInv
#if         Z_ENABLED
            hubbleParamNormedInv = 1._RKC / sqrt(getHubbleParamNormedSq(zplus1))
#elif       ZML_ENABLED
            hubbleParamNormedInv = 1._RKC / sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL))
#elif       ZMLR_ENABLED
            hubbleParamNormedInv = 1._RKC / sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR))
#elif       ZMLRK_ENABLED
            hubbleParamNormedInv = 1._RKC / sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK))
#else
#error      "Unrecognized interface."
#endif
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisComTransNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Z_ENABLED
        disComTransNormed = getDisComNormed(zplus1, reltol, neval, err)
#elif   ZML_ENABLED
        disComTransNormed = getDisComNormed(zplus1, omegaM, omegaL, reltol, neval, err)
#elif   ZMLR_ENABLED
        disComTransNormed = getDisComNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
#elif   ZMLRK_ENABLED
        real(RKC)   , parameter :: POS_EPS = epsilon(0._RKC)
        real(RKC)   , parameter :: NEG_EPS = -POS_EPS
        disComTransNormed = getDisComNormed(zplus1, omegaM, omegaL, omegaR, omegaK, reltol, neval, err)
        if (omegaK > POS_EPS) then
            disComTransNormed = sinh(disComTransNormed * sqrtAbsOmegaK) / sqrtAbsOmegaK
        elseif (omegaK < NEG_EPS) then
            disComTransNormed =  sin(disComTransNormed * sqrtAbsOmegaK) / sqrtAbsOmegaK
        end if
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisAngNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

#if     Z_ENABLED
        disAngNormed = getDisComTransNormed(zplus1, reltol, neval, err) / zplus1
#elif   ZML_ENABLED
        disAngNormed = getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err) / zplus1
#elif   ZMLR_ENABLED
        disAngNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) / zplus1
#elif   ZMLRK_ENABLED
        disAngNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) / zplus1
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisLumNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

#if     Z_ENABLED
        disLumNormed = getDisComTransNormed(zplus1, reltol, neval, err) * zplus1
#elif   ZML_ENABLED
        disLumNormed = getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err) * zplus1
#elif   ZMLR_ENABLED
        disLumNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err) * zplus1
#elif   ZMLRK_ENABLED
        disLumNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err) * zplus1
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisComTransNormedWU10_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Z_ENABLED
        real(RKC), parameter :: omegaL = real(OMEGA_L, RKC)
        real(RKC), parameter :: omegaM = real(OMEGA_M, RKC)
#define SET_PARAMETER(x) parameter(x)
#elif   ZML_ENABLED
#define SET_PARAMETER(x) x
#else
#error  "Unrecognized interface."
#endif
        real(RKC) :: alpha1, x1, x1Sq, psix1
        real(RKC), parameter :: PSI_COEF1 = 2._RKC**(2._RKC/3._RKC)
        real(RKC), parameter :: PSI_COEF2 = -PSI_COEF1 / 252._RKC
        real(RKC), parameter :: PSI_COEF3 = +PSI_COEF1 / 21060._RKC
        real(RKC), parameter :: ONE_THIRD = 1._RKC / 3._RKC
        real(RKC), parameter :: ONE_SIXTH = 1._RKC / 6._RKC
        real(RKC) :: twiceOmegaL2OmegaM
        real(RKC) :: alpha0
        real(RKC) :: x0
        real(RKC) :: x0Sq
        real(RKC) :: psiX0
        SET_PARAMETER(twiceOmegaL2OmegaM = 2._RKC * real(omegaL / omegaM, RKC)) ! fpp
        SET_PARAMETER(alpha0 = 1._RKC + twiceOmegaL2OmegaM) ! fpp
        SET_PARAMETER(x0 = log(alpha0 + sqrt(alpha0**2 - 1._RKC))) ! fpp
        SET_PARAMETER(x0Sq = x0**2) ! fpp
        SET_PARAMETER(psiX0 = x0**ONE_THIRD * (PSI_COEF1 + x0Sq * (PSI_COEF2 + x0Sq * PSI_COEF3))) ! fpp
        CHECK_ASSERTION(__LINE__, 1._RKC <= zplus1, SK_"@getlogDisLum(): The condition `zplus1 >= 0.` must hold. zplus1 = "//getStr(zplus1)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(1._RKC - omegaM - omegaL) <= epsilon(0._RKC), SK_"@getlogDisLum(): The condition `1._RKC - omegaM - omegaL == 0._RKC` must hold. omegaM, omegaL = "//getStr([omegaM, omegaL])) ! fpp
        alpha1  = 1._RKC + twiceOmegaL2OmegaM / zplus1**3
        x1      = log(alpha1 + sqrt(alpha1**2 - 1._RKC))
        x1Sq    = x1**2
        psix1   = x1**ONE_THIRD * (PSI_COEF1 + x1Sq * (PSI_COEF2 + x1Sq * PSI_COEF3))
        disComTransNormedWU10  = (psiX0 - psix1) / real(omegaL**ONE_SIXTH * omegaM**ONE_THIRD, RKC)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getlogDisLookback_ENABLED || getlogDisLookback_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#error  "Unrecognized interface."

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getHubbleParamNormedSq_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Z_ENABLED
        real(RKC)   , parameter :: omegaL = real(OMEGA_L, RKC)
        real(RKC)   , parameter :: omegaM = real(OMEGA_M, RKC)
        hubbleParamNormedSq = omegaL + zplus1**3 * omegaM
#elif   ZML_ENABLED
        CHECK_ASSERTION(__LINE__, 1._RKC <= zplus1, SK_"@getHubbleParamNormedSq(): The condition `1._RKC <= zplus1` must hold. zplus1 = "//getStr(zplus1)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(1._RKC - omegaM - omegaL) < EPS, SK_"@getHubbleParamNormedSq(): The condition `omegaM + omegaL = 1._RKC` must hold. omegaM, omegaL = "//getStr([omegaM, omegaL])) ! fpp
        hubbleParamNormedSq = omegaL + zplus1**3 * omegaM
#elif   ZMLR_ENABLED
        CHECK_ASSERTION(__LINE__, 1._RKC <= zplus1, SK_"@getHubbleParamNormedSq(): The condition `1._RKC <= zplus1` must hold. zplus1 = "//getStr(zplus1)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(1._RKC - omegaM - omegaL - omegaR) < EPS, SK_"@getHubbleParamNormedSq(): The condition `omegaM + omegaL + omegaR = 1._RKC` must hold. omegaM, omegaL, omegaR = "//getStr([omegaM, omegaL, omegaR])) ! fpp
        hubbleParamNormedSq = omegaL + zplus1**2 * (zplus1 * (omegaM + zplus1 * omegaR))
#elif   ZMLRK_ENABLED
        CHECK_ASSERTION(__LINE__, 1._RKC <= zplus1, SK_"@getHubbleParamNormedSq(): The condition `1._RKC <= zplus1` must hold. zplus1 = "//getStr(zplus1)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(1._RKC - omegaM - omegaL - omegaR - omegaK) < EPS, SK_"@getHubbleParamNormedSq(): The condition `omegaM + omegaL + omegaR + omegaK = 1._RKC` must hold. omegaM, omegaL, omegaR, omegaK = "//getStr([omegaM, omegaL, omegaR, omegaK])) ! fpp
        hubbleParamNormedSq = omegaL + zplus1**2 * (omegaK + zplus1 * (omegaM + zplus1 * omegaR))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getVolComDiffNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Z_ENABLED
        volComDiffNormed = getDisComTransNormed(zplus1, reltol, neval, err)
        volComDiffNormed = volComDiffNormed**2 / sqrt(getHubbleParamNormedSq(zplus1))
#elif   ZML_ENABLED
        volComDiffNormed = getDisComTransNormed(zplus1, omegaM, omegaL, reltol, neval, err)
        volComDiffNormed = volComDiffNormed**2 / sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL))
#elif   ZMLR_ENABLED
        volComDiffNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, reltol, neval, err)
        volComDiffNormed = volComDiffNormed**2 / sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR))
#elif   ZMLRK_ENABLED
        volComDiffNormed = getDisComTransNormed(zplus1, omegaM, omegaL, omegaR, omegaK, sqrtAbsOmegaK, reltol, neval, err)
        volComDiffNormed = volComDiffNormed**2 / sqrt(getHubbleParamNormedSq(zplus1, omegaM, omegaL, omegaR, omegaK))
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVolComDiffNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%

        volComDiffNormed = disComTransNormedSq / hubbleParamNormed

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   getVolComNormed_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

#if     D_ENABLED
        real(RKC), parameter :: COEF = 4._RKC * acos(-1._RKC) / 3._RKC
        volComNormed = COEF * disComTransNormed**3
#elif   DOS_ENABLED
        real(RKC), parameter :: COEF = 2._RKC * acos(-1._RKC)
        if (omegaK > 0._RKC) then
            volComNormed = COEF * (disComTransNormed * sqrt(1._RKC + omegaK * disComTransNormed**2) - asinh(sqrtAbsOmegaK * disComTransNormed) / sqrtAbsOmegaK) / omegaK
        elseif (omegaK < 0._RKC) then
            volComNormed = COEF * (disComTransNormed * sqrt(1._RKC + omegaK * disComTransNormed**2) - asin (sqrtAbsOmegaK * disComTransNormed) / sqrtAbsOmegaK) / omegaK
        else
            volComNormed = getVolComNormed(disComTransNormed)
        end if
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  SET_PARAMETER
!        use pm_cosmology, only: LOG_LIGHT_SPEED
!#if     getLogDVDZ_ENABLED || getLogDVDZ_D1_ENABLED
!        use pm_cosmology, only: omegaM => OMEGA_M
!        use pm_cosmology, only: logH0 => LOG_HUBBLE_CONST
!#define SET_PARAMETER(x) parameter(x)
!#elif   getLogDVDZ_HM_ENABLED || getLogDVDZ_HM_D1_ENABLED
!#define SET_PARAMETER(x) x
!#else
!#error  "Unrecognized interface."
!#endif
!
!#if     getLogDVDZ_D1_ENABLED || getLogDVDZ_HM_D1_ENABLED
!        integer(IK) :: i
!#define GET_ELEMENT(Array) Array(i)
!#elif   getLogDVDZ_ENABLED || getLogDVDZ_HM_ENABLED
!#define GET_ELEMENT(Array) Array
!#define ALL
!#else
!#error  "Unrecognized interface."
!#endif
!
!        real(RKC) :: logCoef, omegaDE
!        SET_PARAMETER(omegaDE = 1._RKC - real(omegaM,RKC))
!        SET_PARAMETER(logCoef = log(4._RKC * acos(-1._RKC)) + real(LOG_LIGHT_SPEED,RKC) - real(logH0,RKC))
!
!#if     CHECK_ENABLED
!        block
!            use pm_err, only: setAsserted
!            use pm_val2str, only: getStr
!            intrinsic :: size
!            call setAsserted(ALL(zplus1 >= 1._RKC), MODULE_NAME//SK_"@getLogDVDZ(): `zplus1 >= 0.` must hold. zplus1 = "//getStr(zplus1))
!            call setAsserted(ALL(logzplus1 >= 0._RKC), MODULE_NAME//SK_"@getLogDVDZ(): `logzplus1 >= 0.` must hold. logzplus1 = "//getStr(logzplus1))
!#if         getLogDVDZ_D1_ENABLED || getLogDVDZ_HM_D1_ENABLED
!            call setAsserted(size(zplus1,kind=IK) == size(logzplus1,kind=IK), MODULE_NAME//SK_"@getLogDVDZ(): `size(zplus1) == size(logzplus1)` must hold. size(zplus1), size(logzplus1) = "//getStr([size(zplus1,kind=IK), size(logzplus1,kind=IK)]))
!            call setAsserted(size(zplus1,kind=IK) == size(TwiceLogLumDisMpc,kind=IK), MODULE_NAME//SK_"@getLogDVDZ(): `size(zplus1) == size(TwiceLogLumDisMpc)` must hold. size(zplus1), size(TwiceLogLumDisMpc) = "//getStr([size(zplus1,kind=IK), size(TwiceLogLumDisMpc,kind=IK)]))
!#endif
!        end block
!#endif
!
!
!#if     getLogDVDZ_D1_ENABLED || getLogDVDZ_HM_D1_ENABLED
!        do concurrent (i = 1 : size(zplus1, kind = IK))
!#endif
!            GET_ELEMENT(LogDVDZ) = logCoef + GET_ELEMENT(TwiceLogLumDisMpc) - 3 * GET_ELEMENT(logzplus1) - 0.5_RKC * log(real(omegaM,RKC) * GET_ELEMENT(zplus1)**3 + omegaDE)
!#if     getLogDVDZ_D1_ENABLED || getLogDVDZ_HM_D1_ENABLED
!        end do
!#endif
!
!#undef  SET_PARAMETER
!#undef  GET_ELEMENT
!#undef  ALL
