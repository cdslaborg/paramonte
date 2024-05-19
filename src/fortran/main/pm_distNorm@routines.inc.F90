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
!>  This include file contains procedure implementations of [pm_distNorm](@ref pm_distNorm).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the RNG.
#if     setNormRand_ENABLED
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#else
#error  "Unrecognized interface."
#endif
#endif
        !%%%%%%%%%%%%%%%%%%%%
#if     getNormLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: mu_def, invSigma
        if (present(mu)) then
            mu_def = mu
        else
            mu_def = 0._RKG
        end if
        if (present(sigma)) then
            invSigma = 1._RKG / sigma
        else
            invSigma = 1._RKG
        end if
        call setNormLogPDF(logPDF, x, mu_def, invSigma, log(invSigma))

        !%%%%%%%%%%%%%%%%%%%%
#elif   setNormLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: PI = acos(-1._RKG)
        real(RKG)   , parameter :: INVERSE_SQRT_TWO_PI = 1._RKG / sqrt(2._RKG * PI)
        real(RKG)   , parameter :: LOG_INVERSE_SQRT_TWO_PI = log(INVERSE_SQRT_TWO_PI)
#if     DS_ENABLED || MS_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKG < invSigma, SK_"@setNormLogPDF(): The condition `0. < invSigma` must hold. invSigma = "//getStr(invSigma))
        CHECK_ASSERTION(__LINE__, abs(logInvSigma - log(invSigma)) <= epsilon(0._RKG)*100, SK_"@setNormLogPDF(): The condition `abs(logInvSigma - log(invSigma)) <= epsilon(0._RKG)*100` must hold. logInvSigma, log(invSigma) = "//getStr([logInvSigma, log(invSigma)]))
#endif
#if     DD_ENABLED
        logPDF = LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * x**2
#elif   DS_ENABLED
        logPDF = logInvSigma + LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * (x * invSigma)**2
#elif   MD_ENABLED
        logPDF = LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * (x - mu)**2
#elif   MS_ENABLED
        logPDF = logInvSigma + LOG_INVERSE_SQRT_TWO_PI - 0.5_RKG * ((x - mu) * invSigma)**2
#else
#error  "Unrecognized interface"
#endif

        !%%%%%%%%%%%%%%%%%
#elif   getNormCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        if (present(mu) .and. present(sigma)) then
            call setNormCDF(cdf, x, mu, 1._RKG / sigma)
        elseif (present(sigma)) then
            call setNormCDF(cdf, x, 0._RKG, 1._RKG / sigma)
        elseif (present(mu)) then
            call setNormCDF(cdf, x, mu)
        else
            call setNormCDF(cdf, x)
        end if

        !%%%%%%%%%%%%%%%%%
#elif   setNormCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: INVERSE_SQRT_TWO = 1._RKG / sqrt(2._RKG)
#if     DD_ENABLED
        cdf = 0.5_RKG * (1._RKG + erf(x * INVERSE_SQRT_TWO))
#elif   MD_ENABLED
        cdf = 0.5_RKG * (1._RKG + erf((x - mu) * INVERSE_SQRT_TWO))
#elif   MS_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKG < invSigma, SK_"@setNormCDF(): The condition `0. < invSigma` must hold. invSigma = "//getStr(invSigma))
        cdf = 0.5_RKG * (1._RKG + erf((x - mu) * invSigma * INVERSE_SQRT_TWO))
#else
#error  "Unrecognized interface"
#endif

        !%%%%%%%%%%%%%%%%%%
#elif   getNormQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%

        if (present(mu) .and. present(sigma)) then
            call setNormQuan(quantile, cdf, mu, sigma)
        elseif (present(sigma)) then
            call setNormQuan(quantile, cdf, 0._RKG, sigma)
        elseif (present(mu)) then
            call setNormQuan(quantile, cdf, mu)
        else
            call setNormQuan(quantile, cdf)
        end if

        !%%%%%%%%%%%%%%%%%%
#elif   setNormQuan_ENABLED
        !%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: SQRT_TWO = sqrt(2._RKG)
        if (0._RKG < cdf .and. cdf < 1._RKG) then
            call setErfInv(quantile, 2 * cdf - 1._RKG, abserr = 100 * epsilon(cdf))
#if         DD_ENABLED
            quantile = quantile * SQRT_TWO
#elif       MD_ENABLED
            quantile = quantile * SQRT_TWO + mu
#elif       MS_ENABLED
            CHECK_ASSERTION(__LINE__, 0._RKG < sigma, SK_"@setNormQuan(): The condition `0. < sigma` must hold. sigma = "//getStr(sigma))
            quantile = quantile * SQRT_TWO * sigma + mu
#else
#error      "Unrecognized interface"
#endif
        else
            quantile = huge(quantile)
            if (0._RKG == cdf) quantile = -quantile
            CHECK_ASSERTION(__LINE__, 0._RKG <= cdf .and. cdf <= 1._RKG, SK_"@setNormQuan(): The condition `0. <= cdf .and. cdf <= 1` must hold. cdf = "//getStr(cdf))
        end if

        !%%%%%%%%%%%%%%%%%
#elif   getZigNorm_ENABLED
        !%%%%%%%%%%%%%%%%%

        zig = getZig( nlay = nlay & ! LCOV_EXCL_LINE
                    , getFunc = getFuncNorm_RKG & ! LCOV_EXCL_LINE
                    , getFuncInv = getFuncInvNorm_RKG & ! LCOV_EXCL_LINE
                    , getZigArea = getZigAreaNorm_RKG & ! LCOV_EXCL_LINE
                    , abserr = abserr & ! LCOV_EXCL_LINE
                    , abstol = abstol & ! LCOV_EXCL_LINE
                    )

    contains

        PURE function getFuncNorm_RKG(x) result(func)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: func
            !func = getFuncNorm(x)
            call setNormLogPDF(func, x)
            func = exp(func)
        end function

        pure function getFuncInvNorm_RKG(func) result(x)
            real(RKG)   , parameter     :: LOGSQRT2PI = log(sqrt(2 * acos(-1._RKG)))
            real(RKG)   , intent(in)    :: func
            real(RKG)                   :: x
            x = sqrt(-2._RKG * (LOGSQRT2PI + log(func)))
            !x = getFuncInvNorm(func)
        end function

        PURE function getZigAreaNorm_RKG(r) result(area)
            real(RKG)   , intent(in)    :: r
            real(RKG)                   :: area
            CHECK_ASSERTION(__LINE__, 0._RKG <= r, SK_"@getZigAreaNorm(): The condition `0. <= r` must hold. r = "//getStr(r))
            !area = getZigAreaNorm(r)
            call setNormCDF(area, r)
            area = 1._RKG - area + r * getFuncNorm_RKG(r)
        end function

        !pure function getGradNorm_RKG(x) result(grad)
        !    real(RKG)   , intent(in)    :: x
        !    real(RKG)                   :: grad
        !    grad = getGradNorm(x)
        !end function

!        !%%%%%%%%%%%%%%%%%%
!#elif   getFuncNorm_ENABLED
!        !%%%%%%%%%%%%%%%%%%
!
!        real(RKG), parameter :: INVSQRT2PI = 1._RKG / sqrt(2 * acos(-1._RKG))
!        func = INVSQRT2PI * exp(-0.5_RKG * x**2)
!
!        !%%%%%%%%%%%%%%%%%%
!#elif   getGradNorm_ENABLED
!        !%%%%%%%%%%%%%%%%%%
!
!        grad = -x * getFuncNorm(x)
!
!        !%%%%%%%%%%%%%%%%%%%%%
!#elif   getFuncInvNorm_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%
!
!        x = sqrt(-2._RKG * log(func))
!
!        !%%%%%%%%%%%%%%%%%%%%%
!#elif   getZigAreaNorm_ENABLED
!        !%%%%%%%%%%%%%%%%%%%%%
!
!        real(RKG), parameter :: SQRT2PI = sqrt(2 * acos(-1._RKG))
!        CHECK_ASSERTION(__LINE__, 0._RKG <= r, SK_"@getZigAreaNorm(): The condition `0. <= r` must hold. r = "//getStr(r))
!        call setNormCDF(area, r)
!        area = r * getFuncNorm(r) + SQRT2PI * (1._RKG - area)

        !%%%%%%%%%%%%%%%%%%
#elif   getNormRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        call setNormRand(rand)
        if (present(std)) rand = rand * std
        rand = rand + mean

        !%%%%%%%%%%%%%%%%%%
#elif   setNormRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: index
        real(RKG)   :: dumm, func
#if     ZD_ENABLED
        real(RKG)   , parameter :: zig(1 : 2, 0 : 256) = real(ZIG_RKB, RKG)
        integer(IK) , parameter :: nlay = ubound(zig, 2, IK)
        CHECK_ASSERTION(__LINE__, precision(rand) <= ZIG_PRECISION, SK_"@setNormRand(): The condition `precision(rand) <= ZIG_PRECISION` must hold. precision(rand) <= ZIG_PRECISION = "//getStr([int(precision(rand), IK), ZIG_PRECISION]))
#elif   ZA_ENABLED
        integer(IK) :: nlay
        nlay = ubound(zig, 2, IK)
        CHECK_ASSERTION(__LINE__, all(0._RKG <= zig), SK_"@setNormRand(): The condition `all(0. <= zig)` must hold. zig = "//getStr(zig))
        CHECK_ASSERTION(__LINE__, size(zig, 1, IK) == 2_IK, SK_"@setNormRand(): The condition `size(zig, 1) == 2` must hold. shape(zig) = "//getStr(shape(zig, IK)))
#else
#error  "Unrecognized interface."
#endif
#if     D1_ENABLED
#define GET_RAND(i)rand(i)
        block
            integer(IK) :: irand
            loopOverRand: do irand = 1_IK, size(rand, 1, IK)
#elif           D0_ENABLED
#define         GET_RAND(i)rand
#else
#error          "Unrecognized interface."
#endif
                loopTry: do
                    call setUnifRand(RNG index, 0_IK, nlay - 1)
#if                 RNGD_ENABLED
#define             SET_RAND(X)call random_number(X)
                    call random_number(GET_RAND(irand))
                    !if (GET_RAND(irand) == 0 .or. GET_RAND(irand) == 1) print *, "GET_RAND(irand)", GET_RAND(irand)
                    GET_RAND(irand) = GET_RAND(irand) * 2 - 1._RKG
#elif               RNGF_ENABLED || RNGX_ENABLED
#define             SET_RAND(X)call setUnifRand(rng, X)
                    call setUnifRand(RNG GET_RAND(irand), -1._RKG, +1._RKG)
#else
#error              "Unrecognized interface."
#endif
                    GET_RAND(irand) = GET_RAND(irand) * zig(1, index)
                    if (zig(1, index + 1) < abs(GET_RAND(irand))) then ! rejection involved.
                        SET_RAND(dumm)
                        if (0_IK < index) then
                            dumm = zig(2, index) + dumm * (zig(2, index + 1) - zig(2, index)) ! height
                            call setNormLogPDF(func, GET_RAND(irand)); func = exp(func)
                            if (func < dumm) cycle loopTry
                        else
                            ! For the normal tail, we follow the method of Marsaglia, 1964:
                            ! generate x = -ln(U1) = r; y = -ln(U2) until y + y > x * x, then return r + x.
                            index = int(sign(1._RKG, GET_RAND(irand)), IK)
                            loopTail: do
                                !print *, "loopTail", dumm, zig(1)
                                GET_RAND(irand) = -log(1._RKG - dumm) / zig(1, 1)
                                SET_RAND(dumm)
                                !if (dumm == 0 .or. dumm == 1) print *, "dumm2", dumm
                                dumm = -log(1._RKG - dumm)
                                if (GET_RAND(irand)**2 < 2 * dumm) exit loopTail
                                SET_RAND(dumm)
                                !if (dumm == 0 .or. dumm == 1) print *, "dumm3", dumm
                            end do loopTail
                            GET_RAND(irand) = GET_RAND(irand) + zig(1, 1)
                            if (index < 0_IK) GET_RAND(irand) = -GET_RAND(irand)
                        end if
                    end if
                    exit loopTry
                end do loopTry
#if         D1_ENABLED
            end do loopOverRand
        end block
#endif
#undef  SET_RAND
#undef  GET_RAND

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getNormRandBox_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        ! Box-Muller rejection approach.
        real(RKG)           :: rand1
        real(RKG)   , save  :: rand2 = -1._RKG
        logical(LK) , save  :: failed = .true._LK
        if (failed) then
            do
                call random_number(rand1)
                call random_number(rand2)
                call setNormRandBox(rand1, rand2, failed)
                if (.not. failed) exit
            end do
            rand = rand1
        else
            rand = rand2
            failed = .true._LK
        end if
#if     MD_ENABLED
        rand = rand + mean
#elif   MS_ENABLED
        rand = rand * std + mean
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setNormRandBox_ENABLED && D0_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Perform checks.
#define CHECK_RAND1 \
CHECK_ASSERTION(__LINE__, 0._RKG <= rand1 .and. rand1 < 1._RKG, SK_"@setNormRand(): The condition `0. <= rand1 .and. rand1 < 1.` must hold. rand1 = "//getStr(rand1));
#define CHECK_RAND2 \
CHECK_ASSERTION(__LINE__, 0._RKG <= rand2 .and. rand2 < 1._RKG, SK_"@setNormRand(): The condition `0. <= rand2 .and. rand2 < 1.` must hold. rand2 = "//getStr(rand2));
        ! Define the factor.
#if     Basic_ENABLED
#define GET_RE(x) x%re
#define GET_IM(x) x%im
        complex(RKG) :: factor
        real(RKG), parameter :: TWP_PI = 2 * acos(-1._RKG), EPS = epsilon(0._RKG)
        CHECK_RAND1
        CHECK_RAND2
        factor = sqrt(-2._RKG * log(rand1 + EPS)) * exp(cmplx(0._RKG, -TWP_PI * rand2, RKG)) ! compute sin and cos, hopefully in one instruction.
#elif   Polar_ENABLED
#define GET_RE(x) rand1 * x
#define GET_IM(x) rand2 * x
        real(RKG) :: factor, rsq
        CHECK_RAND1
        CHECK_RAND2
        rand1 = 2._RKG * rand1 - 1._RKG
        rand2 = 2._RKG * rand2 - 1._RKG
        rsq = rand1**2 + rand2**2
        if (0._RKG < rsq .and. rsq < 1._RKG) then
            failed = .false._LK
            factor = sqrt(-2._RKG * log(rsq) / rsq)
        else
            failed = .true._LK
            return
        end if
#else
#error  "Unrecognized interface."
#endif
        ! Scale the numbers.
#if     MS_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKG < std, SK_"@setNormRand(): The condition `0. < std` must hold. std = "//getStr(std))
        rand1 = GET_RE(factor) * std + mean
        rand2 = GET_IM(factor) * std + mean
#elif   MD_ENABLED
        rand1 = GET_RE(factor) + mean
        rand2 = GET_IM(factor) + mean
#elif   DD_ENABLED
        rand1 = GET_RE(factor)
        rand2 = GET_IM(factor)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getNormEntropy_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        real(RKG), parameter :: LOG_TWO_PI = log(2 * acos(-1._RKG))
        entropy = 0.5_RKG + 0.5_RKG * (LOG_TWO_PI + logVar)

        !%%%%%%%%%%%%%%%%%%%%
#elif   getNormFisher_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKG <= varinv, SK_"@getNormFisher(): The input `varInv` must be positive. varInv = "//getStr(varInv))
        Fisher(1,1) = varInv
        Fisher(2,1) = 0._RKG
        Fisher(1,2) = 0._RKG
        Fisher(2,2) = 2._RKG * varInv

        !%%%%%%%%%%%%%%%%%
#elif   getNormKLD_ENABLED
        !%%%%%%%%%%%%%%%%%

#if     DV_ENABLED || MV_ENABLED
        real(RKG) :: varP2varQ
        CHECK_ASSERTION(__LINE__, 0._RKG < varP, SK_"@getNormKLD(): The condition `0. < varP` must hold. varP = "//getStr(varP))
        CHECK_ASSERTION(__LINE__, 0._RKG < varQ, SK_"@getNormKLD(): The condition `0. < varQ` must hold. varP = "//getStr(varQ))
#endif
#if     MD_ENABLED || MV_ENABLED
        CHECK_ASSERTION(__LINE__, 0._RKG <= meanDiffSq, SK_"@getNormKLD(): The condition `0. <= meanDiffSq` must hold. meanDiffSq = "//getStr(meanDiffSq))
#endif
#if     MV_ENABLED
        varP2varQ = varP / varQ
        kld = 0.5_RKG * (varP2varQ - log(varP2varQ) - 1._RKG + meanDiffSq / varQ)
#elif   MD_ENABLED
        kld = 0.5_RKG * meanDiffSq
#elif   DV_ENABLED
        varP2varQ = varP / varQ
        kld = 0.5_RKG * (varP2varQ - log(varP2varQ) - 1._RKG)
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  CHECK_RAND1
#undef  CHECK_RAND2
#undef  GET_RE
#undef  GET_IM
#undef  RNG