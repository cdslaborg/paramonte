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
!>  This include file contains the implementation of procedures in [pm_distPois](@ref pm_distPois).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%
#if     getPoisLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        call setPoisLogPMF(logPMF, count, lambda)

        !%%%%%%%%%%%%%%%%%%%%
#elif   setPoisLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        ! Validate the input.
        CHECK_ASSERTION(__LINE__, 0_IK <= count, SK_"@setPoisLogPMF(): The condition `0 < count` must hold. count = "//getStr(count)) ! fpp
        CHECK_ASSERTION(__LINE__, 0._RKC < lambda, SK_"@setPoisLogPMF(): The condition `0. < lambda` must hold. lambda = "//getStr(lambda)) ! fpp
#if     Log_ENABLED
        CHECK_ASSERTION(__LINE__, abs(logLambda - log(lambda)) < epsilon(0._RKC) * 10, \
        SK_"@setPoisLogPMF(): The condition `logLambda == log(lambda)` must hold. logLambda, log(lambda) = "//getStr([logLambda, log(lambda)])) ! fpp
#elif   Def_ENABLED
#define LOGLAMBDA log(lambda)
#else
#error  "Unrecognized interface."
#endif
        ! Compute the PMF.
        logPMF = real(count, RKC) * LOGLAMBDA - lambda - log_gamma(real(count, RKC) + 1._RKC)

        !%%%%%%%%%%%%%%%%%
#elif   getPoisCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        real(RKC) :: countP1
        CHECK_ASSERTION(__LINE__, 0_IK <= count, SK_"@getPoisCDF(): The condition `0 <= count` must hold. count = "//getStr(count)) ! fpp
        countP1 = real(count, RKC) + 1._RKC
        call setPoisCDF(cdf, countP1, log_gamma(countP1), lambda, info)
        if (info < 0_IK) error stop MODULE_NAME//SK_"@getPoisCDF(): Call to setPoisCDF failed."

        !%%%%%%%%%%%%%%%%%
#elif   setPoisCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        ! Validate the input.
        CHECK_ASSERTION(__LINE__, 0._RKC < lambda, SK_"@setPoisCDF(): The condition `0. < lambda` must hold. lambda = "//getStr(lambda)) ! fpp
        CHECK_ASSERTION(__LINE__, 1._RKC <= countP1, SK_"@setPoisCDF(): The condition `1._RKC <= countP1` must hold. countP1 = "//getStr(countP1)) ! fpp
        CHECK_ASSERTION(__LINE__, abs(logGammaCountP1 - log_gamma(countP1)) <= epsilon(1._RKC) * 100, \
        SK_"@setPoisCDF(): The condition `logGammaCountP1 == log(countP1)` must hold. logGammaCountP1, log_gamma(countP1) = "\
        //getStr([logGammaCountP1, log_gamma(countP1)])) ! fpp
        CHECK_ASSERTION(__LINE__, mod(countP1, 1._RKC) == 0._RKC, \
        SK_"@setPoisCDF(): The condition `mod(countP1, 1._RKC) == 0._RKC` must hold. countP1, mod(countP1, 1._RKC) = "\
        //getStr([countP1, mod(countP1, 1._RKC)])) ! fpp
        call setGammaIncUpp(cdf, x = lambda, logGammaKappa = logGammaCountP1, kappa = countP1, info = info, tol = tol)

        !%%%%%%%%%%%%%%%%%%
#elif   getPoisRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        if (lambda < 10) then
            call setPoisRand(rand, exp(-lambda))
        else
            call setPoisRand(rand, lambda, log(lambda), sqrt(lambda))
        end if

        !%%%%%%%%%%%%%%%%%%
#elif   setPoisRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        ! Set the URNG.
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#else
#error  "Unrecognized interface."
#endif
        ! Set the dimension of `rand`.
#if     D0_ENABLED
#define GET_RAND(i) rand
#elif   D1_ENABLED
#define GET_RAND(i) rand(i)
        integer(IK) :: irand
#else
#error  "Unrecognized interface."
#endif
        real(RKC), parameter :: LAMBDA_LIMIT_RKC = real(LAMBDA_LIMIT, RKC)
#if     Exp_ENABLED
        ! Use only when lambda < 10.
        real(RKC) :: prod, unifrnd
        CHECK_ASSERTION(__LINE__, expNegLambda < 1._RKC, \
        SK_"@setPoisRand(): The condition `expNegLambda < 1.` must hold. expNegLambda = "//getStr(expNegLambda)) ! fpp
        CHECK_ASSERTION(__LINE__, exp(-LAMBDA_LIMIT_RKC) < expNegLambda, \
        SK_"@setPoisRand(): The condition `exp(-LAMBDA_LIMIT) < expNegLambda` must hold. exp(-LAMBDA_LIMIT), expNegLambda = "\
        //getStr([exp(-LAMBDA_LIMIT_RKC), expNegLambda])) ! fpp
#if     D1_ENABLED
        do irand = 1_IK, size(rand, 1, IK)
#endif
            GET_RAND(irand) = 0_IK
            prod = 1._RKC
            do
                call setUnifRand(RNG unifrnd)
                prod = prod * unifrnd
                if (prod <= expNegLambda) exit
                GET_RAND(irand) = GET_RAND(irand) + 1_IK
            end do
#if     D1_ENABLED
        end do
#endif
#elif   Rej_ENABLED
        real(RKB)   , parameter :: UR = .43_RKC
        real(RKC)               :: unifrnd, rndunif, afac, bfac, invAlpha, vr, us, randr
        CHECK_ASSERTION(__LINE__, LAMBDA_LIMIT_RKC <= lambda, \
        SK_"@setPoisRand(): The condition `LAMBDA_LIMIT <= lambda` must hold. lambda = "\
        //getStr([LAMBDA_LIMIT_RKC, lambda])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(logLambda - log(lambda)) < epsilon(lambda) * 100, \
        SK_"@setPoisRand(): The condition `logLambda == sqrt(lambda)` must hold. logLambda, log(lambda) = "\
        //getStr([logLambda, log(lambda)])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(sqrtLambda - sqrt(lambda)) < epsilon(lambda) * 100, \
        SK_"@setPoisRand(): The condition `sqrtLambda == sqrt(lambda)` must hold. sqrtLambda, sqrt(lambda) = "\
        //getStr([sqrtLambda, sqrt(lambda)])) ! fpp
        bfac = 0.931_RKC + 2.53_RKC * sqrtLambda
        afac = -0.059_RKC + 0.02483_RKC * bfac
        invAlpha = 1.1239 + 1.1328_RKC / (bfac - 3.4_RKC)
        vr = 0.9277_RKC - 3.6224_RKC / (bfac - 2._RKC)
#if     D1_ENABLED
        do irand = 1_IK, size(rand, 1, IK)
#endif
            do
                call setUnifRand(RNG unifrnd)
                unifrnd = unifrnd - 0.5_RKC
                call setUnifRand(RNG rndunif)
                us = 0.5_RKC - abs(unifrnd)
                GET_RAND(irand) = floor((2 * afac / us + bfac) * unifrnd + lambda + UR, kind = IK)
                if (0.07_RKC <= us .and. rndunif <= vr) exit
                if (GET_RAND(irand) < 0_IK .or. (us < 0.013_RKC .and. us < rndunif)) cycle
                randr = real(GET_RAND(irand), RKC)
                if (log(rndunif) + log(invAlpha) - log(afac / (us * us) + bfac) <= randr * logLambda - lambda - log_gamma(randr + 1._RKC)) exit
            end do
#if     D1_ENABLED
        end do
#endif
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_RAND
#undef  RNG