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
        CHECK_ASSERTION(__LINE__, 0._RKG < lambda, SK_"@setPoisLogPMF(): The condition `0. < lambda` must hold. lambda = "//getStr(lambda)) ! fpp
#if     Log_ENABLED
        CHECK_ASSERTION(__LINE__, abs(logLambda - log(lambda)) < epsilon(0._RKG) * 10, \
        SK_"@setPoisLogPMF(): The condition `logLambda == log(lambda)` must hold. logLambda, log(lambda) = "//getStr([logLambda, log(lambda)])) ! fpp
#elif   Def_ENABLED
#define LOGLAMBDA log(lambda)
#else
#error  "Unrecognized interface."
#endif
        ! Compute the PMF.
        logPMF = real(count, RKG) * LOGLAMBDA - lambda - log_gamma(real(count, RKG) + 1._RKG)

        !%%%%%%%%%%%%%%%%%
#elif   getPoisCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        real(RKG) :: countP1
        CHECK_ASSERTION(__LINE__, 0_IK <= count, SK_"@getPoisCDF(): The condition `0 <= count` must hold. count = "//getStr(count)) ! fpp
        countP1 = real(count, RKG) + 1._RKG
        call setPoisCDF(cdf, countP1, lambda, info)
        !call setPoisCDF(cdf, countP1, log_gamma(countP1), lambda, info)
        if (info < 0_IK) error stop MODULE_NAME//SK_"@getPoisCDF(): Call to setPoisCDF failed."

        !%%%%%%%%%%%%%%%%%
#elif   setPoisCDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        real(RKG) :: gaminclow
        ! Validate the input.
        CHECK_ASSERTION(__LINE__, 0._RKG < lambda, SK_"@setPoisCDF(): The condition `0. < lambda` must hold. lambda = "//getStr(lambda)) ! fpp
        CHECK_ASSERTION(__LINE__, 1._RKG <= countP1, SK_"@setPoisCDF(): The condition `1._RKG <= countP1` must hold. countP1 = "//getStr(countP1)) ! fpp
        !check_assertion(__LINE__, abs(logGammaCountP1 - log_gamma(countP1)) <= epsilon(1._RKG) * 100, \
        !SK_"@setPoisCDF(): The condition `logGammaCountP1 == log(countP1)` must hold. logGammaCountP1, log_gamma(countP1) = "\
        !//getStr([logGammaCountP1, log_gamma(countP1)])) ! fpp
        CHECK_ASSERTION(__LINE__, mod(countP1, 1._RKG) == 0._RKG, \
        SK_"@setPoisCDF(): The condition `mod(countP1, 1._RKG) == 0._RKG` must hold. countP1, mod(countP1, 1._RKG) = "\
        //getStr([countP1, mod(countP1, 1._RKG)])) ! fpp
        !call setGammaIncUpp(cdf, x = lambda, logGammaKappa = logGammaCountP1, kappa = countP1, info = info, tol = tol)
        call setGammaInc(gaminclow, cdf, x = lambda, kappa = countP1, info = info)

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
        real(RKG), parameter :: LAMBDA_LIMIT_RKG = real(LAMBDA_LIMIT, RKG)
#if     Exp_ENABLED
        ! Use only when lambda < 10.
        real(RKG) :: prod, unifrnd
        CHECK_ASSERTION(__LINE__, expNegLambda < 1._RKG, \
        SK_"@setPoisRand(): The condition `expNegLambda < 1.` must hold. expNegLambda = "//getStr(expNegLambda)) ! fpp
        CHECK_ASSERTION(__LINE__, exp(-LAMBDA_LIMIT_RKG) < expNegLambda, \
        SK_"@setPoisRand(): The condition `exp(-LAMBDA_LIMIT) < expNegLambda` must hold. exp(-LAMBDA_LIMIT), expNegLambda = "\
        //getStr([exp(-LAMBDA_LIMIT_RKG), expNegLambda])) ! fpp
#if     D1_ENABLED
        do irand = 1_IK, size(rand, 1, IK)
#endif
            GET_RAND(irand) = 0_IK
            prod = 1._RKG
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
        real(RKB)   , parameter :: UR = .43_RKG
        real(RKG)               :: unifrnd, rndunif, afac, bfac, invAlpha, vr, us, randr
        CHECK_ASSERTION(__LINE__, LAMBDA_LIMIT_RKG <= lambda, \
        SK_"@setPoisRand(): The condition `LAMBDA_LIMIT <= lambda` must hold. lambda = "\
        //getStr([LAMBDA_LIMIT_RKG, lambda])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(logLambda - log(lambda)) < epsilon(lambda) * 100, \
        SK_"@setPoisRand(): The condition `logLambda == sqrt(lambda)` must hold. logLambda, log(lambda) = "\
        //getStr([logLambda, log(lambda)])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(sqrtLambda - sqrt(lambda)) < epsilon(lambda) * 100, \
        SK_"@setPoisRand(): The condition `sqrtLambda == sqrt(lambda)` must hold. sqrtLambda, sqrt(lambda) = "\
        //getStr([sqrtLambda, sqrt(lambda)])) ! fpp
        bfac = 0.931_RKG + 2.53_RKG * sqrtLambda
        afac = -0.059_RKG + 0.02483_RKG * bfac
        invAlpha = 1.1239 + 1.1328_RKG / (bfac - 3.4_RKG)
        vr = 0.9277_RKG - 3.6224_RKG / (bfac - 2._RKG)
#if     D1_ENABLED
        do irand = 1_IK, size(rand, 1, IK)
#endif
            do
                call setUnifRand(RNG unifrnd)
                unifrnd = unifrnd - 0.5_RKG
                call setUnifRand(RNG rndunif)
                us = 0.5_RKG - abs(unifrnd)
                GET_RAND(irand) = floor((2 * afac / us + bfac) * unifrnd + lambda + UR, kind = IK)
                if (0.07_RKG <= us .and. rndunif <= vr) exit
                if (GET_RAND(irand) < 0_IK .or. (us < 0.013_RKG .and. us < rndunif)) cycle
                randr = real(GET_RAND(irand), RKG)
                if (log(rndunif) + log(invAlpha) - log(afac / (us * us) + bfac) <= randr * logLambda - lambda - log_gamma(randr + 1._RKG)) exit
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