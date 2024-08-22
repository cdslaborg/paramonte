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
!>  This include file contains the implementation of procedures in [pm_distBinom](@ref pm_distBinom).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%
#if     getBinomLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        call setBinomLogPMF(logPMF, nsuc, ntry, log(psuc), log(1._TKG - psuc))

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setBinomLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        ! Validate the input.
        CHECK_ASSERTION(__LINE__, 0 <= nsuc, SK_"@setBinomLogPMF(): The condition `0 <= nsuc` must hold. nsuc = "//getStr(nsuc)) ! fpp
        CHECK_ASSERTION(__LINE__, 0 <= ntry, SK_"@setBinomLogPMF(): The condition `0 <= ntry` must hold. ntry = "//getStr(ntry)) ! fpp
        CHECK_ASSERTION(__LINE__, nsuc <= ntry, SK_"@setBinomLogPMF(): The condition `nsuc <= ntry` must hold. nsuc, ntry = "//getStr([nsuc, ntry])) ! fpp
        CHECK_ASSERTION(__LINE__, abs(1 - exp(logp) - exp(logq)) <= epsilon(logp), SK_"@setBinomLogPMF(): The condition `abs(1 - exp(logp) - exp(logq)) <= epsilon(logp)` must hold. logp, logq = "//getStr([logp, logq])) ! fpp
        CHECK_ASSERTION(__LINE__, logp <= 0, SK_"@setBinomLogPMF(): The condition `logp <= 0` must hold. logp = "//getStr(logp)) ! fpp
        CHECK_ASSERTION(__LINE__, logq <= 0, SK_"@setBinomLogPMF(): The condition `logq <= 0` must hold. logq = "//getStr(logq)) ! fpp
        logPMF = log_gamma(real(ntry + 1_IK, TKG)) - log_gamma(real(nsuc + 1_IK, TKG)) - log_gamma(real(ntry - nsuc + 1_IK, TKG))
        logPMF = logPMF + nsuc * logp + (ntry - nsuc) * logq

        !%%%%%%%%%%%%%%%%%%
#elif   getBinomCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        call setBinomCDF(cdf, nsuc, ntry, psuc, info)
        if (info /= 0_IK) error stop MODULE_NAME//SK_"@getBinomCDF(): Call to setBinomCDF() failed."

        !%%%%%%%%%%%%%%%%%%
#elif   setBinomCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%

        real(TKG) :: alpha, beta
        ! Validate the input.
        CHECK_ASSERTION(__LINE__, 0 <= nsuc, SK_"@setBinomCDF(): The condition `0 <= nsuc` must hold. nsuc = "//getStr(nsuc)) ! fpp
        CHECK_ASSERTION(__LINE__, 0 <= ntry, SK_"@setBinomCDF(): The condition `0 <= ntry` must hold. ntry = "//getStr(ntry)) ! fpp
        CHECK_ASSERTION(__LINE__, nsuc <= ntry, SK_"@setBinomCDF(): The condition `nsuc <= ntry` must hold. nsuc, ntry = "//getStr([nsuc, ntry])) ! fpp
        CHECK_ASSERTION(__LINE__, 0 <= psuc .and. psuc <= 1, SK_"@setBinomCDF(): The condition `0 <= psuc .and. psuc <= 1` must hold. psuc = "//getStr(psuc)) ! fpp
        beta = real(1_IK + nsuc, TKG)
        alpha = real(ntry - nsuc, TKG)
        call setBetaInc ( cdf, x = 1._TKG - psuc &
                        , alpha = alpha &
                        , beta = beta &
                        , logFuncBeta = getLogBeta(alpha = alpha, beta = beta) &
                        , signed = .false._LK &
                        , info = info &
                        )

        !%%%%%%%%%%%%%%%%%%%
#elif   getBinomRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%
#elif   setBinomRand_ENABLED
        !%%%%%%%%%%%%%%%%%%%

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

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_RAND
#undef  RNG