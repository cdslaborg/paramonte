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
!>  This include file contains the implementation of procedures in [pm_distGeom](@ref pm_distGeom).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%
#if     getGeomLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKC < probSuccess, SK_"@setGeomLogPMF(): The condition `0 <= probSuccess` must hold. probSuccess = "//getStr(probSuccess)) ! fpp
        CHECK_ASSERTION(__LINE__, probSuccess <= 1._RKC, SK_"@setGeomLogPMF(): The condition `probSuccess <= 1` must hold. probSuccess = "//getStr(probSuccess)) ! fpp
        call setGeomLogPMF(logPMF, stepSuccess, log(probSuccess))

        !%%%%%%%%%%%%%%%%%%%%
#elif   setGeomLogPMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        ! Validate the input.
        CHECK_ASSERTION(__LINE__, 0_IK < stepSuccess, SK_"@setGeomLogPMF(): The condition `0 < stepSuccess` must hold. stepSuccess = "//getStr(stepSuccess)) ! fpp
        CHECK_ASSERTION(__LINE__, logProbSuccess <= 0._RKC, SK_"@setGeomLogPMF(): The condition `logProbSuccess <= 0.` must hold. logProbSuccess = "//getStr(logProbSuccess)) ! fpp
#if     Log_ENABLED
        CHECK_ASSERTION(__LINE__, abs(1._RKC - exp(logProbSuccess) - exp(logProbFailure)) < epsilon(0._RKC) * 100, \
        SK_"@setGeomLogPMF(): The condition `exp(logProbFailure) + exp(logProbSuccess) == 1.` must hold. logProbFailure, logProbSuccess = "\
        //getStr([logProbFailure, logProbSuccess])) ! fpp
#elif   Def_ENABLED
#define logProbFailure log(get1mexp(logProbSuccess))
#else
#error  "Unrecognized interface."
#endif
        ! Compute the PMF.
        if (logProbSuccess /= 0._RKC) then ! imperfect probability of success.
            logPMF = (real(stepSuccess, RKC) - 1_IK) * logProbFailure + logProbSuccess
        else ! 100% probability of success.
            if (1_IK < stepSuccess) then
                logPMF = -huge(logPMF)
            else
                logPMF = 0._RKC
            end if
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGeomCDF_ENABLED || setGeomCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Validate the input.
        CHECK_ASSERTION(__LINE__, 0._RKC < probSuccess, SK_"@setGeomCDF(): The condition `0. < probSuccess` must hold. probSuccess = "//getStr(probSuccess)) ! fpp
        CHECK_ASSERTION(__LINE__, probSuccess <= 1._RKC, SK_"@setGeomCDF(): The condition `probSuccess <= 1.` must hold. probSuccess = "//getStr(probSuccess)) ! fpp
        CHECK_ASSERTION(__LINE__, 0_IK < stepSuccess, SK_"@setGeomCDF(): The condition `0 <= stepSuccess` must hold. stepSuccess = "//getStr(stepSuccess)) ! fpp
        cdf = 1._RKC - (1._RKC - probSuccess)**stepSuccess

        !%%%%%%%%%%%%%%%%%%
#elif   getGeomRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKC < probSuccess, SK_"@getGeomRand(): The condition `0 <= probSuccess` must hold. probSuccess = "//getStr(probSuccess)) ! fpp
        if (probSuccess < 1._RKC) then
            call setGeomRand(rand, log(1._RKC - probSuccess))
        else
            rand = 1_IK
            CHECK_ASSERTION(__LINE__, probSuccess <= 1._RKC, SK_"@getGeomRand(): The condition `probSuccess <= 1` must hold. probSuccess = "//getStr(probSuccess)) ! fpp
        end if

        !%%%%%%%%%%%%%%%%%%
#elif   setGeomRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     D0_ENABLED
        real(RKC) :: unifrnd
#elif   D1_ENABLED
        real(RKC) :: unifrnd(size(rand, 1, IK)), logProbFailureInverse
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, logProbFailure < 0._RKC, SK_"@setGeomRand(): The condition `logProbFailure < 0.` must hold. logProbFailure = "//getStr(logProbFailure)) ! fpp
#if     RNGD_ENABLED || RNGF_ENABLED
        call random_number(unifrnd)
#elif   RNGX_ENABLED
        call setUnifRand(rng, unifrnd)
#else
#error  "Unrecognized interface."
#endif
#if     D0_ENABLED
        rand = 1_IK + floor(log(1._RKC - unifrnd) / logProbFailure)
#elif   D1_ENABLED
        logProbFailureInverse = 1._RKC / logProbFailure
        rand = 1_IK + floor(log(1._RKC - unifrnd) * logProbFailureInverse)
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  logProbFailure