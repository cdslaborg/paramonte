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
!>  This include file contains the implementation of procedures in [pm_distBeta](@ref pm_distBeta).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%
#if     getBetaPDF_ENABLED
        !%%%%%%%%%%%%%%%%%

        if (x == 0._RKC) then
            if (alpha < 1._RKC) then
                pdf = +huge(pdf)
            else
                pdf = 0._RKC
            end if
        elseif (x == 1._RKC) then
            if (beta < 1._RKC) then
                pdf = +huge(pdf)
            else
                pdf = 0._RKC
            end if
        else
            call setBetaLogPDF(pdf, x, alpha, beta)
            pdf = exp(pdf)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getBetaLogPDF_ENABLED || setBetaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0._RKC < x .and. x < 1._RKC, SK_"@setBetaLogPDF(): The condition `0. < x .and. x < 1.` must hold. x = "//getStr(x))
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@setBetaLogPDF(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@setBetaLogPDF(): The condition `0. < beta` must hold. alpha = "//getStr(beta))
#if     LOGBETA_ENABLED
        CHECK_ASSERTION(__LINE__, abs(logBeta - getLogBeta(alpha, beta)) < sqrt(epsilon(0._RKC)), \
        SK_"@setBetaLogPDF(): The condition `abs(logBeta - getLogBeta(alpha, beta)) < sqrt(epsilon(0._RKC)` must hold. logBeta, getLogBeta(alpha, beta) = "\
        //getStr([logBeta, getLogBeta(alpha, beta)]))
#endif
        logPDF = (alpha - 1._RKC) * log(x) + (beta - 1._RKC) * log(1._RKC - x) - & ! LCOV_EXCL_LINE
#if     DEFAULT_ENABLED
        getLogBeta(alpha, beta)
#elif   LOGBETA_ENABLED
        logBeta
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%
#elif   getBetaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%
        
        cdf = getBetaInc(x, alpha, beta, signed)
        
        !%%%%%%%%%%%%%%%%%
#elif   setBetaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%
        
        call setBetaInc(cdf, x, alpha, beta, logFuncBeta, signed, info)

        !%%%%%%%%%%%%%%%%%%
#elif   setBetaRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        ! Set the URNG.
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0._RKC < alpha, SK_"@setBetaRand(): The condition `0. < alpha` must hold. alpha = "//getStr(alpha))
        CHECK_ASSERTION(__LINE__, 0._RKC < beta, SK_"@setBetaRand(): The condition `0. < beta` must hold. alpha = "//getStr(beta))
        if (1._RKC < alpha .or. 1._RKC < beta) then ! Use the algorithm of Johnk.
            block
#if             D0_ENABLED
                real(RKC) :: temp
#elif           D1_ENABLED
                real(RKC) :: temp(size(rand, 1, IK))
#else       
#error          "Unrecognized interface."
#endif
                call setGammaRand(RNG rand, alpha, 1._RKC)
                call setGammaRand(RNG temp, beta, 1._RKC)
                rand = rand / (rand + temp)
            end block
        else
            block
#if             D1_ENABLED
                integer(IK) :: irand
#endif
                real(RKC)   :: temp, invShape(2), dumm
                invShape(1) = 1._RKC / alpha
                invShape(2) = 1._RKC / beta
#if             D1_ENABLED
#define         GET_RAND(i) rand(i)
                do irand = 1, size(rand, 1, IK)
#elif               D0_ENABLED
#define             GET_RAND(i) rand
#endif
                    do
                        call setUnifRand(RNG GET_RAND(irand))
                        call setUnifRand(RNG temp)
                        GET_RAND(irand) = (1._RKC - GET_RAND(irand))**invShape(1)
                        temp = (1._RKC - temp)**invShape(2)
                        dumm = GET_RAND(irand) + temp
                        ! Rejection and cycling happens only if any(unifrnd == 0._RKC),
                        ! which is approximately 1 in 10^106 odds of occurrence.
                        if (1._RKC < dumm) cycle
                        if (0._RKC < dumm) then
                            GET_RAND(irand) = GET_RAND(irand) / dumm
                        else
                            GET_RAND(irand) = log(GET_RAND(irand))
                            temp = log(temp)
                            dumm = max(GET_RAND(irand), temp)
                            GET_RAND(irand) = GET_RAND(irand) - dumm
                            temp = temp - dumm
                            GET_RAND(irand) = exp(GET_RAND(irand) - log(exp(GET_RAND(irand)) + exp(temp)))
                        end if
                        exit
                    end do
#if             D1_ENABLED
                end do
#endif
            end block
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  GET_RAND
#undef  RNG
