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
!>  This include file contains the implementation of procedures in [pm_distGamma](@ref pm_distGamma).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the URNG.
#if     RNGD_ENABLED
#define RNG
#elif   RNGF_ENABLED || RNGX_ENABLED
#define RNG rng,
#elif   setGammaRand_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%
#if     getGammaLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@getGammaLogPDFNF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa))
#if     KD_ENABLED
        logPDFNF = -log_gamma(kappa)
#elif   KS_ENABLED
        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@getGammaLogPDFNF(): The condition invSigma > 0.` must hold. invSigma = "//getStr(invSigma))
        logPDFNF = getGammaLogPDFNF(kappa)
        if (invSigma /= 1._RKG) logPDFNF = logPDFNF + log(invSigma)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%
#elif   getGammaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

        real(RKG) :: kappa_def, invSigma_def
        kappa_def = 1._RKG; if (present(kappa)) kappa_def = kappa
        invSigma_def = 1._RKG; if (present(invSigma)) invSigma_def = invSigma
        call setGammaLogPDF(logPDF, x, getGammaLogPDFNF(kappa_def, invSigma_def), kappa_def, invSigma_def)

        !%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%

#if     DDD_ENABLED
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@setGammaLogPDF(): The condition `x > 0.` must hold. x = "//getStr(x))
        logPDF = -x
#elif   NKD_ENABLED
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@setGammaLogPDF(): The condition `x > 0.` must hold. x = "//getStr(x))
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@setGammaLogPDF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa))
        CHECK_ASSERTION(__LINE__, abs(logPDFNF - getGammaLogPDFNF(kappa)) < sqrt(epsilon(0._RKG)), \
        SK_"@setGammaLogPDF(): The condition `abs(logPDFNF - getGammaLogPDFNF(kappa)) < sqrt(epsilon(0._RKG)` must hold. logPDFNF, getGammaLogPDFNF(kappa) = "// \
        getStr([logPDFNF, getGammaLogPDFNF(kappa)]))
        logPDF = logPDFNF + (kappa - 1._RKG) * log(x) - x
#elif   NKS_ENABLED
        real(RKG) :: y
        y = x * invSigma
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@setGammaLogPDF(): The condition `x > 0.` must hold. x = "//getStr(x))
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@setGammaLogPDF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa))
        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@setGammaLogPDF(): The condition `invSigma > 0.` must hold. invSigma = "//getStr(invSigma))
        CHECK_ASSERTION(__LINE__, abs(logPDFNF - getGammaLogPDFNF(kappa, invSigma)) < sqrt(epsilon(0._RKG)), \
        SK_"@setGammaLogPDF(): The condition `abs(logPDFNF - getGammaLogPDFNF(kappa, invSigma)) < sqrt(epsilon(0._RKG)` must hold. logPDFNF, getGammaLogPDFNF(kappa, invSigma) = "// \
        getStr([logPDFNF, getGammaLogPDFNF(kappa, invSigma)]))
        logPDF = logPDFNF + (kappa - 1._RKG) * log(y) - y
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%
#elif   getGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
        real(RKG)   :: xnormed
        if (present(invSigma)) then
            xnormed = x * invSigma
        else
            xnormed = x
        end if
        if (present(kappa)) then
            call setGammaCDF(cdf, xnormed, kappa, info)
        else
            call setGammaCDF(cdf, xnormed, info)
        end if
        if (info < 0_IK) error stop MODULE_NAME//SK_"@getGammaCDF(): The computation of the regularized Lower Incomplete Gamma function failed. This can happen if `kappa` is too large."

        !%%%%%%%%%%%%%%%%%%
#elif   setGammaCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%

        real(RKG) :: gamincupp
#if     DD_ENABLED
        real(RKG), parameter :: kappa = 1._RKG!, logGammaKappa = log_gamma(kappa)
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@setGammaCDF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        !call setGammaIncLow(cdf, x, logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, x, kappa, info)
#else
        CHECK_ASSERTION(__LINE__, x > 0._RKG, SK_"@setGammaCDF(): The condition `x > 0.` must hold. x = "//getStr(x)) ! fpp
        CHECK_ASSERTION(__LINE__, kappa > 0._RKG, SK_"@setGammaCDF(): The condition `kappa > 0.` must hold. kappa = "//getStr(kappa)) ! fpp
        !check_assertion(__LINE__, abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG), \
        !SK_"@setGammaCDF(): The condition `abs(log_gamma(kappa) - logGammaKappa) < 100 * epsilon(0._RKG)` must hold. log_gamma(kappa), logGammaKappa = "//\
        !getStr([log_gamma(kappa), logGammaKappa])) ! fpp
#if     KD_ENABLED
        !call setGammaIncLow(cdf, x, logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, x, kappa, info)
#elif   KS_ENABLED
        CHECK_ASSERTION(__LINE__, invSigma > 0._RKG, SK_"@setGammaCDF(): The condition `invSigma > 0.` must hold. invSigma = "//getStr(invSigma)) ! fpp
        !call setGammaIncLow(cdf, x * invSigma, logGammaKappa, kappa, info)
        call setGammaInc(cdf, gamincupp, x * invSigma, kappa, info)
#else
#error  "Unrecognized interface."
#endif
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaRand_ENABLED && (RNGD_ENABLED || RNGF_ENABLED || RNGX_ENABLED) && KR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG)   , parameter :: ONE_THIRD = 1._RKG / 3._RKG
        real(RKG)               :: invSqrt9KappaMinusOneThird
        real(RKG)               :: kappaMinusOneThird
        real(RKG)               :: normrnd, unifrnd
        CHECK_ASSERTION(__LINE__, 0._RKG < kappa .and. 0._RKG < sigma, SK_"@setGammaRand(): The condition `0 < kappa .and. 0 < sigma` must hold. kappa = "//getStr([kappa, sigma])) ! Must appear here.
#define GET_PARAM(offset) \
kappaMinusOneThird = kappa - ONE_THIRD + offset; invSqrt9KappaMinusOneThird = 1._RKG / sqrt(9 * kappaMinusOneThird);
        if (1._RKG <= kappa) then
            GET_PARAM(0._RKG) ! fpp
#if         D1_ENABLED
#define     GET_RAND(i) rand(i)
            block
                integer(IK) :: irand
                do irand = 1_IK, size(rand, 1, IK)
#elif               D0_ENABLED
#define             GET_RAND(i) rand
#else
#error              "Unrecognized interface."
#endif
                    do
                        !call setRandKappaGe1()
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ! If only macros existed in Fortran...
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        call setNormRand(RNG normrnd)
                        GET_RAND(irand) = 1._RKG + invSqrt9KappaMinusOneThird * normrnd
                        if (GET_RAND(irand) <= 0._RKG) cycle
                        call setUnifRand(RNG unifrnd)
                        unifrnd = 1._RKG - unifrnd
                        if (unifrnd < 1._RKG - 0.0331_RKG * normrnd**4) then
                            GET_RAND(irand) = kappaMinusOneThird * sigma * GET_RAND(irand)**3
                            exit
                        end if
                        GET_RAND(irand) = GET_RAND(irand)**3
                        if (log(unifrnd) < 0.5_RKG * normrnd**2 + kappaMinusOneThird * (1._RKG - GET_RAND(irand) + log(GET_RAND(irand)))) then
                            GET_RAND(irand) = kappaMinusOneThird * sigma * GET_RAND(irand)
                            exit
                        end if
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end do
#if             D1_ENABLED
                end do
            end block
#endif
        else
#if         1
            block
                real(RKG) :: kapinv
#if             D1_ENABLED
                real(RKG) :: unifrnd(size(rand, 1, IK))
#endif
                kapinv = 1._RKG / kappa
                call setGammaRand(RNG rand, kappa = kappa + 1._RKG, sigma = sigma)
                call setUnifRand(RNG unifrnd)
                rand = rand * (1._RKG - unifrnd)**kapinv
            end block
#else
            GET_PARAM(1._RKG) ! fpp
#if         D1_ENABLED
            block
                integer(IK) :: irand
                do irand = 1_IK, size(rand, 1, IK)
#endif
                    ! This approach seems to be buggy and yield incorrect results.
                    do
                        !call setRandKappaGe1()
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ! If only macros existed in Fortran...
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        call setNormRand(RNG normrnd)
                        GET_RAND(irand) = 1._RKG + invSqrt9KappaMinusOneThird * normrnd
                        if (GET_RAND(irand) <= 0._RKG) cycle
                        call setUnifRand(RNG unifrnd); unifrnd = 1._RKG - unifrnd
                        if (unifrnd < 1._RKG - 0.0331_RKG * normrnd**4) then
                            GET_RAND(irand) = kappaMinusOneThird * sigma * GET_RAND(irand)**3
                            exit
                        end if
                        GET_RAND(irand) = GET_RAND(irand)**3
                        if (log(unifrnd) < 0.5_RKG * normrnd**2 + kappaMinusOneThird * (1._RKG - GET_RAND(irand) + log(GET_RAND(irand)))) then
                            GET_RAND(irand) = kappaMinusOneThird * sigma * GET_RAND(irand)
                            exit
                        end if
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        call setUnifRand(RNG unifrnd)
                        unifrnd = (1._RKG - unifrnd)**(1._RKG / kappa)
                        GET_RAND(irand) = GET_RAND(irand) * unifrnd
                    end do
#if             D1_ENABLED
                end do
            end block
#endif
#endif
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setGammaRand_ENABLED && (RNGD_ENABLED || RNGF_ENABLED || RNGX_ENABLED) && KI_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Alternative algorithm when `kappa` is a positive integer.
        ! The corresponding interfaces are yet to be implemented.
        ! However, the benefit of this integer `kappa` is unclear.
        real(RKG) :: unifrnd(7)
        CHECK_ASSERTION(__LINE__, 0._RKG < kappa .and. 0._RKG < sigma, SK_"@setGammaRand(): The condition `0 < kappa .and. 0 < sigma` must hold. kappa = "//getStr([kappa, sigma]))
        if (kappa < 6_IK) then
            call setUnifRand(RNG unifrnd(1:kappa))
            rand = -log(product(unifrnd(1:kappa)))
        else ! use rejection sampling
            do
                call setUnifRand(RNG unifrnd(1:2))
                unifrnd(1:2) = 2 * unifrnd(1:2) - 1._RKG
                if (dot_product(unifrnd(1:2), unifrnd(1:2)) > 1._RKG) cycle
                unifrnd(3) = unifrnd(2) / unifrnd(1)
                unifrnd(4) = kappa - 1._RKG
                unifrnd(5) = sqrt(2 * unifrnd(4) + 1._RKG)
                rand = unifrnd(5) * unifrnd(3) + unifrnd(4)
                if (rand <= 0.0) cycle
                unifrnd(6) = (1._RKG + unifrnd(3)**2) * exp(unifrnd(4) * log(rand / unifrnd(4)) - unifrnd(5) * unifrnd(3))
                call setUnifRand(RNG unifrnd(7)) !call random number(unifrnd(7))
                if (unifrnd(7) <= unifrnd(6)) exit
            end do
        end if
#if     0
        ! Alternative old implementation.
        function getRandGamma(alpha) result(randGamma)
            use pm_distNorm, only: setNormRand
            implicit none
            real(RK), intent(in) :: alpha
            real(RK)             :: randGamma
            real(RK)             :: c,u,v,z
            if (alpha<=0._RK) then  ! illegal value of alpha
                randGamma = -1._RK
                return
            else
                randGamma = alpha
                if (randGamma < 1._RK) randGamma = randGamma + 1._RK
                randGamma = randGamma - 0.3333333333333333_RK
                c = 3._RK*sqrt(randGamma)
                c = 1._RK / c
                do
                    do
                        call setNormRand(z)
                        v = 1._RK + c*z
                        if (v<=0._RK) cycle
                        exit
                    end do
                    v = v**3
                    call setUnifRand(RNG u)
                    if (log(u) >= 0.5_RK * z**2 + randGamma * (1._RK - v + log(v))) cycle
                    randGamma = randGamma * v
                    exit
                end do
                if (alpha < 1._RK) then
                    call setUnifRand(RNG u)
                    randGamma = randGamma * u**(1._RK/alpha)
                end if
            end if
        end function getRandGamma
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  GET_PARAM
#undef  GET_RAND
#undef  RNG