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
!>  This include file contains the implementations of procedures of [pm_distBand](@ref pm_distBand).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CHECK_ALPHA_NEQ_TWO \
CHECK_ASSERTION(__LINE__, alpha /= -2._RKG, SK_": The condition `alpha /= -2.` must hold. alpha = "//getStr(alpha))
#define CHECK_ALPHA_GEQ_BETA \
CHECK_ASSERTION(__LINE__, beta < alpha, SK_": The condition `beta < alpha` must hold. beta, alpha = "//getStr([beta, alpha]))
#define CHECK_EBREAK_GEQ_ZERO \
CHECK_ASSERTION(__LINE__, 0._RKG < ebreak, SK_": The condition `0. < ebreak` must hold. ebreak = "//getStr(ebreak))

        !%%%%%%%%%%%%%%%%%%%
#if     getBandEpeak_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        CHECK_ALPHA_NEQ_TWO
        CHECK_ALPHA_GEQ_BETA
        CHECK_EBREAK_GEQ_ZERO
        epeak = ebreak * (alpha + 2._RKG) / (alpha - beta)

        !%%%%%%%%%%%%%%%%%%%%
#elif   getBandEbreak_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        CHECK_ALPHA_NEQ_TWO
        CHECK_ALPHA_GEQ_BETA
        CHECK_ASSERTION(__LINE__, 0._RKG < epeak, SK_"@getBandEbreak(): The condition `0. < epeak` must hold. epeak = "//getStr(epeak))
        ebreak = epeak * (alpha - beta) / (alpha + 2._RKG)

        !%%%%%%%%%%%%%%%%%%
#elif   getBandZeta_ENABLED
        !%%%%%%%%%%%%%%%%%%

        CHECK_EBREAK_GEQ_ZERO
        zeta = ebreak**(alpha - beta) * exp(beta - alpha)
        !logZeta = (alpha - beta) * (ebreak - 1._RKG)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getBandUDF_ENABLED && Any_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ALPHA_NEQ_TWO
        CHECK_ALPHA_GEQ_BETA
        CHECK_EBREAK_GEQ_ZERO
        CHECK_ASSERTION(__LINE__, 0._RKG < energy, SK_"@getBandUDF(): The condition `0. < energy` must hold. energy = "//getStr(energy))
        if (energy < ebreak) then
            if (present(invEfold)) then
                CHECK_ASSERTION(__LINE__, abs(invEfold - (alpha - beta) / ebreak) < 10 * epsilon(energy), \
                SK_"@getBandUDF(): The condition `abs(invEfold - (alpha - beta) / ebreak) < 10 * epsilon(energy)` must hold. alpha, beta, ebreak, invEfold, (alpha - beta) / ebreak, epsilon(energy) = "//\
                getStr([alpha, beta, ebreak, invEfold, (alpha - beta) / ebreak, epsilon(energy)]))
                udf = energy**alpha * exp(-energy * invEfold)
            else
                udf = energy**alpha * exp((beta - alpha) * (energy / ebreak))
            end if
        elseif (present(zeta)) then
            CHECK_ASSERTION(__LINE__, abs(zeta - getBandZeta(alpha, beta, ebreak)) < 10 * epsilon(energy), \
            SK_"@getBandUDF(): The condition `abs(zeta - getBandZeta(alpha, beta, ebreak)) < 10 * epsilon(energy)` must hold. alpha, beta, ebreak, zeta, getBandZeta(alpha, beta, ebreak), epsilon(energy) = "//\
            getStr([alpha, beta, ebreak, zeta, getBandZeta(alpha, beta, ebreak), epsilon(energy)]))
            udf = zeta * energy**beta
        else
            udf = getBandZeta(alpha, beta, ebreak) * energy**beta
        end if

        !%%%%%%%%%%%%%%%%%%
#elif   setBandUCDF_ENABLED
        !%%%%%%%%%%%%%%%%%%

        real(RKG) :: mb, tmp, invEfold

        CHECK_ALPHA_NEQ_TWO
        CHECK_ALPHA_GEQ_BETA
        CHECK_EBREAK_GEQ_ZERO

        info = 0_IK
        ucdf = 0._RKG
        if (lb < ub) then ! integrate the spectrum
            CHECK_ASSERTION(__LINE__, 0._RKG < lb, SK_"@setBandUCDF(): The condition `0. < lb` must hold. lb = "//getStr(lb))
            if (lb < ebreak) then
                ! First compute the contribution from the lower component via either Gamma CDF or brute-force integration.
                mb = min(ub, ebreak)
                invEfold = (alpha - beta) / ebreak
                if (-1._RKG < alpha) then
                    ! use Gamma CDF.
                    block
                        real(RKG) :: kappa!, dumm
                        kappa = alpha + 1._RKG
                        !dumm = log_gamma(kappa)
                        !call setGammaCDF(ucdf, lb, dumm, kappa, invEfold, info)
                        call setGammaCDF(ucdf, lb, kappa, invEfold, info)
                        if (info < 0_IK) return ! LCOV_EXCL_LINE
                        !call setGammaCDF(tmp, mb, dumm, kappa, invEfold, info)
                        call setGammaCDF(tmp, mb, kappa, invEfold, info)
                        if (info < 0_IK) return ! LCOV_EXCL_LINE
                        ucdf = (tmp - ucdf) * gamma(kappa) / invEfold**kappa
                    end block
                else
                    block
                        integer(IK) :: neval, nint, sindex(1000)
                        real(RKG) :: abstol, reltol, abserr, sinfo(4, 1000)
                        abstol = 0._RKG
                        reltol = epsilon(0._RKG)**(2./3.)
                        info = getQuadErr(getBandLowUDF, lb, mb, abstol, reltol, GK21, weps, ucdf, abserr, sinfo, sindex, neval, nint)
                        if (info /= 0_IK) info = -info ! LCOV_EXCL_LINE
                        if (info < 0_IK) return ! LCOV_EXCL_LINE
                    end block
                end if
                if (mb == ub) return
            else
                mb = lb
            end if
            ! Add the remaining part of the ucdf from the high-energy component.
            tmp = beta + 1._RKG
            if (tmp /= 0._RKG) then
                ucdf = ucdf + getBandZeta(alpha, beta, ebreak) * (ub**tmp - mb**tmp) / tmp
            else
                ucdf = ucdf + getBandZeta(alpha, beta, ebreak) * log(ub / mb)
            end if
        end if

    contains

        PURE function getBandLowUDF(energy) result(udf)
            implicit none
            real(RKG), intent(in)    :: energy
            real(RKG)                :: udf
            CHECK_ASSERTION(__LINE__, lb <= energy .and. energy <= ub, SK_"@setBandUCDF(): The condition `lb <= energy .and. energy <= ub` must hold. lb, energy, ub = "//getStr([lb, energy, ub]))
            udf = energy**alpha * exp(-invEfold * energy)
            !print *, udf, invEfold, energy
        end function

        !%%%%%%%%%%%%%%%%%%
#elif   setBandMean_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     Def_ENABLED
#define LBNEW lb
#define UBNEW ub
#elif   !New_ENABLED
#error  "Unrecognized interface."
#endif
        real(RKG) :: denom
        call setBandUCDF(denom, lb, ub, alpha, beta, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        call setBandUCDF(mean, LBNEW, UBNEW, alpha + 1._RKG, beta + 1._RKG, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        mean = mean / denom
#undef  LBNEW
#undef  UBNEW

        !%%%%%%%%%%%%%%%%%%%%
#elif   setBandPhoton_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

#if     FromPhoton_ENABLED && NewB_ENABLED
        real(RKG) :: denom, numer
        CHECK_ASSERTION(__LINE__, lb <= ub, SK_"@setBandPhoton(): The condition `lb <= ub` must hold. lb, ub = "//getStr([lb, ub]))
        CHECK_ASSERTION(__LINE__, lbnew <= ubnew, SK_"@setBandPhoton(): The condition `lbnew <= ubnew` must hold. lbnew, ubnew = "//getStr([lbnew, ubnew]))
        CHECK_ASSERTION(__LINE__, 0._RKG < photon, SK_"@setBandPhoton(): The condition `0. < photon` must hold. photon = "//getStr(photon))
        call setBandUCDF(denom, lb, ub, alpha, beta, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        call setBandUCDF(numer, lbnew, ubnew, alpha, beta, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        photon = numer * (photon / denom)
#elif   FromEnergy_ENABLED
#if     OldB_ENABLED
#define LBNEW lb
#define UBNEW ub
#elif   !NewB_ENABLED
#error  "Unrecognized interface."
#endif
        real(RKG) :: denom, numer
        CHECK_ASSERTION(__LINE__, 0._RKG < energy, SK_"@setBandPhoton(): The condition `0. < energy` must hold. energy = "//getStr(energy))
        call setBandUCDF(denom, lb, ub, alpha + 1._RKG, beta + 1._RKG, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        call setBandUCDF(numer, LBNEW, UBNEW, alpha, beta, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        photon = numer * (energy / denom)
#undef  LBNEW
#undef  UBNEW
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%
#elif   setBandEnergy_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

#if     FromEnergy_ENABLED && NewB_ENABLED
        real(RKG) :: denom, numer
        CHECK_ASSERTION(__LINE__, lb <= ub, SK_"@setBandPhoton(): The condition `lb <= ub` must hold. lb, ub = "//getStr([lb, ub]))
        CHECK_ASSERTION(__LINE__, lbnew <= ubnew, SK_"@setBandPhoton(): The condition `lbnew <= ubnew` must hold. lbnew, ubnew = "//getStr([lbnew, ubnew]))
        CHECK_ASSERTION(__LINE__, 0._RKG < energy, SK_"@setBandPhoton(): The condition `0. < energy` must hold. energy = "//getStr(energy))
        call setBandUCDF(denom, lb, ub, alpha + 1._RKG, beta + 1._RKG, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        call setBandUCDF(numer, lbnew, ubnew, alpha + 1._RKG, beta + 1._RKG, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        energy = numer * (energy / denom)
#elif   FromPhoton_ENABLED
#if     OldB_ENABLED
#define LBNEW lb
#define UBNEW ub
#elif   !NewB_ENABLED
#error  "Unrecognized interface."
#endif
        real(RKG) :: denom, numer
        CHECK_ASSERTION(__LINE__, 0._RKG < photon, SK_"@setBandPhoton(): The condition `0. < photon` must hold. photon = "//getStr(photon))
        call setBandUCDF(denom, lb, ub, alpha, beta, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        call setBandUCDF(numer, LBNEW, UBNEW, alpha + 1._RKG, beta + 1._RKG, ebreak, info)
        if (info < 0_IK) return ! LCOV_EXCL_LINE
        energy = numer * (photon / denom)
#undef  LBNEW
#undef  UBNEW
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_EBREAK_GEQ_ZERO
#undef  CHECK_ALPHA_GEQ_BETA
#undef  CHECK_ALPHA_NEQ_TWO