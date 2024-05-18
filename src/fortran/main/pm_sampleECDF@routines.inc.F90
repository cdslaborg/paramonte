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
!>  This file contains the implementation details of the routines under the generic interface [pm_sampleECDF](@ref pm_sampleECDF).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! The accumulation method accumulates significant error for large arrays.
! Either multiplication must be used, or the accurate summation method of fabcomp().
! The current accumulation workaround is ugly and needs a revamp using fabcomp().
#define ACCUMULATION_ENABLED 1
#define MULTIPLICATION_ENABLED 0

        !%%%%%%%%%%%%%%
#if     setECDF_ENABLED
        !%%%%%%%%%%%%%%

#define CHECK_LEN_LCDF \
CHECK_ASSERTION(__LINE__, size(lcdf, 1, IK) == size(ecdf, 1, IK), \
SK_"@setECDF(): The condition `size(lcdf) == size(ecdf)` must hold. size(lcdf), size(ecdf) = "//getStr([size(lcdf, 1, IK), size(ecdf, 1, IK)]))
#define CHECK_LEN_UCDF \
CHECK_ASSERTION(__LINE__, size(ucdf, 1, IK) == size(ecdf, 1, IK), \
SK_"@setECDF(): The condition `size(ucdf) == size(ecdf)` must hold. size(ucdf), size(ecdf) = "//getStr([size(ucdf, 1, IK), size(ecdf, 1, IK)]))

        integer(IK) :: nsam, isam
        real(TKC), parameter :: ZERO = 0._TKC, UNIT = 1._TKC
        real(TKC):: nsamInv
        nsam = size(ecdf, 1, IK)
        if (nsam == 0_IK) return
#if     ONE_ENABLED
#define GET_WEIGHTED(W,X)X
        nsamInv = UNIT / real(nsam, TKC)
#elif   WIK_ENABLED || WRK_ENABLED
#define GET_WEIGHTED(W,X)real(W, TKC) * X
        CHECK_ASSERTION(__LINE__, all(0 <= weight), SK_"@setECDF(): The condition `all(0 <= weight)` must hold. weight = "//getStr(weight))
        CHECK_ASSERTION(__LINE__, abs(weisum - sum(weight)) < epsilon(0._TKC) * 100, SK_"@setECDF(): The condition `abs(weisum - sum(weight)) < epsilon(0._TKC) * 100` must hold. weisum, sum(weight) = "//getStr([weisum, sum(weight)]))
        nsamInv = UNIT / real(weisum, TKC)
#else
#error  "Unrecognized interface."
#endif
        if (1_IK < nsam) then
            ecdf(1) = GET_WEIGHTED(weight(1),nsamInv)
#if         !WRK_ENABLED
            blockExpectedLowerUpper: if (.not. (present(lcdf) .or. present(ucdf))) then
#endif
#if             ACCUMULATION_ENABLED || WIK_ENABLED || WRK_ENABLED
#define         ISAM isam
                do isam = 2, nsam - 1
                    ecdf(isam) = ecdf(isam - 1) + GET_WEIGHTED(weight(isam),nsamInv)
                end do
                if (1._TKC < ecdf(nsam - 1)) then
                    nsamInv = 1._TKC / (ecdf(nsam - 1) + GET_WEIGHTED(weight(nsam),nsamInv))
                    do concurrent(isam = 1 : nsam - 1)
                        ecdf(isam) = ecdf(isam) * nsamInv
                    end do
                end if
#elif           MULTIPLICATION_ENABLED
#define         ISAM nsam
                do concurrent(isam = 2 : nsam - 1)
                    ecdf(isam) = real(isam, TKC) * nsamInv
                end do
#else
#error          "Unrecognized interface."
#endif
                ecdf(ISAM) = UNIT
#if         !WRK_ENABLED
            else blockExpectedLowerUpper
                block
                    real(TKC):: sqrt_nsamInvTwice_logTwoOverAlpha
                    if (present(alpha)) then
                        CHECK_ASSERTION(__LINE__, ZERO < alpha .and. alpha < UNIT, SK_"@setECDF(): The condition `0 < alpha .and. alpha < 1` must hold. alpha = "//getStr(alpha))
                        sqrt_nsamInvTwice_logTwoOverAlpha = sqrt(0.5_TKC * nsamInv * log(2._TKC / alpha))
                    else
                        sqrt_nsamInvTwice_logTwoOverAlpha = sqrt(0.5_TKC * nsamInv * log(2._TKC / 0.05_TKC))
                    end if
                    if (present(lcdf) .and. present(ucdf)) then
                        CHECK_LEN_LCDF ! fpp
                        CHECK_LEN_UCDF ! fpp
                        lcdf(1) = max(ZERO, GET_WEIGHTED(weight(1),nsamInv) - sqrt_nsamInvTwice_logTwoOverAlpha)
                        ucdf(1) = min(UNIT, GET_WEIGHTED(weight(1),nsamInv) + sqrt_nsamInvTwice_logTwoOverAlpha)
#if                     ACCUMULATION_ENABLED || WIK_ENABLED || WRK_ENABLED
                        do isam = 2, nsam - 1
                            ecdf(isam) = ecdf(isam - 1) + GET_WEIGHTED(weight(isam),nsamInv)
                            lcdf(isam) = max(ZERO, ecdf(isam) - sqrt_nsamInvTwice_logTwoOverAlpha)
                            ucdf(isam) = min(UNIT, ecdf(isam) + sqrt_nsamInvTwice_logTwoOverAlpha)
                        end do
#elif                   MULTIPLICATION_ENABLED
                        do concurrent(isam = 2 : nsam - 1)
                            ecdf(isam) = real(isam, TKC) * nsamInv
                            lcdf(isam) = max(ZERO, ecdf(isam) - sqrt_nsamInvTwice_logTwoOverAlpha)
                            ucdf(isam) = min(UNIT, ecdf(isam) + sqrt_nsamInvTwice_logTwoOverAlpha)
                        end do
#endif
                        ecdf(ISAM) = UNIT
                        lcdf(ISAM) = max(ZERO, UNIT - sqrt_nsamInvTwice_logTwoOverAlpha)
                        ucdf(ISAM) = UNIT
                    elseif (present(lcdf)) then
                        CHECK_LEN_LCDF ! fpp
                        lcdf(1) = max(ZERO, GET_WEIGHTED(weight(1),nsamInv) - sqrt_nsamInvTwice_logTwoOverAlpha)
#if                     ACCUMULATION_ENABLED || WIK_ENABLED || WRK_ENABLED
                        do isam = 2, nsam - 1
                            ecdf(isam) = ecdf(isam - 1) + GET_WEIGHTED(weight(isam),nsamInv)
                            lcdf(isam) = max(ZERO, ecdf(isam) - sqrt_nsamInvTwice_logTwoOverAlpha)
                        end do
#elif                   MULTIPLICATION_ENABLED
                        do concurrent(isam = 2 : nsam - 1)
                            ecdf(isam) = real(isam, TKC) * nsamInv
                            lcdf(isam) = max(ZERO, ecdf(isam) - sqrt_nsamInvTwice_logTwoOverAlpha)
                        end do
#endif
                        ecdf(ISAM) = UNIT
                        lcdf(ISAM) = max(ZERO, UNIT - sqrt_nsamInvTwice_logTwoOverAlpha)
                    elseif (present(ucdf)) then
                        CHECK_LEN_UCDF ! fpp
                        ucdf(1) = min(UNIT, GET_WEIGHTED(weight(1),nsamInv) + sqrt_nsamInvTwice_logTwoOverAlpha)
#if                     ACCUMULATION_ENABLED || WIK_ENABLED || WRK_ENABLED
                        do isam = 2, nsam - 1
                            ecdf(isam) = ecdf(isam - 1) + GET_WEIGHTED(weight(isam),nsamInv)
                            ucdf(isam) = min(UNIT, ecdf(isam) + sqrt_nsamInvTwice_logTwoOverAlpha)
                        end do
#elif                   MULTIPLICATION_ENABLED
                        do concurrent(isam = 2 : nsam - 1)
                            ecdf(isam) = real(isam, TKC) * nsamInv
                            ucdf(isam) = min(UNIT, ecdf(isam) + sqrt_nsamInvTwice_logTwoOverAlpha)
                        end do
#endif
                        ecdf(ISAM) = UNIT
                        ucdf(ISAM) = UNIT
                    end if
                end block
            end if blockExpectedLowerUpper
#endif
        else
            ! nsam == 1
            ecdf = UNIT
#if         !WRK_ENABLED
            if (present(lcdf)) lcdf = UNIT
            if (present(ucdf)) ucdf = UNIT
#endif
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_LEN_LCDF
#undef  CHECK_LEN_UCDF
#undef  GET_WEIGHTED
#undef  ISAM