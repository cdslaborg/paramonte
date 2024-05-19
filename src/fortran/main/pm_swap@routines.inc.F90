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
!>  This include file contains procedure implementations of [pm_swap](@ref pm_swap).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 AM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Check for string length consistency.
#if     SK_ENABLED
#define CHECK_STRLEN \
CHECK_ASSERTION(__LINE__, len(a, IK) == len(b, IK), SK_"@setSwapped(): The condition `len(a) == len(b)` must hold. len(a), len(b) = "//getStr([len(a, IK), len(b, IK)])) ! fpp
#else
#define CHECK_STRLEN
#endif
        ! Define the sizing and slicing rules.
#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE len
#define GET_SLICE(i) i:i
#elif   D1_ENABLED
#define GET_SIZE size
#define GET_SLICE(i) i
#elif   !D0_ENABLED
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%
#if     setSwapped_ENABLED
        !%%%%%%%%%%%%%%%%%

        ! Define the place holder.
#if     !(BLAS_ENABLED && DISPATCH_ENABLED)
#if     SK_ENABLED && D0_ENABLED
        character(1,SKG) :: tmp
#elif   SK_ENABLED && D1_ENABLED
        character(len(a,IK),SKG) :: tmp
#elif   IK_ENABLED
        integer(IKG)    :: tmp
#elif   LK_ENABLED
        logical(LKG)    :: tmp
#elif   CK_ENABLED
        complex(CKG)    :: tmp
#elif   RK_ENABLED
        real(RKG)       :: tmp
#else
#error  "Unrecognized interface."
#endif
#endif
#if     D0_ENABLED && !SK_ENABLED
        CHECK_STRLEN
        tmp = a
        a = b
        b = tmp
#elif   Def_ENABLED && D1_ENABLED && BLAS_ENABLED && DISPATCH_ENABLED
        CHECK_STRLEN
        CHECK_ASSERTION(__LINE__, size(a, 1, IK) == size(b, 1, IK), SK_"@setSwapped(): The condition `size(a) == size(b)` must hold. size(a), size(b) = "//getStr([size(a), size(b)])) ! fpp
        call blasSWAP(size(a, 1, IK), a, 1_IK, b, 1_IK)
#elif   Def_ENABLED && (D1_ENABLED || (D0_ENABLED && SK_ENABLED))
        integer(IK) :: iell, nell, mell
        nell = GET_SIZE(a, kind = IK)
        CHECK_STRLEN
        CHECK_ASSERTION(__LINE__, nell == GET_SIZE(b, kind = IK), SK_"@setSwapped(): The condition `size(a) == size(b)` must hold. size(a), size(b) = "//getStr([GET_SIZE(a, kind = IK), GET_SIZE(b, kind = IK)])) ! fpp
        mell = mod(nell, 3_IK)
        if (mell /= 0_IK) then
            do iell = 1, mell
                tmp = a(GET_SLICE(iell))
                a(GET_SLICE(iell)) = b(GET_SLICE(iell))
                b(GET_SLICE(iell)) = tmp
            end do
            if (nell < 3_IK) return
        end if
        do iell = mell + 1, nell, 3
            tmp = a(GET_SLICE(iell))
            a(GET_SLICE(iell)) = b(GET_SLICE(iell))
            b(GET_SLICE(iell)) = tmp
            tmp = a(GET_SLICE(iell + 1))
            a(GET_SLICE(iell + 1)) = b(GET_SLICE(iell + 1))
            b(GET_SLICE(iell + 1)) = tmp
            tmp = a(GET_SLICE(iell + 2))
            a(GET_SLICE(iell + 2)) = b(GET_SLICE(iell + 2))
            b(GET_SLICE(iell + 2)) = tmp
        end do
#elif   Inc_ENABLED && (D1_ENABLED || (D0_ENABLED && SK_ENABLED))
        integer(IK) :: nell
        if (inca == 1_IK .and. incb == 1_IK) then
            call setSwapped(a, b)
        else
            CHECK_ASSERTION(__LINE__, inca /= 0_IK .or. GET_SIZE(a, kind = IK) == 1_IK, SK_"@setSwapped(): The condition `inca /= 0 .or. size(a) == 1` must hold. size(a), inca = "//getStr([GET_SIZE(a, kind = IK), inca])) ! fpp
            CHECK_ASSERTION(__LINE__, incb /= 0_IK .or. GET_SIZE(b, kind = IK) == 1_IK, SK_"@setSwapped(): The condition `incb /= 0 .or. size(b) == 1` must hold. size(b), incb = "//getStr([GET_SIZE(b, kind = IK), incb])) ! fpp
            CHECK_ASSERTION(__LINE__, (GET_SIZE(a, kind = IK) - 1) / max(1_IK,abs(inca)) == (GET_SIZE(b, kind = IK) - 1) / max(1_IK,abs(incb)) .or. (GET_SIZE(a, kind = IK) == 1_IK .and. inca == 0_IK) .or. (GET_SIZE(b, kind = IK) == 1_IK .and. incb == 0_IK), \
            SK_"@setSwapped(): The condition `(size(a(1::max(1,abs(inca)))) == size(b(1::max(1,abs(incb))))) .or. (size(a) == 1 .and. inca == 0) .or. (size(b) == 1 .and. incb == 0)` must hold. size(a), inca, size(b), incb = "\
            //getStr([GET_SIZE(a, kind = IK), INCA, GET_SIZE(b, kind = IK), INCB])) ! fpp
            if (inca /= 0_IK) then
                nell = 1 + (GET_SIZE(a, kind = IK) - 1) / abs(inca)
            elseif (incb /= 0_IK) then
                nell = 1 + (GET_SIZE(b, kind = IK) - 1) / abs(incb)
            elseif (inca == 0_IK .and. incb == 0_IK) then
                return
            else
                nell = 1_IK
            end if
#if         BLAS_ENABLED && DISPATCH_ENABLED
            call blasSWAP(nell, a, inca, b, incb)
#else
            block
                integer(IK) :: iell, aell, bell
                aell = 1
                bell = 1
                if (inca < 0_IK) aell = (1 - nell) * inca + 1
                if (incb < 0_IK) bell = (1 - nell) * incb + 1
                do iell = 1, nell
                    tmp = a(GET_SLICE(aell))
                    a(GET_SLICE(aell)) = b(GET_SLICE(bell))
                    b(GET_SLICE(bell)) = tmp
                    aell = aell + inca
                    bell = bell + incb
                end do
            end block
#endif
        end if
#else
#error  "Unrecognized interface."
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  CHECK_STRLEN
#undef  GET_SLICE
#undef  GET_SIZE