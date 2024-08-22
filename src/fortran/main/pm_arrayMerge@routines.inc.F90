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
!>  This include file contcounterArray1ns the procedure implementation of [pm_arrayMerge](@ref pm_arrayMerge).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Object components.
#if     CK_ENABLED
#define COMPONENT %re
#elif   PSSK_ENABLED || BSSK_ENABLED
#define COMPONENT %val
#else
#define COMPONENT
#endif
        ! Array subsetting.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        ! Logical comparison.
#if     LK_ENABLED
        use pm_logicalCompare, only: operator(<)
#endif
        integer(IK) :: counterArray1, counterArray2, counterMergedArray, i
        integer(IK) :: lenSortedArray1, lenSortedArray2
        lenSortedArray1 = GET_SIZE(sortedArray1, kind = IK)
        lenSortedArray2 = GET_SIZE(sortedArray2, kind = IK)
#if     setMerged_ENABLED
        CHECK_ASSERTION(__LINE__, logical(GET_SIZE(mergedSortedArray, kind = IK) == lenSortedArray1 + lenSortedArray2, LK), \
        SK_"@setMerged(): The output array size must equal the sum of the input array sizes. size(sortedArray1), size(sortedArray2), size(mergedSortedArray) = "//\
        getStr([lenSortedArray1, lenSortedArray2, GET_SIZE(mergedSortedArray, kind = IK)])) ! fpp
#endif
        !   \todo
        !   This runtime check must be extended to container arrays.
        !   Currently these tests are only performed for non-container arrays because `getStr()` cannot handle container arrays.
        !   This must be fixed in the future.
#if     DefCom_ENABLED 
#if     !(PSSK_ENABLED || BSSK_ENABLED)
        CHECK_ASSERTION(__LINE__, isAscending(sortedArray1), SK_": The input argument `sortedArray1` must be ascending sorted. sortedArray1 = "//getStr(sortedArray1)) ! fpp
        CHECK_ASSERTION(__LINE__, isAscending(sortedArray2), SK_": The input argument `sortedArray2` must be ascending sorted. sortedArray2 = "//getStr(sortedArray2)) ! fpp
#endif
#define IS_SORTED(i, j) i < j
#elif   CusCom_ENABLED
#define IS_SORTED(i, j) isSorted(i, j)
        !   \todo
        !   This runtime check must be extended to container arrays.
        !   Custom check can become problematic when `isSorted` is passed from `pm_arraySort`. This may further look in the future.
#if     !(PSSK_ENABLED || BSSK_ENABLED)
        CHECK_ASSERTION(__LINE__, isSortedCheck(sortedArray1, isSortedEqual), SK_": The input argument `sortedArray1` must be sorted. sortedArray1 = "//getStr(sortedArray1)) ! fpp
        CHECK_ASSERTION(__LINE__, isSortedCheck(sortedArray2, isSortedEqual), SK_": The input argument `sortedArray2` must be sorted. sortedArray2 = "//getStr(sortedArray2)) ! fpp
#endif
#else
#error  "Unrecognized interface."
#endif
        counterArray1 = 1_IK
        counterArray2 = 1_IK
        counterMergedArray = 1_IK
        do
            if (counterArray1 > lenSortedArray1) then
                do i = counterArray2, lenSortedArray2
                    mergedSortedArray(GET_INDEX(counterMergedArray)) = sortedArray2(GET_INDEX(i))
                    counterMergedArray = counterMergedArray + 1_IK
                end do
                return
            end if
            if (counterArray2 > lenSortedArray2) then
                do i = counterArray1, lenSortedArray1
                    mergedSortedArray(GET_INDEX(counterMergedArray)) = sortedArray1(GET_INDEX(i))
                    counterMergedArray = counterMergedArray + 1_IK
                end do
                return
            end if
            if (IS_SORTED(sortedArray1(GET_INDEX(counterArray1))COMPONENT, sortedArray2(GET_INDEX(counterArray2))COMPONENT)) then
                mergedSortedArray(GET_INDEX(counterMergedArray)) = sortedArray1(GET_INDEX(counterArray1))
                counterArray1 = counterArray1 + 1_IK
            else
                mergedSortedArray(GET_INDEX(counterMergedArray)) = sortedArray2(GET_INDEX(counterArray2))
                counterArray2 = counterArray2 + 1_IK
            end if
            counterMergedArray = counterMergedArray + 1_IK
        end do
#if CHECK_ENABLED && CusCom_ENABLED && !(PSSK_ENABLED || BSSK_ENABLED)
    contains
        function isSortedEqual(lhs, rhs) result(sorted)
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG)        , intent(in) :: lhs, rhs
#elif       SK_ENABLED && D1_ENABLED
            character(*,SKG)        , intent(in) :: lhs, rhs
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)            , intent(in) :: lhs, rhs
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)            , intent(in) :: lhs, rhs
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)            , intent(in) :: lhs, rhs
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)               , intent(in) :: lhs, rhs
#elif       PSSK_ENABLED && D1_ENABLED
            use pm_container, only: css_pdt
            type(css_pdt(SKG)), intent(in) :: lhs, rhs
#elif       BSSK_ENABLED && D1_ENABLED
            use pm_container, only: css_type
            type(css_type), intent(in) :: lhs, rhs
#else
#error      "Unrecognized interface."
#endif
            logical(LK) :: sorted
            sorted = .not. isSorted(rhs, lhs)
        end function
#endif

#undef  COMPONENT
#undef  IS_SORTED
#undef  GET_INDEX
#undef  GET_SIZE