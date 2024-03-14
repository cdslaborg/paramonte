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
!>  This include file contains the procedure implementation of Non-Recursive QuickSort selecting the smallest `rank`th value in the input array.
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the runtime debugging parameters.
#if     CHECK_ENABLED
#if     getSelected_ENABLED
        character(*,SK), parameter :: PROCEDURE_NAME = SK_"@getSelected"
#elif   setSelected_ENABLED
        character(*,SK), parameter :: PROCEDURE_NAME = SK_"@setSelected"
#else
#error  "Unrecognized interface."
#endif
#endif
        ! Define the auxiliary variables.
#if     SK_ENABLED && D0_ENABLED
        character(1,SKC) :: pivot, temp
#elif   SK_ENABLED && D1_ENABLED
        character(len(array),SKC) :: pivot, temp
#elif   IK_ENABLED && D1_ENABLED
        integer(IKC) :: pivot, temp
#elif   LK_ENABLED && D1_ENABLED
        logical(LKC) :: pivot, temp
#elif   CK_ENABLED && D1_ENABLED
        complex(CKC) :: pivot, temp
#elif   RK_ENABLED && D1_ENABLED
        real(RKC) :: pivot, temp
#elif   PSSK_ENABLED && D1_ENABLED
        type(css_pdt(SKC)) :: pivot, temp
#elif   BSSK_ENABLED && D1_ENABLED
        type(css_type) :: pivot, temp
#else
#error  "Unrecognized interface."
#endif
        ! Define the auxiliary variables for the functional interface.
#if     getSelected_ENABLED
#if     SK_ENABLED && D0_ENABLED
        character(len(array),SKC) :: arrayCopy
#elif   SK_ENABLED && D1_ENABLED
        character(len(array),SKC) :: arrayCopy(size(array))
#elif   IK_ENABLED && D1_ENABLED
        integer(IKC) :: arrayCopy(size(array))
#elif   LK_ENABLED && D1_ENABLED
        logical(LKC) :: arrayCopy(size(array))
#elif   CK_ENABLED && D1_ENABLED
        complex(CKC) :: arrayCopy(size(array))
#elif   RK_ENABLED && D1_ENABLED
        real(RKC) :: arrayCopy(size(array))
#elif   PSSK_ENABLED && D1_ENABLED
        type(css_pdt(SKC)) :: arrayCopy(size(array))
#elif   BSSK_ENABLED && D1_ENABLED
        type(css_type) :: arrayCopy(size(array))
#else
#error  "Unrecognized interface."
#endif
#elif   !setSelected_ENABLED
#error  "Unrecognized interface."
#endif
        ! Set the custom vs. default sorting criterion.
#if     CusCom_ENABLED
#define IS_SORTED(i,j) isSorted(i,j)
#elif   DefCom_ENABLED && D1_ENABLED && (PSSK_ENABLED || BSSK_ENABLED)
#define IS_SORTED(i,j) i%val < j%val
#elif   DefCom_ENABLED && D1_ENABLED && LK_ENABLED
#define IS_SORTED(i,j) j .and. .not. i
#elif   DefCom_ENABLED && D1_ENABLED && CK_ENABLED
#define IS_SORTED(i,j) i%re < j%re
#elif   DefCom_ENABLED
#define IS_SORTED(i,j) i < j
#else
#error  "Unrecognized interface."
#endif
        ! Define the indexing rules.
#if     D0_ENABLED && SK_ENABLED
#define GET_SIZE(array) len(array, kind = IK)
#define GET_INDEX(i) i:i
#elif   D1_ENABLED && (SK_ENABLED || IK_ENABLED || LK_ENABLED || CK_ENABLED || RK_ENABLED || PSSK_ENABLED || BSSK_ENABLED)
#define GET_SIZE(array) size(array, kind = IK)
#define GET_INDEX(i) i
#else
#error  "Unrecognized interface."
#endif
#if     indexing_ENABLED
        integer(IK), allocatable :: arrayIndex(:)
#define SELECTION index
#endif
        integer(IK) :: mid, start, low, high, lenArray
        lenArray = GET_SIZE(array)
        if (present(lb)) then
            low = lb
            CHECK_ASSERTION(__LINE__, 1_IK <= lb, PROCEDURE_NAME//SK_": The condition `1 <= lb` must hold. lb = "//getStr(lb))
        else
            low = 1_IK
        end if
        if (present(ub)) then
            high = ub
            CHECK_ASSERTION(__LINE__, ub <= lenArray, PROCEDURE_NAME//SK_": The condition `ub <= lenArray` must hold. ub, lenArray = "//getStr([ub, lenArray]))
        else
            high = lenArray
        end if
        ! This condition together with the previous ones also guarantees that the input array length is non-zero.
        CHECK_ASSERTION(__LINE__, low <= rank .and. rank <= high, PROCEDURE_NAME//SK_": The condition `low <= rank .and. rank <= high` must hold. low, rank, high = "//getStr([low, rank, high]))
#if     indexing_ENABLED
#define GET_VALUE(i) arrayIndex(i)
#define PIVOT array(GET_INDEX(pivot))
        allocate(arrayIndex(low:high))
        do concurrent(mid = low:high)
            arrayIndex(mid) = mid
        end do
#elif   getSelected_ENABLED
#define ARRAY arrayCopy
#define GET_VALUE(i) ARRAY(GET_INDEX(i))
        arrayCopy(low:high) = array(low:high)
#elif   setSelected_ENABLED
#define GET_VALUE(i) ARRAY(GET_INDEX(i))
#else
#error  "Unrecognized interface."
#endif
        do
            if (high - low <= 1_IK) then
                if (high - low == 1_IK) then
                    if (IS_SORTED(ARRAY(GET_INDEX(high)), ARRAY(GET_INDEX(low)))) then
                        temp = GET_VALUE(low)
                        GET_VALUE(low) = GET_VALUE(high)
                        GET_VALUE(high) = temp
                    end if
                end if
                SELECTION = GET_VALUE(rank) ! \warning `SELECTION` is a preprocessor macro.
                return
            else
                mid = (low + high) / 2_IK
                temp = GET_VALUE(mid)
                GET_VALUE(mid) = GET_VALUE(low+1_IK)
                GET_VALUE(low+1_IK) = temp
                if (IS_SORTED(ARRAY(GET_INDEX(high)), ARRAY(GET_INDEX(low)))) then
                    temp = GET_VALUE(low)
                    GET_VALUE(low) = GET_VALUE(high)
                    GET_VALUE(high) = temp
                end if
                if (IS_SORTED(ARRAY(GET_INDEX(high)), ARRAY(GET_INDEX(low+1_IK)))) then
                    temp = GET_VALUE(low+1_IK)
                    GET_VALUE(low+1_IK) = GET_VALUE(high)
                    GET_VALUE(high) = temp
                end if
                if (IS_SORTED(ARRAY(GET_INDEX(low+1_IK)), ARRAY(GET_INDEX(low)))) then
                    temp = GET_VALUE(low)
                    GET_VALUE(low) = GET_VALUE(low+1_IK)
                    GET_VALUE(low+1_IK) = temp
                end if
                start = high
                mid = low + 1_IK
                pivot = GET_VALUE(mid)
                do
                    do
                        mid = mid + 1_IK
                        if (IS_SORTED(ARRAY(GET_INDEX(mid)), PIVOT)) cycle ! fpp
                        exit
                    end do
                    do
                        start = start - 1_IK
                        if (IS_SORTED( PIVOT, ARRAY(GET_INDEX(start)))) cycle
                        exit
                    end do
                    if (start < mid) exit
                    temp = GET_VALUE(mid)
                    GET_VALUE(mid) = GET_VALUE(start)
                    GET_VALUE(start) = temp
                end do
                GET_VALUE(low + 1_IK) = GET_VALUE(start)
                GET_VALUE(start) = pivot
                if (start >= rank) high = start - 1_IK
                if (start <= rank) low = mid
            end if
        end do
#undef  indexing_ENABLED
#undef  SELECTION
#undef  GET_VALUE
#undef  GET_INDEX
#undef  IS_SORTED
#undef  GET_SIZE
#undef  PIVOT
#undef  ARRAY