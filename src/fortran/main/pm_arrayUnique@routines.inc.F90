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
!>  This file contains the implementations of the procedures in [pm_arrayUnique](@ref pm_arrayUnique).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the equivalence checks.
#if     LK_ENABLED
#define IS_NEQ .neqv.
#define IS_EQ .eqv.
#else
#define IS_NEQ /=
#define IS_EQ ==
#endif
        ! Define the comparison operation.
#if     DefCom_ENABLED
#define NOT_COMPARABLE(i, j) i IS_NEQ j
#define COMPARABLE(i, j) i IS_EQ j
#elif   CusCom_ENABLED
#define NOT_COMPARABLE(i, j) .not. iseq(i, j)
#define COMPARABLE(i, j) iseq(i, j)
#else
#error  "Unrecognized interface."
#endif
        ! Define the indexing and length rules.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE(x) len(x, kind = IK)
#elif   D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE(x) size(x, kind = IK)
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%
#if     isUnique_ENABLED
        !%%%%%%%%%%%%%%%

        integer(IK) :: iarray, jarray, lenArray
        lenArray = GET_SIZE(array)
        unique = .true._LK
        loopOuter: do iarray = 1, lenArray - 1
            if (unique(iarray)) then
                loopInner: do jarray = iarray + 1, lenArray
                    if (COMPARABLE(array(GET_INDEX(iarray)), array(GET_INDEX(jarray)))) then
                        unique(iarray) = .false._LK
                        unique(jarray) = .false._LK
                    end if
                end do loopInner
            end if
        end do loopOuter

        !%%%%%%%%%%%%%%%%%%
#elif   isUniqueAll_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: iarray, jarray, lenArray
        lenArray = GET_SIZE(array)
        uniqueAll = .true._LK
        loopOuter: do iarray = 1, lenArray - 1
            loopInner: do jarray = iarray + 1, lenArray
                if (NOT_COMPARABLE(array(GET_INDEX(iarray)), array(GET_INDEX(jarray)))) cycle loopInner
                uniqueAll = .false._LK
                return
            end do loopInner
        end do loopOuter

        !%%%%%%%%%%%%%%%%%%
#elif   isUniqueAny_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK) :: iarray, jarray, lenArray
        lenArray = GET_SIZE(array)
        uniqueAny = .false._LK
        loopOuter: do iarray = 1, lenArray
            loopInner1: do jarray = 1, iarray - 1
                if (COMPARABLE(array(GET_INDEX(iarray)), array(GET_INDEX(jarray)))) cycle loopOuter
            end do loopInner1
            loopInner2: do jarray = iarray + 1, lenArray
                if (COMPARABLE(array(GET_INDEX(iarray)), array(GET_INDEX(jarray)))) cycle loopOuter
            end do loopInner2
            uniqueAny = .true._LK
            return
        end do loopOuter

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getUnique_ENABLED || setUnique_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setUnique_ENABLED
        integer(IK) , allocatable :: countIndexSorted(:)
        integer(IK) :: counter
#endif
        logical(LK) :: equivalent
        integer(IK) :: iarray, iuniq, lenArray
#if     !(setUnique_ENABLED && UniFix_ENABLED)
        integer(IK) :: lenUnique
        allocate(unique, mold = array)
#endif
        lenArray = GET_SIZE(array) ! fpp
#if     setUnique_ENABLED && UniArb_ENABLED
        allocate(count(lenArray))
#elif   !((setUnique_ENABLED && UniFix_ENABLED) || getUnique_ENABLED)
#error  "Unrecognized interface."
#endif
        lenUnique = 0_IK
        loopOverArray: do iarray = 1_IK, lenArray
            equivalent = .false._LK
            loopOverUnique: do iuniq = 1_IK, lenUnique
                equivalent = COMPARABLE(array(GET_INDEX(iarray)), unique(GET_INDEX(iuniq))) ! fpp
                if (equivalent) then
#if                 setUnique_ENABLED
                    count(iuniq) = count(iuniq) + 1_IK
#endif
                    exit loopOverUnique
                end if
            end do loopOverUnique
            if (.not. equivalent) then
                lenUnique = lenUnique + 1_IK
#if             setUnique_ENABLED && UniFix_ENABLED
                CHECK_ASSERTION(__LINE__, lenUnique <= size(count, 1, IK), SK_"@setUnique(): The condition `lenUnique <= size(count)` must hold. lenUnique, size(count) = "//getStr([lenUnique, size(count, 1, IK)]))
                CHECK_ASSERTION(__LINE__, lenUnique <= GET_SIZE(unique), SK_"@setUnique(): The condition `lenUnique <= len/size(count)` must hold. lenUnique, len/size(unique) = "//getStr([lenUnique, GET_SIZE(unique)]))
#elif           !((setUnique_ENABLED && UniArb_ENABLED) || getUnique_ENABLED)
#error          "Unrecognized interface."
#endif
                unique(GET_INDEX(lenUnique)) = array(GET_INDEX(iarray))
#if             setUnique_ENABLED
                count(lenUnique) = 1_IK
#endif
            end if
        end do loopOverArray
#if     setUnique_ENABLED
#endif
#if     !(setUnique_ENABLED && UniFix_ENABLED)
        unique = unique(1:lenUnique)
#endif
        ! This section is relevant only to the subroutine interfaces, to compute the count and index of unique elements.
#if     setUnique_ENABLED
#if     UniArb_ENABLED
        count = count(1:lenUnique)
#endif
        if (present(order)) then
            if (order /= 0_IK) then
                allocate(countIndexSorted(lenUnique))
                call setSorted(count(1:lenUnique), countIndexSorted)
                if (order > 0_IK) then
#if                 UniArb_ENABLED
                    call setRemapped(count  , countIndexSorted)
                    call setRemapped(unique , countIndexSorted)
#elif               UniFix_ENABLED
                    count(1:lenUnique)  = count(countIndexSorted(1:lenUnique))
#if                 D0_ENABLED && SK_ENABLED
                    unique(1:lenUnique) = getRemapped(unique(1:lenUnique) , countIndexSorted)
#elif               D1_ENABLED
                    !   \bug Intel ifort bug: `getRemapped()` cannot assign `unique(1:lenUnique)` correctly.
                    unique(1:lenUnique) = unique(countIndexSorted(1:lenUnique))
#else
#error              "Unrecognized interface."
#endif
#endif
                else
#if                 UniArb_ENABLED
                    call setRemapped(count  , countIndexSorted, action = reverse)
                    call setRemapped(unique , countIndexSorted, action = reverse)
#elif               UniFix_ENABLED
                    count(1:lenUnique)  = count(countIndexSorted(lenUnique:1:-1))
#if                 D0_ENABLED && SK_ENABLED
                    unique(1:lenUnique) = getRemapped(unique(1:lenUnique), countIndexSorted, action = reverse)
#elif               D1_ENABLED
                    !   \bug Intel ifort bug: `getRemapped()` cannot assign `unique(1:lenUnique)` correctly.
                    unique(1:lenUnique) = unique(countIndexSorted(lenUnique:1:-1))
#else
#error              "Unrecognized interface."
#endif
#endif
                end if
                deallocate(countIndexSorted)
            end if
        end if

        if (present(index)) then
#if         UniArb_ENABLED
            allocate(index(lenUnique))
#elif       UniFix_ENABLED
            CHECK_ASSERTION(__LINE__, lenUnique <= size(index, 1, IK), SK_"@setUnique(): The condition `lenUnique <= size(index)` must hold. lenUnique, size(index) = "//getStr([lenUnique, size(index, 1, IK)]))
#else
#error      "Unrecognized interface."
#endif
            loopOverUniqueForIndex: do iuniq = 1_IK, lenUnique
                allocate(index(iuniq)%val(count(iuniq)))
                iarray = 0_IK
                counter = 1_IK
                loopOverArrayForIndex: do
                    iarray = iarray + 1_IK
                    if (NOT_COMPARABLE(unique(GET_INDEX(iuniq)), array(GET_INDEX(iarray)))) cycle loopOverArrayForIndex ! fpp
                    index(iuniq)%val(counter) = iarray
                    counter = counter + 1_IK
                    if (counter > count(iuniq)) exit loopOverArrayForIndex
                end do loopOverArrayForIndex
            end do loopOverUniqueForIndex
        end if
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  NOT_COMPARABLE
#undef  COMPARABLE
#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_NEQ
#undef  IS_EQ
#undef  ALL