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
!>  This file contains procedures implementations of the module [pm_arrayComplement](@ref pm_arrayComplement).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%
#if     getCompRange_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        use pm_arrayRange, only: getRange
        integer(IKC)    :: i, j
        integer(IKC)    :: lenSetA
        integer(IKC)    :: lenComplement
        integer(IKC)    :: complementTemp(max(0_IKC, 1_IKC + floor(real(stop - start) / real(step), kind = IKC)))
#if     Sorted_ENABLED
        integer(IKC)    :: jstart
#elif   !Random_ENABLED
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, step /= 0_IKC, SK_"@getCompRange(): The input `step` must be non-zero. step = "//getStr(step)) ! fpp
        lenComplement = 0_IKC
        lenSetA = size(setA, kind = IKC)
        if (lenSetA > 0_IK) then
#if     Sorted_ENABLED
            if (sorted) then
                jstart = 1_IKC
                if (unique) then
                    if (step == 1_IKC .and. start <= setA(1) .and. setA(lenSetA) <= stop) then
                        allocate(complement(abs(stop - start) + 1_IKC - lenSetA))
                        loopOverSuperSetOrderedUnique: do i = start, stop, step
                            loopOverSubSetOrderedUnique: do j = jstart, lenSetA
                                if (i == setA(j)) then
                                    jstart = j + 1_IKC
                                    cycle loopOverSuperSetOrderedUnique
                                end if
                            end do loopOverSubSetOrderedUnique
                            lenComplement = lenComplement + 1_IKC
                            complement(lenComplement) = i
                        end do loopOverSuperSetOrderedUnique
                    else
                        loopOverRangeOrderedUnique: do i = start, stop, step
                            loopOverSetOrderedUnique: do j = jstart, lenSetA
                                if (i == setA(j)) then
                                    jstart = j + 1_IKC
                                    cycle loopOverRangeOrderedUnique
                                end if
                            end do loopOverSetOrderedUnique
                            lenComplement = lenComplement + 1_IKC
                            complementTemp(lenComplement) = i
                        end do loopOverRangeOrderedUnique
                        complement = complementTemp(1:lenComplement)
                    end if
                else
                    loopOverRangeOrdered: do i = start, stop, step
                        loopOverSetOrdered: do j = jstart, lenSetA
                            if (i == setA(j)) then
                                jstart = j
                                do
                                    jstart = jstart + 1_IKC
                                    if (jstart > lenSetA) then
                                        do jstart = i + step, stop, step
                                            !print *, i, step, i + step, jstart, i == step
                                            !print *, start, stop, step, size(setA), size(complementTemp), lenComplement
                                            !print *, complementTemp
                                            !print *, "setA"
                                            !print *, setA
                                            !print *, "setA"
                                            lenComplement = lenComplement + 1_IKC
                                            complementTemp(lenComplement) = jstart
                                        end do
                                        exit loopOverRangeOrdered
                                    end if
                                    if (i /= setA(jstart)) cycle loopOverRangeOrdered
                                end do
                                cycle loopOverRangeOrdered
                            end if
                        end do loopOverSetOrdered
                        lenComplement = lenComplement + 1_IKC
                        complementTemp(lenComplement) = i
                    end do loopOverRangeOrdered
                    complement = complementTemp(1:lenComplement)
                end if
            else
#endif
                loopOverRange: do i = start, stop, step
                    loopOverSet: do j = 1_IKC, lenSetA
                        if (i == setA(j)) then
                            cycle loopOverRange
                        end if
                    end do loopOverSet
                    lenComplement = lenComplement + 1_IKC
                    complementTemp(lenComplement) = i
                end do loopOverRange
                complement = complementTemp(1:lenComplement)
#if             Sorted_ENABLED
            end if
#endif
        else
            complement = getRange(start, stop, step)
        end if

        !%%%%%%%%%%%%%%%%%%%%
#elif   getComplement_ENABLED
        !%%%%%%%%%%%%%%%%%%%%

        ! Define the equivalence checking method.
#if     DefCom_ENABLED && LK_ENABLED
#define ISEQ(elementA,elementB) elementA .eqv. elementB
#elif   DefCom_ENABLED
#define ISEQ(elementA,elementB) elementA == elementB
#elif   CusCom_ENABLED
#define ISEQ(elementA, elementB) iseq(elementA, elementB)
#else
#error  "Unrecognized interface."
#endif
        ! Define temporary complement storage.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
        character(len(setB,IK),SKC) :: complementTemp
#elif   D1_ENABLED
#define GET_SIZE size
#define GET_INDEX(i) i
#if     SK_ENABLED
        character(len(setB,IK),SKC) &
#elif   IK_ENABLED
        integer(IKC) &
#elif   LK_ENABLED
        logical(LKC) &
#elif   CK_ENABLED
        complex(CKC) &
#elif   RK_ENABLED
        real(RKC) &
#else
#error  "Unrecognized interface."
#endif
        & :: complementTemp(size(setB))
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: lenComplement
        integer(IK) :: lenSetA
        integer(IK) :: lenSetB
        integer(IK) :: i, j

        ! Define the handling method of sorted vs. unsorted sets.
#if     Random_ENABLED
        integer(IK) , parameter :: jstart = 1_IK
#define INCREMENT(jstart)
#elif   Sorted_ENABLED
#define INCREMENT(jstart) jstart = j + 1_IK
        integer(IK) :: jstart
        jstart = 1_IK
        if (sorted) then
            if (unique) then
#else
#error          "Unrecognized interface."
#endif
                lenSetA = GET_SIZE(setA, kind = IK) ! fpp
                lenSetB = GET_SIZE(setB, kind = IK) ! fpp
                lenComplement = 0_IK
                loopOverUniqueSetB: do i = 1_IK, lenSetB
                    loopOverUniqueSetA: do j = jstart, lenSetA
                        if (ISEQ(setA(GET_INDEX(j)), setB(GET_INDEX(i)))) then ! fpp
                            INCREMENT(jstart) ! fpp
                            cycle loopOverUniqueSetB
                        end if
                    end do loopOverUniqueSetA
                    lenComplement = lenComplement + 1_IK
                    complementTemp(GET_INDEX(lenComplement)) = setB(GET_INDEX(i))
                end do loopOverUniqueSetB
                complement = complementTemp(1:lenComplement)
#if         Sorted_ENABLED
            else ! sorted but not unique.
                lenSetA = GET_SIZE(setA, kind = IK) ! fpp
                lenSetB = GET_SIZE(setB, kind = IK) ! fpp
                lenComplement = 0_IK
                if (lenSetB > 0_IK) then
                    loopOverSetB: do i = 1_IK, lenSetB
                        loopOverSetA: do j = jstart, lenSetA
                            if (ISEQ(setA(GET_INDEX(j)), setB(GET_INDEX(i)))) then ! fpp
                                if (i < lenSetB) then
                                    ! go to the next element in setA only if the next element in setB is not the same as the current element in setB.
                                    if (.not. ISEQ(setB(GET_INDEX(i)), setB(GET_INDEX(i+1)))) INCREMENT(jstart) ! fpp
                                else
                                    INCREMENT(jstart) ! fpp
                                end if
                                cycle loopOverSetB
                            end if
                        end do loopOverSetA
                        lenComplement = lenComplement + 1_IK
                        complementTemp(GET_INDEX(lenComplement)) = setB(GET_INDEX(i))
                    end do loopOverSetB
                end if
                complement = complementTemp(1:lenComplement)
            end if
        else ! not sorted
#if         DefCom_ENABLED
            complement = getComplement(setA, setB)
#elif       CusCom_ENABLED
            complement = getComplement(setA, setB, iseq)
#else
#error      "Unrecognized interface."
#endif
        end if
#endif

#undef  INCREMENT
#undef  GET_INDEX
#undef  GET_SIZE
#undef  ISEQ

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
