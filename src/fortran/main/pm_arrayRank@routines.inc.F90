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
!>  This include file contains the procedure implementation of [pm_arrayRank](@ref pm_arrayRank).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getRank_ENABLED && DefCom_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Dense_ENABLED
        call setRankDense(rank, array)
#elif   Fractional_ENABLED
        call setRankFractional(rank, array)
#elif   Modified_ENABLED
        call setRankModified(rank, array)
#elif   Ordinal_ENABLED
        call setRankOrdinal(rank, array)
#elif   Standard_ENABLED
        call setRankStandard(rank, array)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getRank_ENABLED && CusCom_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     Dense_ENABLED
        call setRankDense(rank, array, isSorted)
#elif   Fractional_ENABLED
        call setRankFractional(rank, array, isSorted)
#elif   Modified_ENABLED
        call setRankModified(rank, array, isSorted)
#elif   Ordinal_ENABLED
        call setRankOrdinal(rank, array, isSorted)
#elif   Standard_ENABLED
        call setRankStandard(rank, array, isSorted)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%
#elif   setRank_ENABLED
        !%%%%%%%%%%%%%%

        ! Define the sorting method.
#if     CusCom_ENABLED
#define SET_SORTED_INDEX call setSorted(array, index, isSorted)
#define IS_EQ(i,j) .not. (isSorted(i,j) .or. isSorted(j,i))
#elif   DefCom_ENABLED
#define SET_SORTED_INDEX call setSorted(array, index)
#if     PSSK_ENABLED || BSSK_ENABLED
#define IS_EQ(i,j) i%val == j%val
#elif   LK_ENABLED
#define IS_EQ(i,j) j .eqv. i
#elif   CK_ENABLED
#define IS_EQ(i,j) i%re == j%re
#else
#define IS_EQ(i,j) i == j
#endif
#else
#error  "Unrecognized interface."
#endif
        ! Define the indexing method.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#else
#define GET_INDEX(i) i
#define GET_SIZE size
#endif
        ! Define the runtime check.
#define CHECK_LEN_RANK \
CHECK_ASSERTION(__LINE__, GET_SIZE(array, kind = IK) == size(rank, kind = IK), \
SK_"@setRank(): The input `array` and `rank` must be of the same size: "// \
getStr([GET_SIZE(array, kind = IK), size(rank, kind = IK)])) ! fpp
        ! perform ranking.
#if     Dense_ENABLED
        integer, parameter :: TKR = kind(rank)
        integer(TKR), parameter :: INCREMENT = +1_TKR, FIRST = +1_TKR
        integer(TKR) :: index(size(rank, kind = IK))
        integer(TKR) :: i, last, current
        CHECK_LEN_RANK
        last = size(rank, 1, IK)
        if (size(rank, kind = IK) > 0_IK) then
            SET_SORTED_INDEX ! fpp
            current = first
            rank(index(current)) = current
            loopTie: do
                if (current == last) return
                i = current + INCREMENT
                loopTieSegment: do
                    if (IS_EQ(array(GET_INDEX(index(current))), array(GET_INDEX(index(i))))) then
                        rank(index(i)) = rank(index(i - 1_TKR)) ! This is technically the same as `current`.
                        if (i < last) then
                            i = i + INCREMENT
                        else ! happens only if there is a tied segment at the end.
                            return
                        end if
                    else
                        rank(index(i)) = rank(index(i - 1_TKR)) + 1_TKR
                        current = i
                        cycle loopTie
                    end if
                end do loopTieSegment
            end do loopTie
        end if
#elif   Fractional_ENABLED
        integer     , parameter :: TKR = kind(rank) ! Real Kind of rank.
        integer(IK) , parameter :: INCREMENT = +1_IK, FIRST = +1_IK
        integer(IK) :: i, last, current, index(size(rank, kind = IK))
        real(TKR) :: sumRank
        CHECK_LEN_RANK
        last = size(rank, 1, IK)
        if (size(rank, kind = IK) > 0_IK) then
            SET_SORTED_INDEX ! fpp
            current = first
            rank(index(current)) = current
            sumRank = real(current, TKR)
            loopTie: do
                if (current == last) return
                i = current + INCREMENT
                loopTieSegment: do
                    if (IS_EQ(array(GET_INDEX(index(current))), array(GET_INDEX(index(i))))) then
                        sumRank = sumRank + real(i, TKR)
                        if (i < last) then
                            i = i + INCREMENT
                        else ! happens only if there is a tied segment at the end.
                            rank(index(current : i)) = sumRank / real(i - current + 1_IK, TKR)
                            return
                        end if
                    else
                        rank(index(current : i - 1_IK)) = sumRank / real(i - current, TKR)
                        rank(index(i)) = real(i, TKR)
                        sumRank = real(i, TKR)
                        current = i
                        cycle loopTie
                    end if
                end do loopTieSegment
            end do loopTie
        end if
#elif   Modified_ENABLED
        integer, parameter :: TKR = kind(rank)
        integer(TKR), parameter :: INCREMENT = -1_IK, LAST = +1_IK
        integer(TKR) :: i, first, current, index(size(rank, kind = IK))
        CHECK_LEN_RANK
        first = size(rank, 1, TKR)
        if (size(rank, kind = IK) > 0_IK) then
            SET_SORTED_INDEX ! fpp
            current = first
            rank(index(current)) = current
            loopTie: do
                if (current == last) return
                i = current + INCREMENT
                loopTieSegment: do
                    if (IS_EQ(array(GET_INDEX(index(current))), array(GET_INDEX(index(i))))) then
                        rank(index(i)) = current
                        if (i > last) then
                            i = i + INCREMENT
                        else ! happens only if there is a tied segment at the end.
                            return
                        end if
                    else
                        rank(index(i)) = i
                        current = i
                        cycle loopTie
                    end if
                end do loopTieSegment
            end do loopTie
        end if
#elif   Ordinal_ENABLED
        integer, parameter :: TKR = kind(rank)
        integer(TKR) :: i, index(size(rank, kind = IK))
        CHECK_LEN_RANK
        SET_SORTED_INDEX ! fpp
        do concurrent(i = 1_TKR : size(rank, 1, TKR))
            rank(index(i)) = i
        end do
#elif   Standard_ENABLED
        integer, parameter :: TKR = kind(rank)
        integer(TKR), parameter :: INCREMENT = +1_IK, FIRST = +1_IK
        integer(TKR) :: current, i, last, index(size(rank, kind = IK))
        CHECK_LEN_RANK
        last = size(rank, kind = TKR)
        if (size(rank, kind = IK) > 0_IK) then
            SET_SORTED_INDEX ! fpp
            current = first
            rank(index(current)) = current
            loopTie: do
                if (current == last) return
                i = current + INCREMENT
                loopTieSegment: do
                    if (IS_EQ(array(GET_INDEX(index(current))), array(GET_INDEX(index(i))))) then
                        rank(index(i)) = current
                        if (i < last) then
                            i = i + INCREMENT
                        else ! happens only if there is a tied segment at the end.
                            return
                        end if
                    else
                        rank(index(i)) = i
                        current = i
                        cycle loopTie
                    end if
                end do loopTieSegment
            end do loopTie
        end if
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef SET_SORTED_INDEX
#undef GET_INDEX
#undef GET_SIZE
#undef IS_EQ