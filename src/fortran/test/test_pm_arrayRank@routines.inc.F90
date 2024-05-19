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
!>  This module contains implementations of the tests of the procedures of [pm_arrayRank](@ref pm_arrayRank).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the testing functions.
#if     Dense_ENABLED
#define GETRANK getRankDense
#define SETRANK setRankDense
#elif   Fractional_ENABLED
#define GETRANK getRankFractional
#define SETRANK setRankFractional
#elif   Modified_ENABLED
#define GETRANK getRankModified
#define SETRANK setRankModified
#elif   Ordinal_ENABLED
#define GETRANK getRankOrdinal
#define SETRANK setRankOrdinal
#elif   Standard_ENABLED
#define GETRANK getRankStandard
#define SETRANK setRankStandard
#else
#error  "Unrecognized interface."
#endif
        ! Define the indexing rules.
#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE(array) len(array, kind = IK)
#define GET_INDEX(i) i:i
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || CK_ENABLED || RK_ENABLED || PSSK_ENABLED
#define GET_SIZE(array) size(array, kind = IK)
#define GET_INDEX(i) i
#else
#error  "Unrecognized interface."
#endif
        ! Define the comparison operator.
#if     LK_ENABLED && D1_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
#if     Fractional_ENABLED
        real(RKR), allocatable :: rank(:), rank_ref(:)
#else
        integer(IK) , allocatable :: rank(:), rank_ref(:)
#endif
#if     PSSK_ENABLED
        use pm_container, only: strc => css_pdt, operator(==)
#endif
        integer(IK) :: itest, maxLenArray
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG), allocatable :: array, lower, upper
        lower = SKG_"a"; upper = SKG_"z"
#elif   SK_ENABLED && D1_ENABLED
        character(:,SKG), allocatable :: array(:), lower, upper
        lower = SKG_"aa"; upper = SKG_"zz"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , allocatable :: array(:), lower, upper
        lower = 0_IKG; upper = 9_IKG
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , allocatable :: array(:), lower, upper
        lower = .false._LKG; upper = .true._LKG
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , allocatable :: array(:), lower, upper
        lower = -(1._CKG, 1._CKG); upper = (1._CKG, 1._CKG)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , allocatable :: array(:), lower, upper
        lower = 0._RKG; upper = 1._RKG
#elif   PSSK_ENABLED && D1_ENABLED
        type(strc(SKG)) , allocatable :: array(:), lower, upper
        lower = strc(SKG_"a"); upper = strc(SKG_"zzz")
#else
#error  "Unrecognized Interface."
#endif
        assertion = .true._LK

        call runTestsWith()
        call runTestsWith(isSorted)

    contains

        subroutine runTestsWith(isSorted)
            procedure(logical(LK)), optional :: isSorted

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         SK_ENABLED && D0_ENABLED
            array = SKG_""
#elif       SK_ENABLED && D1_ENABLED
            array = [character(2,SKG) ::]
#elif       IK_ENABLED && D1_ENABLED
            array = [integer(IKG) ::]
#elif       LK_ENABLED && D1_ENABLED
            array = [logical(LKG) ::]
#elif       CK_ENABLED && D1_ENABLED
            array = [complex(CKG) ::]
#elif       RK_ENABLED && D1_ENABLED
            array = [real(RKG) ::]
#elif       PSSK_ENABLED && D1_ENABLED
#else
#error      "Unrecognized Interface."
#endif
            allocate(rank_ref(0))
            call report(int(__LINE__, IK), isSorted)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         SK_ENABLED && D0_ENABLED
            array = SKG_"1122"
#elif       SK_ENABLED && D1_ENABLED
            array = [character(2,SKG) :: "1", "1", "2", "2"]
#elif       IK_ENABLED && D1_ENABLED
            array = [integer(IKG) :: 1, 1, 2, 2]
#elif       LK_ENABLED && D1_ENABLED
            array = [logical(LKG) :: .false., .false., .true., .true.]
#elif       CK_ENABLED && D1_ENABLED
            array = [complex(CKG) :: (1., -1.), (1., -1.), (2., -2.), (2., -2.)]
#elif       RK_ENABLED && D1_ENABLED
            array = [real(RKG) :: 1., 1., 2., 2.]
#else
#error      "Unrecognized Interface."
#endif
#if         Dense_ENABLED
            rank_ref = [integer(IK) :: 1, 1, 2, 2]
            if (present(isSorted)) call setReversed(rank_ref)
#elif       Fractional_ENABLED
            rank_ref = [real(RK) :: 1.5, 1.5, 3.5, 3.5]
            if (present(isSorted)) call setReversed(rank_ref)
#elif       Modified_ENABLED
            rank_ref = [integer(IK) :: 2, 2, 4, 4]
            if (present(isSorted)) call setReversed(rank_ref)
#elif       Ordinal_ENABLED
            rank_ref = [integer(IK) :: 1, 2, 3, 4]
            if (present(isSorted)) rank_ref = [integer(IK) :: 3, 4, 1, 2]
#elif       Standard_ENABLED
            rank_ref = [integer(IK) :: 1, 1, 3, 3]
            if (present(isSorted)) call setReversed(rank_ref)
#else
#error      "Unrecognized Interface."
#endif
            call report(int(__LINE__, IK), isSorted)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         SK_ENABLED && D0_ENABLED
            array = SKG_"123454"
#elif       SK_ENABLED && D1_ENABLED
            array = [character(2,SKG) :: "1", "2", "3", "4", "5", "4"]
#elif       IK_ENABLED && D1_ENABLED
            array = [integer(IKG) :: 1, 2, 3, 4, 5, 4]
#elif       LK_ENABLED && D1_ENABLED
            array = [logical(LKG) :: .false., .false., .true., .true., .false., .true.]
#elif       CK_ENABLED && D1_ENABLED
            array = [complex(CKG) :: (1., -1.), (2., -2.), (3., -3.), (4., -4.), (5., -5.), (4., -4.)]
#elif       RK_ENABLED && D1_ENABLED
            array = [real(RKG) :: 1., 2., 3., 4., 5., 4.]
#else
#error      "Unrecognized Interface."
#endif
            rank_ref = getRank(array, isSorted)
            call report(int(__LINE__, IK), isSorted)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            maxLenArray = 50_IK
            do itest = 1, 100
                call reset()
#if             SK_ENABLED && D0_ENABLED
                array = getUnifRand(repeat(lower, len(array)), repeat(upper, len(array)))
#else
                array = getUnifRand(lower, upper, getUnifRand(0, maxLenArray))
#endif
                rank_ref = getRank(array, isSorted)
                call report(int(__LINE__, IK), isSorted)
            end do

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        subroutine reset()
            if (allocated(array)) deallocate(array)
            if (allocated(rank_ref)) deallocate(rank_ref)
        end subroutine

        function getRank(array, isSorted) result(rank)
            procedure(logical(LK))  , optional  :: isSorted
            logical(LK) :: isSortedPresent
#if         SK_ENABLED && D0_ENABLED
            character(*,SKG) :: array
#elif       SK_ENABLED && D1_ENABLED
            character(*,SKG), contiguous :: array(:)
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG), contiguous :: array(:)
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG), contiguous :: array(:)
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG), contiguous :: array(:)
#elif       RK_ENABLED && D1_ENABLED
            real(RKG), contiguous :: array(:)
#elif       PSSK_ENABLED && D1_ENABLED
            type(strc(SKG)), contiguous :: array(:)
#else
#error      "Unrecognized Interface."
#endif
#if         Dense_ENABLED
            integer(IK) :: rank(GET_SIZE(array))
            integer(IK) :: index(GET_SIZE(array))
            integer(IK) :: i, currentRank
            isSortedPresent = logical(present(isSorted), LK) ! ifort 2021.1 bug: ifort cannot recognize `present`.
            if (isSortedPresent) then
                call setSorted(array, index, isSorted)
            else
                call setSorted(array, index)
            end if
            i = 0_IK
            currentRank = 0_IK
            loopi: do
                i = i + 1_IK
                if (i > GET_SIZE(array)) exit
                currentRank = currentRank + 1_IK
                rank(index(i)) = currentRank
                block
                    integer(IK) :: j
                    logical(LK) :: equivalent
                    do j = i + 1_IK, GET_SIZE(array)
                        if (isSortedPresent) then
                            equivalent = .not. isSorted(array(GET_INDEX(index(i))), array(GET_INDEX(index(j)))) .and. .not. isSorted(array(GET_INDEX(index(j))), array(GET_INDEX(index(i))))
                        else
                            equivalent = logical(array(GET_INDEX(index(i))) IS_EQUAL array(GET_INDEX(index(j))), LK)
                        end if
                        if (.not. equivalent) then
                            rank(index(i+1:j-1)) = currentRank
                            i = j - 1_IK
                            cycle loopi
                        end if
                    end do
                    rank(index(i+1:)) = currentRank
                end block
                exit loopi
            end do loopi
#elif       Fractional_ENABLED
            real(RKR)   :: rank(GET_SIZE(array))
            integer(IK) :: index(GET_SIZE(array))
            integer(IK) :: i, currentRank
            isSortedPresent = logical(present(isSorted), LK) ! ifort 2021.1 bug: ifort cannot recognize `present`.
            if (isSortedPresent) then
                call setSorted(array, index, isSorted)
            else
                call setSorted(array, index)
            end if
            i = 0_IK
            currentRank = 0_IK
            loopi: do
                i = i + 1_IK
                if (i > GET_SIZE(array)) exit
                currentRank = currentRank + 1_IK
                block
                    integer(IK) :: j
                    logical(LK) :: equivalent
                    do j = i + 1_IK, GET_SIZE(array)
                        if (isSortedPresent) then
                            equivalent = .not. isSorted(array(GET_INDEX(index(i))), array(GET_INDEX(index(j)))) .and. .not. isSorted(array(GET_INDEX(index(j))), array(GET_INDEX(index(i))))
                        else
                            equivalent = logical(array(GET_INDEX(index(i))) IS_EQUAL array(GET_INDEX(index(j))), LK)
                        end if
                        if (.not. equivalent) then
                            rank(index(i:j-1)) = real(sum(getRange(i,j-1)), RKR) / real(j - i, RKR)
                            i = j - 1_IK
                            cycle loopi
                        end if
                    end do
                    rank(index(i:)) = real(sum(getRange(i,j-1)), RKR) / real(j - i, RKR)
                end block
                exit loopi
            end do loopi
#elif       Modified_ENABLED
            integer(IK) :: rank(GET_SIZE(array))
            integer(IK) :: index(GET_SIZE(array))
            integer(IK) :: i, currentRank
            isSortedPresent = logical(present(isSorted), LK) ! ifort 2021.1 bug: ifort cannot recognize `present`.
            if (isSortedPresent) then
                call setSorted(array, index, isSorted)
            else
                call setSorted(array, index)
            end if
            i = 0_IK
            currentRank = 0_IK
            loopi: do
                i = i + 1_IK
                if (i > GET_SIZE(array)) exit
                currentRank = currentRank + 1_IK
                block
                    integer(IK) :: j
                    logical(LK) :: equivalent
                    do j = i + 1_IK, GET_SIZE(array)
                        if (isSortedPresent) then
                            equivalent = .not. isSorted(array(GET_INDEX(index(i))), array(GET_INDEX(index(j)))) .and. .not. isSorted(array(GET_INDEX(index(j))), array(GET_INDEX(index(i))))
                        else
                            equivalent = logical(array(GET_INDEX(index(i))) IS_EQUAL array(GET_INDEX(index(j))), LK)
                        end if
                        if (.not. equivalent) then
                            rank(index(i:j-1)) = j - 1_IK
                            i = j - 1_IK
                            cycle loopi
                        end if
                    end do
                    rank(index(i:)) = size(index)
                end block
                exit loopi
            end do loopi
#elif       Standard_ENABLED
            integer(IK) :: rank(GET_SIZE(array))
            integer(IK) :: index(GET_SIZE(array))
            integer(IK) :: i, currentRank
            isSortedPresent = logical(present(isSorted), LK) ! ifort 2021.1 bug: ifort cannot recognize `present`.
            if (isSortedPresent) then
                call setSorted(array, index, isSorted)
            else
                call setSorted(array, index)
            end if
            i = 0_IK
            currentRank = 0_IK
            loopi: do
                i = i + 1_IK
                if (i > GET_SIZE(array)) exit
                currentRank = currentRank + 1_IK
                block
                    integer(IK) :: j
                    logical(LK) :: equivalent
                    do j = i + 1_IK, GET_SIZE(array)
                        if (isSortedPresent) then
                            equivalent = .not. isSorted(array(GET_INDEX(index(i))), array(GET_INDEX(index(j)))) .and. .not. isSorted(array(GET_INDEX(index(j))), array(GET_INDEX(index(i))))
                        else
                            equivalent = logical(array(GET_INDEX(index(i))) IS_EQUAL array(GET_INDEX(index(j))), LK)
                        end if
                        if (.not. equivalent) then
                            rank(index(i:j-1)) = currentRank
                            currentRank = j - 1
                            i = j - 1_IK
                            cycle loopi
                        end if
                    end do
                    rank(index(i:)) = currentRank
                end block
                exit loopi
            end do loopi
#elif       Ordinal_ENABLED
            integer(IK) :: rank(GET_SIZE(array))
            integer(IK) :: index(GET_SIZE(array))
            integer(IK) :: i, currentRank
            isSortedPresent = logical(present(isSorted), LK) ! ifort 2021.1 bug: ifort cannot recognize `present`.
            if (isSortedPresent) then
                call setSorted(array, index, isSorted)
            else
                call setSorted(array, index)
            end if
            i = 0_IK
            currentRank = 0_IK
            loopi: do
                i = i + 1_IK
                if (i > GET_SIZE(array)) exit
                currentRank = currentRank + 1_IK
                rank(index(i)) = currentRank
            end do loopi
#else
#error      "Unrecognized Interface."
#endif
        end function

        PURE function isSorted(a, b) result(sorted)
            logical(LK) :: sorted
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in) :: a, b
#elif       SK_ENABLED && D1_ENABLED
            character(*,SKG), intent(in) :: a, b
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in) :: a, b
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in) :: a, b
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in) :: a, b
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in) :: a, b
#elif       PSSK_ENABLED && D1_ENABLED
            type(strc(SKG)) , intent(in) :: a, b
#else
#error      "Unrecognized Interface."
#endif
            sorted = logical(a > b, LK)
        end function

        subroutine report(line, isSorted)
            integer(IK), intent(in) :: line
            procedure(logical(LK)), optional :: isSorted
#if         getRank_ENABLED
            if (present(isSorted)) then
                rank = GETRANK(array, isSorted)
            else
                rank = GETRANK(array)
            endif
#elif       setRank_ENABLED
            if (allocated(rank)) deallocate(rank)
            allocate(rank, mold = rank_ref)
            if (present(isSorted)) then
                call SETRANK(rank, array, isSorted)
            else
                call SETRANK(rank, array)
            endif
#else
#error      "Unrecognized Interface."
#endif
            assertion = assertion .and. logical(all(rank == rank_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "array", array
                write(test%disp%unit,"(*(g0,:,', '))") "rank", rank
                write(test%disp%unit,"(*(g0,:,', '))") "rank_ref", rank_ref
                write(test%disp%unit,"(*(g0,:,', '))") "rank == rank_ref", rank == rank_ref
                write(test%disp%unit,"(*(g0,:,', '))") "present(isSorted)", present(isSorted)
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_": The rank of the input array must be computed correctly.", line)
        end subroutine

#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_EQUAL
#undef  GETRANK
#undef  SETRANK