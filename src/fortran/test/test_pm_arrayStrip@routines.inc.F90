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
!>  This module contains implementations of the tests of the procedures of [pm_arrayStrip](@ref pm_arrayStrip).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     D1_D0_ENABLED
#define ISEQ iseqVec
#endif
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
#if     D0_D0_ENABLED && SK_ENABLED
#define GET_REPEAT(x, count) repeat(x, count)
#define GET_SIZE len
#define ALL
#else
#define GET_REPEAT(x, count) x
#define GET_SIZE size
#endif
        !%%%%%%%%%%%%%%%%%%
#if     getStripped_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_D0_ENABLED
        character(:,SKG), allocatable   :: arrayStripped, arrayStripped_ref, array, pattern
        character(1,SKG), parameter     :: lower = SKG_"a", upper = SKG_"d"
#elif   SK_ENABLED && D1_D0_ENABLED
        character(2,SKG), allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern
        character(2,SKG), parameter     :: lower = SKG_"aa", upper = SKG_"dd"
#elif   IK_ENABLED && D1_D0_ENABLED
        integer(IKG)    , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern
        integer(IKG)    , parameter     :: lower = 0_IKG, upper = 10_IKG
#elif   LK_ENABLED && D1_D0_ENABLED
        logical(LKG)    , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern
        logical(LKG)    , parameter     :: lower = .false._LKG, upper = .true._LKG
#elif   CK_ENABLED && D1_D0_ENABLED
        complex(CKG)    , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern
        complex(CKG)    , parameter     :: lower = (-1._CKG, -1._CKG), upper = (1._CKG, 1._CKG)
#elif   RK_ENABLED && D1_D0_ENABLED
        real(RKG)       , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern
        real(RKG)       , parameter     :: lower = -1._RKG, upper = 1._RKG
#elif   SK_ENABLED && D1_D1_ENABLED
        character(2,SKG), allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern(:)
        character(2,SKG), parameter     :: lower = SKG_"aa", upper = SKG_"dd"
#elif   IK_ENABLED && D1_D1_ENABLED
        integer(IKG)    , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern(:)
        integer(IKG)    , parameter     :: lower = 0_IKG, upper = 10_IKG
#elif   LK_ENABLED && D1_D1_ENABLED
        logical(LKG)    , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern(:)
        logical(LKG)    , parameter     :: lower = .false._LKG, upper = .true._LKG
#elif   CK_ENABLED && D1_D1_ENABLED
        complex(CKG)    , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern(:)
        complex(CKG)    , parameter     :: lower = (-1._CKG, -1._CKG), upper = (1._CKG, 1._CKG)
#elif   RK_ENABLED && D1_D1_ENABLED
        real(RKG)       , allocatable   :: arrayStripped(:), arrayStripped_ref(:), array(:), pattern(:)
        real(RKG)       , parameter     :: lower = -1._RKG, upper = 1._RKG
#else
#error  "Unrecognized interface."
#endif
#if     SB_ENABLED
#define SIDE_TYPE leftRight_type
#elif   SR_ENABLED
#define SIDE_TYPE right_type
#elif   LR_ENABLED
#define SIDE_TYPE left_type
#else
#error  "Unrecognized interface."
#endif
        type(SIDE_TYPE), parameter :: side = SIDE_TYPE()
        assertion = .true._LK
        call runTestsWith()
        call runTestsWith(iseq)

    contains

        subroutine runTestsWith(iseq)

            logical(LK), external, optional :: iseq
            integer(IK) :: i, lenArray

            do i = 1_IK, 300_IK

                call reset()
                lenArray = getUnifRand(0_IK, 100_IK)
                call setResized(array, size = lenArray)
                call setUnifRand(array, GET_REPEAT(lower, lenArray), GET_REPEAT(upper, lenArray))
                ! set pattern.
                if (getUnifRand()) then
#if                 D0_D0_ENABLED || D1_D1_ENABLED
                    block
                        integer(IK) :: lenpattern
                        lenpattern = getUnifRand(0, 2)
                        call setResized(pattern, lenpattern)
                        call setUnifRand(pattern, GET_REPEAT(lower, lenpattern), GET_REPEAT(upper, lenpattern))
                    end block
#elif               D1_D0_ENABLED
                    pattern = getUnifRand(lower, upper)

#else
#error              "Unrecognized interface."
#endif
                else
#if                 D0_D0_ENABLED || D1_D1_ENABLED
                    block
                        integer(IK) :: lindex, rindex
                        if (getUnifRand()) then
                            lindex = 1_IK
                            rindex = getUnifRand(0_IK, min(2_IK, lenArray))
                        else
                            lindex = lenArray - getUnifRand(0_IK, min(2_IK, lenArray)) + 1_IK
                            rindex = lenArray
                        end if
                        pattern = array(lindex : rindex)
                    end block
#elif               D1_D0_ENABLED
                    if (lenArray > 0_IK) then
                        pattern = merge(array(1), array(lenArray), getUnifRand())
                    else
                        pattern = getUnifRand(lower, upper)
                    end if
#else
#error              "Unrecognized interface."
#endif
                end if

                ! strip.

                if (present(iseq)) then
#if                 SB_ENABLED
                    arrayStripped_ref = array(getSIL(array, pattern, iseq) : getSIR(array, pattern, iseq))
#elif               LR_ENABLED
                    arrayStripped_ref = array(getSIL(array, pattern, iseq) : )
#elif               SR_ENABLED
                    arrayStripped_ref = array( : getSIR(array, pattern, iseq))
#else
#error              "Unrecognized interface."
#endif
                    arrayStripped = getStripped(array, pattern, iseq, side)
                else
#if                 SB_ENABLED
                    arrayStripped_ref = array(getSIL(array, pattern) : getSIR(array, pattern))
#elif               LR_ENABLED
                    arrayStripped_ref = array(getSIL(array, pattern) : )
#elif               SR_ENABLED
                    arrayStripped_ref = array( : getSIR(array, pattern))
#endif
                    arrayStripped = getStripped(array, pattern, side)
                end if
                call report(__LINE__, iseq, side)
#if             SB_ENABLED
                arrayStripped = getStripped(array, pattern)
                call report(__LINE__, iseq)
#endif

            end do

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, iseq, side)
            type(SIDE_TYPE), intent(in), optional :: side
            logical(LK) , external, optional :: iseq
            integer, intent(in) :: line
            assertion = assertion .and. logical(ALL(arrayStripped IS_EQUAL arrayStripped_ref), LK) ! fpp
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("arrayStripped_ref")
                call test%disp%show( arrayStripped_ref )
                call test%disp%show("arrayStripped")
                call test%disp%show( arrayStripped )
                call test%disp%show("array")
                call test%disp%show( array )
                call test%disp%show("pattern")
                call test%disp%show( pattern )
                call test%disp%show("present(iseq)")
                call test%disp%show( present(iseq) )
                call test%disp%show("present(side)")
                call test%disp%show( present(side) )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"@getStripped(): The test array must be stripped correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     D0_D0_ENABLED || D1_D0_ENABLED
        function iseq(Segment, pattern) result(equivalent)
#if         SK_ENABLED && D0_D0_ENABLED
            character(*,SKG), intent(in)    :: Segment, pattern
#elif       SK_ENABLED && D1_D0_ENABLED
            character(*,SKG), intent(in)    :: Segment, pattern
#elif       IK_ENABLED && D1_D0_ENABLED
            integer(IKG)    , intent(in)    :: Segment, pattern
#elif       LK_ENABLED && D1_D0_ENABLED
            logical(LKG)    , intent(in)    :: Segment, pattern
#elif       CK_ENABLED && D1_D0_ENABLED
            complex(CKG)    , intent(in)    :: Segment, pattern
#elif       RK_ENABLED && D1_D0_ENABLED
            real(RKG)       , intent(in)    :: Segment, pattern
#endif
            logical(LK) :: equivalent
            equivalent = Segment IS_EQUAL pattern
        end function
#elif   D1_D1_ENABLED
        function iseq(Segment, pattern, lenpattern) result(equivalent)
            logical(LK)             :: equivalent
            integer(IK), intent(in) :: lenpattern
#if         SK_ENABLED
            character(*,SKG), intent(in)    :: Segment(lenpattern), pattern(lenpattern)
#elif       IK_ENABLED
            integer(IKG)    , intent(in)    :: Segment(lenpattern), pattern(lenpattern)
#elif       LK_ENABLED
            logical(LKG)    , intent(in)    :: Segment(lenpattern), pattern(lenpattern)
#elif       CK_ENABLED
            complex(CKG)    , intent(in)    :: Segment(lenpattern), pattern(lenpattern)
#elif       RK_ENABLED
            real(RKG)       , intent(in)    :: Segment(lenpattern), pattern(lenpattern)
#endif
            equivalent = all(Segment IS_EQUAL pattern)
        end function
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(array)) deallocate(array)
            if (allocated(pattern)) deallocate(pattern)
            if (allocated(arrayStripped)) deallocate(arrayStripped)
            if (allocated(arrayStripped_ref)) deallocate(arrayStripped_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getSIL_ENABLED || getSIR_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && D0_D0_ENABLED
        character(:,SKG), allocatable   :: array, pattern
        character(1,SKG), parameter     :: lower = SKG_"a", upper = SKG_"d"
#elif   SK_ENABLED && D1_D0_ENABLED
        character(2,SKG), allocatable   :: array(:), pattern
        character(2,SKG), parameter     :: lower = SKG_"aa", upper = SKG_"dd"
#elif   IK_ENABLED && D1_D0_ENABLED
        integer(IKG)    , allocatable   :: array(:), pattern
        integer(IKG)    , parameter     :: lower = 0_IKG, upper = 10_IKG
#elif   LK_ENABLED && D1_D0_ENABLED
        logical(LKG)    , allocatable   :: array(:), pattern
        logical(LKG)    , parameter     :: lower = .false._LKG, upper = .true._LKG
#elif   CK_ENABLED && D1_D0_ENABLED
        complex(CKG)    , allocatable   :: array(:), pattern
        complex(CKG)    , parameter     :: lower = (-1._CKG, -1._CKG), upper = (1._CKG, 1._CKG)
#elif   RK_ENABLED && D1_D0_ENABLED
        real(RKG)       , allocatable   :: array(:), pattern
        real(RKG)       , parameter     :: lower = -1._RKG, upper = 1._RKG
#elif   SK_ENABLED && D1_D1_ENABLED
        character(2,SKG), allocatable   :: array(:), pattern(:)
        character(2,SKG), parameter     :: lower = SKG_"aa", upper = SKG_"dd"
#elif   IK_ENABLED && D1_D1_ENABLED
        integer(IKG)    , allocatable   :: array(:), pattern(:)
        integer(IKG)    , parameter     :: lower = 0_IKG, upper = 10_IKG
#elif   LK_ENABLED && D1_D1_ENABLED
        logical(LKG)    , allocatable   :: array(:), pattern(:)
        logical(LKG)    , parameter     :: lower = .false._LKG, upper = .true._LKG
#elif   CK_ENABLED && D1_D1_ENABLED
        complex(CKG)    , allocatable   :: array(:), pattern(:)
        complex(CKG)    , parameter     :: lower = (-1._CKG, -1._CKG), upper = (1._CKG, 1._CKG)
#elif   RK_ENABLED && D1_D1_ENABLED
        real(RKG)       , allocatable   :: array(:), pattern(:)
        real(RKG)       , parameter     :: lower = -1._RKG, upper = 1._RKG
#else
#error  "Unrecognized interface."
#endif

#if     getSIL_ENABLED
#define GETSIX getSIL
#elif   getSIR_ENABLED
#define GETSIX getSIR
#endif
#if     D1_D0_ENABLED
#define ISEQ iseqVec
#endif
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
#if     SK_ENABLED && D0_D0_ENABLED
#define GET_SIZE len
#else
#define GET_SIZE size
#endif
#if     getSIL_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@getSIL()"
#elif   getSIR_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@getSIR()"
#endif
        integer(IK) :: index, index_ref

        assertion = .true._LK
        call runTestsWith()
        call runTestsWith(iseq)

    contains

        subroutine runTestsWith(iseq)
            logical(LK), external, optional :: iseq

#if         D1_D0_ENABLED

            integer(IK) :: i
            do i = 1_IK, 200_IK
                call reset()
                array = getUnifRand(lower, upper, s1 = getUnifRand(0_IK, 100_IK))
                pattern = getUnifRand(lower, upper)
                call report(__LINE__, iseq)
            end do

#elif       D0_D0_ENABLED || D1_D1_ENABLED

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_""
            pattern = SKG_""
#elif       SK_ENABLED && D1_D1_ENABLED
            allocate(character(2,SKG) :: array(0), pattern(0))
#elif       IK_ENABLED && D1_D1_ENABLED
            allocate(array(0), pattern(0))
#elif       LK_ENABLED && D1_D1_ENABLED
            allocate(array(0), pattern(0))
#elif       CK_ENABLED && D1_D1_ENABLED
            allocate(array(0), pattern(0))
#elif       RK_ENABLED && D1_D1_ENABLED
            allocate(array(0), pattern(0))
#endif

#if         getSIL_ENABLED
            index_ref = 1_IK
#elif       getSIR_ENABLED
            index_ref = 0_IK
#endif
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_"aaabb"
            pattern = SKG_""
#elif       SK_ENABLED && D1_D1_ENABLED
            array = [character(2,SKG) :: "aa", "aa", "aa", "bb", "bb"]
            pattern = [character(2,SKG) ::]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [integer(IKG) :: 1, 1, 1, 2, 2]
            pattern = [integer(IKG) ::]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [logical(LKG) :: .false., .false., .false., .true., .true.]
            pattern = [logical(LKG) ::]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [complex(CKG) :: 1, 1, 1, 2, 2]
            pattern = [complex(CKG) ::]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [real(RKG) :: 1, 1, 1, 2, 2]
            pattern = [real(RKG) ::]
#endif

#if         getSIL_ENABLED
            index_ref = 1_IK
#elif       getSIR_ENABLED
            index_ref = 5_IK
            call setReversed(array)
#endif
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_"aaabb"
            pattern = SKG_"a"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = [character(2,SKG) :: "aa", "aa", "aa", "bb", "bb"]
            pattern = [character(2,SKG) :: "aa"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [integer(IKG) :: 1, 1, 1, 2, 2]
            pattern = [integer(IKG) :: 1]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [logical(LKG) :: .false., .false., .false., .true., .true.]
            pattern = [logical(LKG) :: .false.]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [complex(CKG) :: 1, 1, 1, 2, 2]
            pattern = [complex(CKG) :: 1]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [real(RKG) :: 1, 1, 1, 2, 2]
            pattern = [real(RKG) :: 1]
#endif

#if         getSIL_ENABLED
            index_ref = 4_IK
#elif       getSIR_ENABLED
            index_ref = 2_IK
            call setReversed(array)
#endif
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_"aaabb"
            pattern = SKG_"aa"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = [character(2,SKG) :: "aa", "aa", "aa", "bb", "bb"]
            pattern = [character(2,SKG) :: "aa", "aa"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [integer(IKG) :: 1, 1, 1, 2, 2]
            pattern = [integer(IKG) :: 1, 1]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [logical(LKG) :: .false., .false., .false., .true., .true.]
            pattern = [logical(LKG) :: .false., .false.]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [complex(CKG) :: 1, 1, 1, 2, 2]
            pattern = [complex(CKG) :: 1, 1]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [real(RKG) :: 1, 1, 1, 2, 2]
            pattern = [real(RKG) :: 1, 1]
#endif

#if         getSIL_ENABLED
            index_ref = 3_IK
#elif       getSIR_ENABLED
            index_ref = 3_IK
            call setReversed(array)
#endif
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_"aaabb"
            pattern = SKG_"b"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = [character(2,SKG) :: "aa", "aa", "aa", "bb", "bb"]
            pattern = [character(2,SKG) :: "bb"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [integer(IKG) :: 1, 1, 1, 2, 2]
            pattern = [integer(IKG) :: 2]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [logical(LKG) :: .false., .false., .false., .true., .true.]
            pattern = [logical(LKG) :: .true.]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [complex(CKG) :: 1, 1, 1, 2, 2]
            pattern = [complex(CKG) :: 2]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [real(RKG) :: 1, 1, 1, 2, 2]
            pattern = [real(RKG) :: 2]
#endif

#if         getSIL_ENABLED
            index_ref = 1_IK
#elif       getSIR_ENABLED
            index_ref = 5_IK
            call setReversed(array)
#endif
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_"aaabb"
            pattern = SKG_"bb"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = [character(2,SKG) :: "aa", "aa", "aa", "bb", "bb"]
            pattern = [character(2,SKG) :: "bb", "bb"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [integer(IKG) :: 1, 1, 1, 2, 2]
            pattern = [integer(IKG) :: 2, 2]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [logical(LKG) :: .false., .false., .false., .true., .true.]
            pattern = [logical(LKG) :: .true., .true.]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [complex(CKG) :: 1, 1, 1, 2, 2]
            pattern = [complex(CKG) :: 2, 2]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [real(RKG) :: 1, 1, 1, 2, 2]
            pattern = [real(RKG) :: 2, 2]
#endif

#if         getSIL_ENABLED
            index_ref = 1_IK
#elif       getSIR_ENABLED
            index_ref = 5_IK
            call setReversed(array)
#endif
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_"aaabb"
#elif       SK_ENABLED && D1_D1_ENABLED
            array = [character(2,SKG) :: "aa", "aa", "aa", "bb", "bb"]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [integer(IKG) :: 1, 1, 1, 2, 2]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [logical(LKG) :: .false., .false., .false., .true., .true.]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [complex(CKG) :: 1, 1, 1, 2, 2]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [real(RKG) :: 1, 1, 1, 2, 2]
#endif

#if         getSIL_ENABLED
            index_ref = 6_IK
#elif       getSIR_ENABLED
            index_ref = 0_IK
            call setReversed(array)
#endif
            pattern = array
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         SK_ENABLED && D0_D0_ENABLED
            array = SKG_"aaabb"
            pattern = array//array
#elif       SK_ENABLED && D1_D1_ENABLED
            array = [character(2,SKG) :: "aa", "aa", "aa", "bb", "bb"]
            pattern = [array, array]
#elif       IK_ENABLED && D1_D1_ENABLED
            array = [integer(IKG) :: 1, 1, 1, 2, 2]
            pattern = [array, array]
#elif       LK_ENABLED && D1_D1_ENABLED
            array = [logical(LKG) :: .false., .false., .false., .true., .true.]
            pattern = [array, array]
#elif       CK_ENABLED && D1_D1_ENABLED
            array = [complex(CKG) :: 1, 1, 1, 2, 2]
            pattern = [array, array]
#elif       RK_ENABLED && D1_D1_ENABLED
            array = [real(RKG) :: 1, 1, 1, 2, 2]
            pattern = [array, array]
#endif

#if         getSIL_ENABLED
            index_ref = 1_IK
#elif       getSIR_ENABLED
            index_ref = 5_IK
            call setReversed(array)
#endif
            call report(__LINE__, iseq)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
#error      "Unrecognized interface."
#endif

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, iseq)
            integer     , intent(in)         :: line
            logical(LK) , external, optional :: iseq

            if (present(iseq)) then

#if             D1_D0_ENABLED
                index_ref = GETSIX(array, [pattern], ISEQ) ! fpp
#endif
                index = GETSIX(array, pattern, iseq)
            else
#if             D1_D0_ENABLED
                index_ref = GETSIX(array, [pattern])
#endif
                index = GETSIX(array, pattern)
            end if
            assertion = assertion .and. logical(index == index_ref, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("index_ref")
                call test%disp%show( index_ref )
                call test%disp%show("index")
                call test%disp%show( index )
                call test%disp%show("array")
                call test%disp%show( array )
                call test%disp%show("pattern")
                call test%disp%show( pattern )
                call test%disp%show("present(iseq)")
                call test%disp%show( present(iseq) )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, PROCEDURE_NAME//SK_": The `index` of the stripped array must be computed correctly.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     D0_D0_ENABLED
        function iseq(Segment, pattern) result(equivalent)
            character(*,SKG), intent(in) :: pattern, Segment
            logical(LK) :: equivalent
            equivalent = Segment == pattern
        end function
#elif   D1_D0_ENABLED || D1_D1_ENABLED
#if     D1_D0_ENABLED
        function iseq(Segment, pattern) result(equivalent)
#if         SK_ENABLED
            character(*,SKG), intent(in)    :: Segment, pattern
#elif       IK_ENABLED
            integer(IKG)    , intent(in)    :: Segment, pattern
#elif       LK_ENABLED
            logical(LKG)    , intent(in)    :: Segment, pattern
#elif       CK_ENABLED
            complex(CKG)    , intent(in)    :: Segment, pattern
#elif       RK_ENABLED
            real(RKG)       , intent(in)    :: Segment, pattern
#endif
            logical(LK) :: equivalent
            equivalent = Segment IS_EQUAL pattern
        end function
#endif
        function ISEQ(Segment, pattern, lenPattern) result(equivalent) ! fpp
            logical(LK)             :: equivalent
            integer(IK), intent(in) :: lenPattern
#if         SK_ENABLED
            character(*,SKG), intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       IK_ENABLED
            integer(IKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       LK_ENABLED
            logical(LKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       CK_ENABLED
            complex(CKG)    , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#elif       RK_ENABLED
            real(RKG)       , intent(in)    :: Segment(lenPattern), pattern(lenPattern)
#endif
            equivalent = all(Segment IS_EQUAL pattern)
        end function
#else
#error  "Unrecognized interface."
#endif
        subroutine reset()
            if (allocated(array)) deallocate(array)
            if (allocated(pattern)) deallocate(pattern)
        end subroutine
#undef  IS_EQUAL
#undef  GET_SIZE
#undef  ISEQ
#undef  ALL

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#undef  GET_REPEAT
#undef  SIDE_TYPE
#undef  IS_EQUAL
#undef  GET_SIZE
#undef  GETSIX
#undef  ISEQ
#undef  ALL