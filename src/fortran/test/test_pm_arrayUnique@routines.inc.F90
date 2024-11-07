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
!>  This module contains implementations of the tests of the procedures under the generic interfaces
!>  [getUnique](@ref pm_arrayUnique::getUnique),
!>  [setUnique](@ref pm_arrayUnique::setUnique).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the comparison operator.
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#elif   SK_ENABLED || IK_ENABLED || LK_ENABLED || CK_ENABLED || RK_ENABLED
#define IS_EQUAL ==
#else
#error  "Unrecognized interface."
#endif
        ! Define the slicing rule.
#if     SK_ENABLED && D0_ENABLED
#define GET_INDEX(i) i:i
#define GET_SIZE len
#elif   D1_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%
#if     isUnique_ENABLED
        !%%%%%%%%%%%%%%%

        character(*, SK), parameter :: PROCEDURE_NAME = "@isUnique()"
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG), allocatable :: array
        character(1,SKG), parameter :: lb = "a", ub = "i"
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: array
        character(2,SKG), parameter :: lb = "aa", ub = "ii"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: array
        integer(IKG)    , parameter :: lb = 0, ub = 9
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: array
        logical(LKG)    , parameter :: lb = .false., ub = .true.
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: array
        complex(CKG)    , parameter :: lb = (0., -9.), ub = (+9., 0.)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: array
        real(RKG)       , parameter :: lb = 0., ub = 9.
#else
#error  "Unrecognized interface."
#endif
        type(display_type) :: disp
        logical(LK), allocatable :: unique(:)
        integer(IK) :: lenArray, itry, iell, jell, repetition
        assertion = .true._LK
        do itry = 1, 100

#if         SK_ENABLED && D0_ENABLED
            iell = getUnifRand(0_IK, 9_IK)
            array = getUnifRand(repeat(lb, iell), repeat(ub, iell))
#else
            array = getUnifRand(lb, ub, getUnifRand(0_IK, 9_IK))
#endif
            call report(__LINE__, iseq)
            call report(__LINE__)

        end do

    contains

        subroutine report(line, iseq)
            integer     , intent(in)            :: line
            logical(LK) , external  , optional  :: iseq
            if (present(iseq)) then
                unique = isUnique(array, iseq)
            else
                unique = isUnique(array)
            end if
            lenArray = GET_SIZE(array, kind = IK)
            assertion = assertion .and. size(unique, 1, IK) == lenArray
            call test%assert(assertion, PROCEDURE_NAME//SK_": The length of the output `logical` array must match that of the input sequence.", line)
            if (0_IK < lenArray) then
                do iell = 1, lenArray
                    if (present(iseq)) then
                        repetition = 0
                        do jell = 1, lenArray
                            if (iseq(array(GET_INDEX(iell)), array(GET_INDEX(jell)))) repetition = repetition + 1
                        end do
                        assertion = assertion .and. (unique(iell) .eqv. 1_IK == repetition)
                    else
#if                     SK_ENABLED && D0_ENABLED
                        assertion = assertion .and. (unique(iell) .eqv. 1_IK == count(array(iell:iell) == getCharVec(array)))
#else
                        assertion = assertion .and. (unique(iell) .eqv. 1_IK == count(array(iell) IS_EQUAL array))
#endif
                    end if
                end do
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    call disp%skip
                    call disp%show("array")
                    call disp%show( array )
                    call disp%show("unique")
                    call disp%show( unique )
                    call disp%skip
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, PROCEDURE_NAME//SK_": The output argument `unique` must be correctly set.", line)
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%
#elif   isUniqueAll_ENABLED
        !%%%%%%%%%%%%%%%%%%

        character(*, SK), parameter :: PROCEDURE_NAME = "@isUniqueAll()"
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG), allocatable :: array
        character(1,SKG), parameter :: lb = "a", ub = "i"
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: array
        character(2,SKG), parameter :: lb = "aa", ub = "ii"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: array
        integer(IKG)    , parameter :: lb = 0, ub = 9
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: array
        logical(LKG)    , parameter :: lb = .false., ub = .true.
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: array
        complex(CKG)    , parameter :: lb = (0., -9.), ub = (+9., 0.)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: array
        real(RKG)       , parameter :: lb = 0., ub = 9.
#else
#error  "Unrecognized interface."
#endif
        type(display_type) :: disp
        logical(LK) :: uniqueAll, allUnique
        integer(IK) :: lenArray, itry
        assertion = .true._LK
        do itry = 1, 100

#if         SK_ENABLED && D0_ENABLED
            lenArray = getUnifRand(0_IK, 9_IK)
            array = getUnifRand(repeat(lb, lenArray), repeat(ub, lenArray))
#else
            array = getUnifRand(lb, ub, getUnifRand(0_IK, 9_IK))
#endif
            call report(__LINE__, iseq)
            call report(__LINE__)

        end do

    contains

        subroutine report(line, iseq)
            integer     , intent(in)            :: line
            logical(LK) , external  , optional  :: iseq
            if (present(iseq)) then
                uniqueAll = isUniqueAll(array, iseq)
                allUnique = all(isUnique(array, iseq))
            else
                uniqueAll = isUniqueAll(array)
                allUnique = all(isUnique(array))
            end if
            lenArray = GET_SIZE(array, kind = IK)
            assertion = assertion .and. (0_IK < lenArray .or. uniqueAll)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty sequence is all-unique elements.", line)
            if (0_IK < lenArray) then
                assertion = assertion .and. (uniqueAll .eqv. allUnique)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    call disp%skip
                    call disp%show("array")
                    call disp%show( array )
                    call disp%show("uniqueAll")
                    call disp%show( uniqueAll )
                    call disp%skip
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, PROCEDURE_NAME//SK_": The output argument `unique` must be correctly set.", line)
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%
#elif   isUniqueAny_ENABLED
        !%%%%%%%%%%%%%%%%%%

        character(*, SK), parameter :: PROCEDURE_NAME = "@isUniqueAny()"
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG), allocatable :: array
        character(1,SKG), parameter :: lb = "a", ub = "i"
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), dimension(:), allocatable :: array
        character(2,SKG), parameter :: lb = "aa", ub = "ii"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: array
        integer(IKG)    , parameter :: lb = 0, ub = 9
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: array
        logical(LKG)    , parameter :: lb = .false., ub = .true.
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: array
        complex(CKG)    , parameter :: lb = (0., -9.), ub = (+9., 0.)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: array
        real(RKG)       , parameter :: lb = 0., ub = 9.
#else
#error  "Unrecognized interface."
#endif
        type(display_type) :: disp
        logical(LK) :: uniqueAny, anyUnique
        integer(IK) :: lenArray, itry
        assertion = .true._LK
        do itry = 1, 100

#if         SK_ENABLED && D0_ENABLED
            lenArray = getUnifRand(0_IK, 9_IK)
            array = getUnifRand(repeat(lb, lenArray), repeat(ub, lenArray))
#else
            array = getUnifRand(lb, ub, getUnifRand(0_IK, 9_IK))
#endif
            call report(__LINE__, iseq)
            call report(__LINE__)

        end do

    contains

        subroutine report(line, iseq)
            integer     , intent(in)            :: line
            logical(LK) , external  , optional  :: iseq
            if (present(iseq)) then
                uniqueAny = isUniqueAny(array, iseq)
                anyUnique = any(isUnique(array, iseq))
            else
                uniqueAny = isUniqueAny(array)
                anyUnique = any(isUnique(array))
            end if
            lenArray = GET_SIZE(array, kind = IK)
            assertion = assertion .and. (0_IK < lenArray .or. .not. uniqueAny)
            call test%assert(assertion, PROCEDURE_NAME//SK_": An empty sequence has non-unique elements.", line)
            if (0_IK < lenArray) then
                assertion = assertion .and. (uniqueAny .eqv. anyUnique)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    call disp%skip
                    call disp%show("array")
                    call disp%show( array )
                    call disp%show("uniqueAny")
                    call disp%show( uniqueAny )
                    call disp%skip
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, PROCEDURE_NAME//SK_": The output argument `unique` must be correctly set.", line)
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getUnique_ENABLED || setUnique_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getUnique_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@getUnique()"
#elif   setUnique_ENABLED
        integer(IK) , allocatable :: Count(:), Count_ref(:)
        type(cvi_type), allocatable :: index(:), Index_ref(:)
        character(*, SK), parameter :: PROCEDURE_NAME = "@setUnique()"
#endif
        integer(IK) :: lenUnique

#if     SK_ENABLED && D0_ENABLED
#define ALL
        character(:,SKG), allocatable :: array, unique, unique_ref
#elif   SK_ENABLED && D1_ENABLED && getUnique_ENABLED
        character(:,SKG), dimension(:), allocatable :: array, unique, unique_ref
#elif   SK_ENABLED && D1_ENABLED && setUnique_ENABLED
        character(2,SKG), dimension(:), allocatable :: array, unique, unique_ref
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)    , dimension(:), allocatable :: array, unique, unique_ref
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)    , dimension(:), allocatable :: array, unique, unique_ref
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)    , dimension(:), allocatable :: array, unique, unique_ref
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)       , dimension(:), allocatable :: array, unique, unique_ref
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)
#if     setUnique_ENABLED
        call runTestsWith(order = +0_IK)
        call runTestsWith(order = +1_IK)
        call runTestsWith(order = -1_IK)
        call runTestsWith(index = index)
        call runTestsWith(index = index, order =  0_IK)
        call runTestsWith(index = index, order =  1_IK)
        call runTestsWith(index = index, order = -1_IK)
        call runTestsWith(iseq = iseq, order =  0_IK)
        call runTestsWith(iseq = iseq, order =  1_IK)
        call runTestsWith(iseq = iseq, order = -1_IK)
        call runTestsWith(iseq = iseq, index = index)
        call runTestsWith(iseq = iseq, index = index, order =  0_IK)
        call runTestsWith(iseq = iseq, index = index, order =  1_IK)
        call runTestsWith(iseq = iseq, index = index, order = -1_IK)
        call runTestsWith(fixed = .true., order = +0_IK)
        call runTestsWith(fixed = .true., order = +1_IK)
        call runTestsWith(fixed = .true., order = -1_IK)
        call runTestsWith(fixed = .true., index = index)
        call runTestsWith(fixed = .true., index = index, order =  0_IK)
        call runTestsWith(fixed = .true., index = index, order =  1_IK)
        call runTestsWith(fixed = .true., index = index, order = -1_IK)
        call runTestsWith(fixed = .true., iseq = iseq, order =  0_IK)
        call runTestsWith(fixed = .true., iseq = iseq, order =  1_IK)
        call runTestsWith(fixed = .true., iseq = iseq, order = -1_IK)
        call runTestsWith(fixed = .true., iseq = iseq, index = index)
        call runTestsWith(fixed = .true., iseq = iseq, index = index, order =  0_IK)
        call runTestsWith(fixed = .true., iseq = iseq, index = index, order =  1_IK)
        call runTestsWith(fixed = .true., iseq = iseq, index = index, order = -1_IK)
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith ( iseq & ! LCOV_EXCL_LINE
#if         setUnique_ENABLED
                                , fixed & ! LCOV_EXCL_LINE
                                , index & ! LCOV_EXCL_LINE
                                , order & ! LCOV_EXCL_LINE
                                )
            logical                             , intent(in)   , optional   :: fixed
            type(cvi_type), allocatable   , intent(inout), optional   :: index(:)
            integer(IK)                         , intent(in)   , optional   :: order
            integer(IK)         , allocatable                               :: RemapIndex(:)
            integer(IK) :: order_def
#else
                                )
#endif
            logical(LK), external, optional :: iseq
#if         setUnique_ENABLED
            order_def = 0_IK
            if (present(order)) order_def = order
            if (allocated(Count_ref)) deallocate(Count_ref)
            if (allocated(Index_ref)) deallocate(Index_ref)
#endif
            if (allocated(array)) deallocate(array)
            if (allocated(unique)) deallocate(unique)
            if (allocated(unique_ref)) deallocate(unique_ref)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            array = ""
            unique_ref = ""
#elif       SK_ENABLED && D1_ENABLED && getUnique_ENABLED
            allocate(character(2,SKG) :: array(0), unique_ref(0))
#elif       SK_ENABLED && D1_ENABLED && setUnique_ENABLED
            allocate(array(0), unique_ref(0))
#elif       IK_ENABLED && D1_ENABLED
            allocate(array(0), unique_ref(0))
#elif       LK_ENABLED && D1_ENABLED
            allocate(array(0), unique_ref(0))
#elif       CK_ENABLED && D1_ENABLED
            allocate(array(0), unique_ref(0))
#elif       RK_ENABLED && D1_ENABLED
            allocate(array(0), unique_ref(0))
#endif

#if         getUnique_ENABLED
            call report(__LINE__, iseq)
#elif       setUnique_ENABLED
            allocate(Count_ref(0), Index_ref(0))
            call report(__LINE__, iseq, fixed, index, order)
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            array = SKG_" "
            unique_ref = SKG_" "
#elif       SK_ENABLED && D1_ENABLED
            array = [character(2,SKG) :: " "]
            unique_ref = [character(2,SKG) :: " "]
#elif       IK_ENABLED && D1_ENABLED
            array = [1_IKG]
            unique_ref = [1_IKG]
#elif       LK_ENABLED && D1_ENABLED
            array = [logical(LKG) :: .false.]
            unique_ref = [logical(LKG) :: .false.]
#elif       CK_ENABLED && D1_ENABLED
            array = [(+1._CKG, -1._CKG)]
            unique_ref = [(+1._CKG, -1._CKG)]
#elif       RK_ENABLED && D1_ENABLED
            array = [1._RKG]
            unique_ref = [1._RKG]
#endif

#if         getUnique_ENABLED
            call report(__LINE__, iseq)
#elif       setUnique_ENABLED
            Count_ref = [1_IK]
            deallocate(Index_ref)
            allocate(Index_ref(1))
            Index_ref(1)%val = [1_IK]
            call report(__LINE__, iseq, fixed, index, order)
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            array = SKG_"  "
#elif       SK_ENABLED && D1_ENABLED
            array = [character(2,SKG) :: " ", " "]
#elif       IK_ENABLED && D1_ENABLED
            array = [1_IKG, 1_IKG]
#elif       LK_ENABLED && D1_ENABLED
            array = [logical(LKG) :: .false., .false.]
#elif       CK_ENABLED && D1_ENABLED
            array = [complex(CKG) :: (+1._CKG, -1._CKG), (+1._CKG, -1._CKG)]
#elif       RK_ENABLED && D1_ENABLED
            array = [real(RKG) :: 1._RKG, 1._RKG]
#endif

            unique_ref = array(GET_INDEX(1))
#if         getUnique_ENABLED
            call report(__LINE__, iseq)
#elif       setUnique_ENABLED
            Count_ref = [2_IK]
            deallocate(Index_ref)
            allocate(Index_ref(1))
            Index_ref(1)%val = [1_IK, 2_IK]
            call report(__LINE__, iseq, fixed, index, order)
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            array = SKG_"ABaB "
            unique_ref = SKG_"ABa "
#elif       SK_ENABLED && D1_ENABLED
            array = [character(2,SKG) :: "A", "B", "a", "B", " "]
            unique_ref = [character(2,SKG) :: "A", "B", "a", " "]
#elif       IK_ENABLED && D1_ENABLED
            array = [1_IKG, 2_IKG, 0_IKG, 2_IKG, -1_IKG]
            unique_ref = [1_IKG, 2_IKG, 0_IKG, -1_IKG]
#elif       LK_ENABLED && D1_ENABLED
            array = [logical(LKG) :: .false., .false., .true., .false.]
            unique_ref = [logical(LKG) :: .false., .true.]
#elif       CK_ENABLED && D1_ENABLED
            array = [(+1._CKG, -1._CKG), (+2._CKG, -2._CKG), (+0._CKG, 0._CKG), (+2._CKG, -2._CKG), (-1._CKG, +1._CKG)]
            unique_ref = [(+1._CKG, -1._CKG), (+2._CKG, -2._CKG), (+0._CKG, 0._CKG), (-1._CKG, +1._CKG)]
#elif       RK_ENABLED && D1_ENABLED
            array = [1._RKG, 2._RKG, 0._RKG, 2._RKG, -1._RKG]
            unique_ref = [1._RKG, 2._RKG, 0._RKG, -1._RKG]
#endif

#if         getUnique_ENABLED
            call report(__LINE__, iseq)
#elif       setUnique_ENABLED
            deallocate(Index_ref)
#if         LK_ENABLED
            allocate(Index_ref(2))
            Count_ref = [3_IK, 1_IK]
            Index_ref(1)%val = [integer(IK) :: 1, 2, 4]
            Index_ref(2)%val = [integer(IK) :: 3]
            if (order_def == 0_IK) then
                RemapIndex = [integer(IK) :: 1, 2]
            elseif (order_def > 0_IK) then
                RemapIndex = [integer(IK) :: 2, 1]
            elseif (order_def < 0_IK) then
                RemapIndex = [integer(IK) :: 1, 2]
            end if
#else
            allocate(Index_ref(4))
            Count_ref = [1_IK, 2_IK, 1_IK, 1_IK]
            Index_ref(1)%val = [1_IK]
            Index_ref(2)%val = [2_IK, 4_IK]
            Index_ref(3)%val = [3_IK]
            Index_ref(4)%val = [5_IK]
            if (order_def == 0_IK) then
                RemapIndex = [1_IK, 2_IK, 3_IK, 4_IK]
            elseif (order_def > 0_IK) then
                RemapIndex = [1_IK, 3_IK, 4_IK, 2_IK]
            elseif (order_def < 0_IK) then
                RemapIndex = [2_IK, 4_IK, 3_IK, 1_IK]
            end if
#endif
            Index_ref = Index_ref(RemapIndex)
            call setRemapped(Count_ref, RemapIndex)
            call setRemapped(unique_ref, RemapIndex)
            call report(__LINE__, iseq, fixed, index, order)
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D0_ENABLED
            array = "ABaX "
#elif       SK_ENABLED && D1_ENABLED
            array = ["A", "B", "a", "X", " "]
#elif       IK_ENABLED && D1_ENABLED
            array = [1_IKG, 2_IKG, 0_IKG, -1_IKG, -2_IKG]
#elif       LK_ENABLED && D1_ENABLED
            array = [logical(LKG) :: .true., .false.]
#elif       CK_ENABLED && D1_ENABLED
            array = [(+1._CKG, -1._CKG), (+0._CKG, 0._CKG), (+2._CKG, -2._CKG), (-1._CKG, +1._CKG), (-2._CKG, +2._CKG)]
#elif       RK_ENABLED && D1_ENABLED
            array = [1._RKG, 0._RKG, 2._RKG, -1._RKG, -2._RKG]
#endif

            unique_ref = array
#if         getUnique_ENABLED
            call report(__LINE__, iseq)
#elif       setUnique_ENABLED && LK_ENABLED
            deallocate(Index_ref)
            allocate(Index_ref(2))
            Count_ref = [integer(IK) :: 1, 1]
            Index_ref(1)%val = [integer(IK) :: 1]
            Index_ref(2)%val = [integer(IK) :: 2]
            if (order_def == 0_IK) then
                RemapIndex = [integer(IK) :: 1, 2]
            elseif (order_def > 0_IK) then
                RemapIndex = [integer(IK) :: 1, 2]
            elseif (order_def < 0_IK) then
                RemapIndex = [integer(IK) :: 2, 1]
            end if
            Index_ref = Index_ref(RemapIndex)
            call setRemapped(Count_ref, RemapIndex)
            call setRemapped(unique_ref, RemapIndex)
            call report(__LINE__, iseq, fixed, index, order)
#elif       setUnique_ENABLED
            deallocate(Index_ref)
            allocate(Index_ref(5))
            Count_ref = [1_IK, 1_IK, 1_IK, 1_IK, 1_IK]
            Index_ref(1)%val = [1_IK]
            Index_ref(2)%val = [2_IK]
            Index_ref(3)%val = [3_IK]
            Index_ref(4)%val = [4_IK]
            Index_ref(5)%val = [5_IK]
            if (order_def == 0_IK) then
                RemapIndex = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK]
            elseif (order_def > 0_IK) then
                RemapIndex = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK]
            elseif (order_def < 0_IK) then
                RemapIndex = [5_IK, 4_IK, 3_IK, 2_IK, 1_IK]
            end if
            Index_ref = Index_ref(RemapIndex)
            call setRemapped(Count_ref, RemapIndex)
            call setRemapped(unique_ref, RemapIndex)
            call report(__LINE__, iseq, fixed, index, order)
#endif

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         SK_ENABLED && D1_ENABLED
            array = ["", "", "", ""]
            unique_ref = [""]
#if         getUnique_ENABLED
            call report(__LINE__)
#elif       setUnique_ENABLED
            Count_ref = [4_IK]
            deallocate(Index_ref)
            allocate(Index_ref(5))
            Index_ref(1)%val = [1_IK, 2_IK, 3_IK, 4_IK]
            call report(__LINE__, iseq, fixed, index, order)
#endif
#endif

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getUnique_ENABLED
        subroutine report(line, iseq)
            integer     , intent(in)            :: line
            logical(LK) , external  , optional  :: iseq
            if (present(iseq)) then
                unique = getUnique(array, iseq)
            else
                unique = getUnique(array)
            end if
            lenUnique = GET_SIZE(unique, kind = IK)
            assertion = assertion .and. ALL(unique(1:lenUnique) IS_EQUAL unique_ref)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "array      ", array
                write(test%disp%unit,"(*(g0,:,', '))") "unique     ", unique(1:lenUnique)
                write(test%disp%unit,"(*(g0,:,', '))") "unique_ref ", unique_ref
                write(test%disp%unit,"(*(g0,:,', '))") "lenUnique  ", lenUnique
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, PROCEDURE_NAME//SK_": The output argument `unique` must be correctly set.", line)
        end subroutine
#elif   setUnique_ENABLED
        subroutine report   ( line & ! LCOV_EXCL_LINE
                            , iseq & ! LCOV_EXCL_LINE
                            , fixed & ! LCOV_EXCL_LINE
                            , index & ! LCOV_EXCL_LINE
                            , order & ! LCOV_EXCL_LINE
                            )
            use pm_option, only: getOption
            integer             , intent(in)                                :: line
            logical(LK)         , external      , optional                  :: iseq
            logical(LK)         , intent(in)    , optional                  :: fixed
            type(cvi_type), intent(inout) , optional , allocatable    :: index(:)
            integer(IK)         , intent(in)    , optional                  :: order
            integer(IK) :: i
            if (getOption(.false._LK, fixed)) then ! Test the contiguous array interfaces.
                if (allocated(Count)) deallocate(Count)
                if (allocated(unique)) deallocate(unique)
                lenUnique = GET_SIZE(array, kind = IK)
                allocate(unique, mold = array)
                allocate(Count(lenUnique))
                if (present(index)) then
                    if (allocated(index)) deallocate(index)
                    allocate(index(lenUnique))
                    if (present(iseq)) then
                        call setUnique(array, unique, lenUnique, Count, iseq = iseq, index = index, order = order)
                    else
                        call setUnique(array, unique, lenUnique, Count, index = index, order = order)
                    end if
                else
                    if (present(iseq)) then
                        call setUnique(array, unique, lenUnique, Count, iseq = iseq, order = order)
                    else
                        call setUnique(array, unique, lenUnique, Count, order = order)
                    end if
                end if
            else
                if (present(index)) then
                    if (present(iseq)) then
                        call setUnique(array, unique, Count, iseq = iseq, index = index, order = order)
                    else
                        call setUnique(array, unique, Count, index = index, order = order)
                    end if
                else
                    if (present(iseq)) then
                        call setUnique(array, unique, Count, iseq = iseq, order = order)
                    else
                        call setUnique(array, unique, Count, order = order)
                    end if
                end if
                lenUnique = GET_SIZE(unique, kind = IK)
            end if

            assertion = assertion .and. lenUnique == GET_SIZE(unique_ref, kind = IK)
            call outputSpec(lenUnique, iseq, fixed, index, order)
            call test%assert(assertion, PROCEDURE_NAME//SK_": The size of the output argument `unique` must be correctly set.", line)

            assertion = assertion .and. all(Count(1:lenUnique) == Count_ref)
            call outputSpec(lenUnique, iseq, fixed, index, order)
            call test%assert(assertion, PROCEDURE_NAME//SK_": The output argument `Count` must be correctly set.", line)

            if (present(index)) then
                do i = 1, lenUnique
                    assertion = assertion .and. all(index(i)%val == Index_ref(i)%val)
                    call outputSpec(lenUnique, iseq, fixed, index, order)
                    call test%assert(assertion, PROCEDURE_NAME//SK_": The output argument `index` must be correctly set.", line)
                end do
            end if

            assertion = assertion .and. ALL(unique(1:lenUnique) IS_EQUAL unique_ref)
            call outputSpec(lenUnique, iseq, fixed, index, order)
            call test%assert(assertion, PROCEDURE_NAME//SK_": The output argument `unique` must be correctly set.", line)

        end subroutine

        subroutine outputSpec(lenUnique, iseq, fixed, index, order)
            integer(IK)         , intent(in)                            :: lenUnique
            logical(LK)         , external  , optional                  :: iseq
            logical(LK)         , intent(in), optional                  :: fixed
            type(cvi_type), intent(in), optional , allocatable    :: index(:)
            integer(IK)         , intent(in), optional                  :: order
            integer(IK) :: i
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "array          ", array
                write(test%disp%unit,"(*(g0,:,', '))") "lenArray       ", GET_SIZE(array, kind = IK)
                write(test%disp%unit,"(*(g0,:,', '))") "unique         ", unique(1:lenUnique)
                write(test%disp%unit,"(*(g0,:,', '))") "unique_ref     ", unique_ref
                write(test%disp%unit,"(*(g0,:,', '))") "lenUnique      ", lenUnique
                write(test%disp%unit,"(*(g0,:,', '))") "Count          ", Count(1:lenUnique)
                write(test%disp%unit,"(*(g0,:,', '))") "Count_ref      ", Count_ref
                write(test%disp%unit,"(*(g0,:,', '))") "present(iseq)  ", present(iseq)
                write(test%disp%unit,"(*(g0,:,', '))") "present(fixed) ", present(fixed)
                write(test%disp%unit,"(*(g0,:,', '))")
                if (present(fixed)) then
                write(test%disp%unit,"(*(g0,:,', '))") "fixed          ", fixed
                end if
                write(test%disp%unit,"(*(g0,:,', '))") "present(order) ", present(order)
                if (present(order)) then
                write(test%disp%unit,"(*(g0,:,', '))") "order          ", order
                end if
                write(test%disp%unit,"(*(g0,:,', '))") "present(index) ", present(index)
                if (present(index)) then
                do i = 1, lenUnique
                write(test%disp%unit,"(*(g0,:,', '))") "index          ", index(i)%val
                write(test%disp%unit,"(*(g0,:,', '))") "Index_ref      ", Index_ref(i)%val
                write(test%disp%unit,"(*(g0,:,', '))")
                end do
                end if
                ! LCOV_EXCL_STOP
            end if
        end subroutine
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

        pure function iseq(element1, element2) result(equivalent)
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG), intent(in) :: element1, element2
#elif       SK_ENABLED && D1_ENABLED
            character(*,SKG), intent(in) :: element1, element2
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)    , intent(in) :: element1, element2
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)    , intent(in) :: element1, element2
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)    , intent(in) :: element1, element2
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)       , intent(in) :: element1, element2
#endif
            logical(LK) :: equivalent
            equivalent = element1 IS_EQUAL element2
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GET_INDEX
#undef GET_SIZE
#undef IS_EQUAL
#undef ALL