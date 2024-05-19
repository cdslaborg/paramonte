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
!>  This module contains implementations of the tests of the procedures under the generic interface [setSorted](@ref pm_arraySelect::setSorted).
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setSelected_D1_PSSK_ENABLED || getSelected_D1_PSSK_ENABLED
#define GET_COMP(X)X%val
#else
#define GET_COMP(X)X
#endif

#if     setSelected_D0_SK_ENABLED || getSelected_D0_SK_ENABLED
#define GET_INDEX(i) i:i
#else
#define GET_INDEX(i) i
#endif

#if     setSelected_D1_CK_ENABLED || getSelected_D1_CK_ENABLED
        use pm_complexCompareLex, only: operator(>), operator(<)
#endif

#if     setSelected_D1_LK_ENABLED || getSelected_D1_LK_ENABLED
        use pm_logicalCompare, only: operator(>), operator(<)
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
        integer(IK) , parameter     :: NDATA = 1000_IK
#if     setSelected_D0_SK_ENABLED   || getSelected_D0_SK_ENABLED
        character(NDATA,SKG)        :: dataUnsorted, DataUnsorted_ref
        character(1,SKG)            :: selection
        call setUnifRand(DataUnsorted_ref)
#elif   setSelected_D1_SK_ENABLED   || getSelected_D1_SK_ENABLED
        character(2,SKG)            :: dataUnsorted(NDATA), DataUnsorted_ref(NDATA), selection
        call setUnifRand(DataUnsorted_ref)
#elif   setSelected_D1_IK_ENABLED   || getSelected_D1_IK_ENABLED
        integer(IKG)                :: dataUnsorted(NDATA), DataUnsorted_ref(NDATA), selection
        call setUnifRand(DataUnsorted_ref, 0_IKG, huge(DataUnsorted_ref) - 1_IKG)
#elif   setSelected_D1_LK_ENABLED   || getSelected_D1_LK_ENABLED
        logical(LKG)                :: dataUnsorted(NDATA), DataUnsorted_ref(NDATA), selection
        call setUnifRand(DataUnsorted_ref)
#elif   setSelected_D1_CK_ENABLED   || getSelected_D1_CK_ENABLED
        complex(CKG)                :: dataUnsorted(NDATA), DataUnsorted_ref(NDATA), selection
        call setUnifRand(DataUnsorted_ref)
#elif   setSelected_D1_RK_ENABLED   || getSelected_D1_RK_ENABLED
        real(RKG)                   :: dataUnsorted(NDATA), DataUnsorted_ref(NDATA), selection
        call setUnifRand(DataUnsorted_ref)
#elif   setSelected_D1_PSSK_ENABLED || getSelected_D1_PSSK_ENABLED
        type(css_pdt(SKG))          :: dataUnsorted(NDATA), DataUnsorted_ref(NDATA), selection
        integer(IK)                 :: stringSize
        integer(IK)                 :: i
        do i = 1_IK, NDATA
            call setUnifRand(stringSize, 1_IK, 100_IK)
            allocate(character(stringSize,SKG) :: DataUnsorted_ref(i)%val)
            call setUnifRand(DataUnsorted_ref(i)%val)
        end do
#else
#error  "Unrecognized Interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call runWith()
        call runWith(isAscending_local)
        call runWith(isDescending_local)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runWith(isSorted_local)

            logical(LK), external, optional :: isSorted_local

            dataUnsorted = DataUnsorted_ref; call runTestsWith(rank = 1_IK)
            dataUnsorted = DataUnsorted_ref; call runTestsWith(rank = NDATA / 2)
            dataUnsorted = DataUnsorted_ref; call runTestsWith(rank = NDATA)

            !> \warning The following tests are ordered.
            dataUnsorted = DataUnsorted_ref
            call runTestsWith(isSorted_local = isSorted_local, rank = 1_IK, lb = 1_IK)
            call runTestsWith(isSorted_local = isSorted_local, rank = NDATA / 4_IK - 1_IK)
            call runTestsWith(isSorted_local = isSorted_local, rank = NDATA / 2_IK, lb = NDATA / 4_IK - 1_IK)
            call runTestsWith(isSorted_local = isSorted_local, rank = NDATA, lb = NDATA / 2_IK + 1_IK)

            !> \warning The following tests are ordered.
            dataUnsorted = DataUnsorted_ref
            call runTestsWith(isSorted_local = isSorted_local, rank = NDATA, ub = NDATA)
            call runTestsWith(isSorted_local = isSorted_local, rank = 1_IK, ub = 3 * NDATA / 4_IK + 1_IK)
            call runTestsWith(isSorted_local = isSorted_local, rank = NDATA / 2_IK, ub = 3_IK * NDATA / 4_IK - 1_IK)

            !> \warning The following tests are ordered.
            dataUnsorted = DataUnsorted_ref
            call runTestsWith(isSorted_local = isSorted_local, rank = 1_IK, lb = 1_IK, ub = NDATA)
            call runTestsWith(isSorted_local = isSorted_local, rank = NDATA / 2_IK, lb = NDATA / 4_IK - 1_IK, ub = 3_IK * NDATA / 4_IK - 1_IK)
            call runTestsWith(isSorted_local = isSorted_local, rank = NDATA, lb = NDATA / 4_IK + 1_IK, ub = NDATA)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(isSorted_local, rank, lb, ub)

            logical(LK), external, optional :: isSorted_local
            integer(IK), optional           :: rank, lb, ub

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if         setSelected_D0_SK_ENABLED || setSelected_D1_SK_ENABLED || setSelected_D1_IK_ENABLED || setSelected_D1_LK_ENABLED || setSelected_D1_CK_ENABLED || setSelected_D1_RK_ENABLED || setSelected_D1_PSSK_ENABLED
            if (present(isSorted_local)) then
                call setSelected(selection, dataUnsorted, rank, isSorted_local, lb, ub)
            else
                call setSelected(selection, dataUnsorted, rank, lb, ub)
            end if
#elif       getSelected_D0_SK_ENABLED || getSelected_D1_SK_ENABLED || getSelected_D1_IK_ENABLED || getSelected_D1_LK_ENABLED || getSelected_D1_CK_ENABLED || getSelected_D1_RK_ENABLED || getSelected_D1_PSSK_ENABLED
            if (present(isSorted_local)) then
                selection = getSelected(dataUnsorted, rank, isSorted_local, lb, ub)
            else
                selection = getSelected(dataUnsorted, rank, lb, ub)
            end if
#else
#error      "Unrecognized interface."
#endif
            if (present(isSorted_local)) then
                call setSorted(dataUnsorted, isSorted = isSorted_local)
            else
                call setSorted(dataUnsorted)
            end if
            assertion = assertion .and. GET_COMP(dataUnsorted(GET_INDEX(rank))) IS_EQUAL GET_COMP(selection)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("GET_COMP(dataUnsorted(GET_INDEX(rank)))")
                call test%disp%show( GET_COMP(dataUnsorted(GET_INDEX(rank))) )
                call test%disp%show("GET_COMP(selection)")
                call test%disp%show( GET_COMP(selection) )
                call test%disp%show("rank")
                call test%disp%show( rank )
                if (present(lb)) then
                    call test%disp%show("lb")
                    call test%disp%show( lb )
                end if
                if (present(ub)) then
                    call test%disp%show("ub")
                    call test%disp%show( ub )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"sort() must be able to sort input `contiguous` array of rank 1.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine runTestsWith

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function isAscending_local(a, b) result(sorted)
#if         setSelected_D0_SK_ENABLED || getSelected_D0_SK_ENABLED
            character(1,SKG), intent(in) :: a, b
#elif       setSelected_D1_SK_ENABLED || getSelected_D1_SK_ENABLED
            character(*,SKG), intent(in) :: a, b
#elif       setSelected_D1_IK_ENABLED || getSelected_D1_IK_ENABLED
            integer(IKG)    , intent(in) :: a, b
#elif       setSelected_D1_LK_ENABLED || getSelected_D1_LK_ENABLED
            logical(LKG)    , intent(in) :: a, b
#elif       setSelected_D1_CK_ENABLED || getSelected_D1_CK_ENABLED
            complex(CKG)    , intent(in) :: a, b
#elif       setSelected_D1_RK_ENABLED || getSelected_D1_RK_ENABLED
            real(RKG)       , intent(in) :: a, b
#elif       setSelected_D1_PSSK_ENABLED || getSelected_D1_PSSK_ENABLED
            type(css_pdt(SKG)) , intent(in) :: a, b
#else
#error      "Unrecognized interface."
#endif
            logical(LK) :: sorted
            sorted = GET_COMP(a) < GET_COMP(b)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function isDescending_local(a, b) result(sorted)
#if         setSelected_D0_SK_ENABLED || getSelected_D0_SK_ENABLED
            character(1,SKG), intent(in) :: a, b
#elif       setSelected_D1_SK_ENABLED || getSelected_D1_SK_ENABLED
            character(*,SKG), intent(in) :: a, b
#elif       setSelected_D1_IK_ENABLED || getSelected_D1_IK_ENABLED
            integer(IKG)    , intent(in) :: a, b
#elif       setSelected_D1_LK_ENABLED || getSelected_D1_LK_ENABLED
            logical(LKG)    , intent(in) :: a, b
#elif       setSelected_D1_CK_ENABLED || getSelected_D1_CK_ENABLED
            complex(CKG)    , intent(in) :: a, b
#elif       setSelected_D1_RK_ENABLED || getSelected_D1_RK_ENABLED
            real(RKG)       , intent(in) :: a, b
#elif       setSelected_D1_PSSK_ENABLED || getSelected_D1_PSSK_ENABLED
            type(css_pdt(SKG)) , intent(in) :: a, b
#else
#error      "Unrecognized interface."
#endif
            logical(LK) :: sorted
            sorted = GET_COMP(a) > GET_COMP(b)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  GET_INDEX
#undef  IS_EQUAL