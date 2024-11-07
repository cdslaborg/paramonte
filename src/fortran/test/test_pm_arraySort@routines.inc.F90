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
!>  This module contains implementations of the tests of the procedures of [pm_arraySort](@ref pm_arraySort).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%
#if     setSorted_ENABLED
        !%%%%%%%%%%%%%%%%

#if     Ind_ENABLED
        integer(IK), allocatable :: index(:)
#else
#if     Arr_ENABLED && Qsorti_ENABLED
#define METHOD_TYPE qsorti_type
#elif   Arr_ENABLED && Qsortr_ENABLED
#define METHOD_TYPE qsortr_type
#elif   Arr_ENABLED && Qsortrdp_ENABLED
#define METHOD_TYPE qsortrdp_type
#elif   Arr_ENABLED && Bubble_ENABLED
#define METHOD_TYPE bubble_type
#elif   Arr_ENABLED && Heapi_ENABLED
#define METHOD_TYPE heapi_type
#elif   Arr_ENABLED && Heapr_ENABLED
#define METHOD_TYPE heapr_type
#elif   Arr_ENABLED && Insertionl_ENABLED
#define METHOD_TYPE insertionl_type
#elif   Arr_ENABLED && Insertionb_ENABLED
#define METHOD_TYPE insertionb_type
#elif   Arr_ENABLED && Merger_ENABLED
#define METHOD_TYPE merger_type
#elif   Arr_ENABLED && Selection_ENABLED
#define METHOD_TYPE selection_type
#elif   Arr_ENABLED && Shell_ENABLED
#define METHOD_TYPE shell_type
#else
#error  "Unrecognized Interface."
#endif
        type(METHOD_TYPE), parameter :: method = METHOD_TYPE()
#endif
        ! test data.
        integer(IK), parameter :: NDATA = 1000_IK
#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE(array) len(array, kind = IK)
#define GET_INDEX(i) i:i
        character(0,SKG)                        :: empty
        character(:,SKG)        , allocatable   :: array
        character(1,SKG)        , parameter     :: LOWER = SKG_"a", UPPER = SKG_"z"
#else
#define GET_INDEX(i) i
#define GET_SIZE(array) size(array, kind = IK)
#if     SK_ENABLED && D1_ENABLED
        character(2,SKG)                        :: empty(0)
        character(2,SKG)        , allocatable   :: array(:)
        character(2,SKG)        , parameter     :: LOWER = SKG_"aA", UPPER = SKG_"zZ"
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)                            :: empty(0)
        integer(IKG)            , allocatable   :: array(:)
        integer(IKG)            , parameter     :: LOWER = 1_IKG, UPPER = huge(1_IKG)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)                            :: empty(0)
        logical(LKG)            , allocatable   :: array(:)
        logical(LKG)            , parameter     :: LOWER = .false._LKG, UPPER = .true._LKG
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)                            :: empty(0)
        complex(CKG)            , allocatable   :: array(:)
        complex(CKG)            , parameter     :: LOWER = cmplx(1._CKG, -huge(1._CKG), CKG), UPPER = cmplx(huge(1._CKG), -1._CKG, CKG)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)                               :: empty(0)
        real(RKG)               , allocatable   :: array(:)
        real(RKG)               , parameter     :: LOWER = 1._RKG, UPPER = huge(1._RKG)
#elif   PSSK_ENABLED && D1_ENABLED
        integer(IK) :: i
        type(css_pdt(SKG))                      :: empty(0)
        type(css_pdt(SKG))      , allocatable   :: array(:)
        do i = 1, NDATA
            allocate(character(SKG,2) :: array(i)%val)
            call setUnifRand(array(i)%val, SKG_"AA", SKG_"ZZ")
        end do
        !write(*,"(1(g0,:,' '))") array
        !error stop
#else
#error  "Unrecognized Interface."
#endif
#endif
        assertion = .true._LK
        call runTestsWith()
        call runTestsWith(isSortedElement)

    contains

        subroutine runTestsWith(isSortedElement)
            procedure(logical(LK)), optional :: isSortedElement
            logical(LK) :: isPresentMethod
            integer(IK) :: i, lenArray
            do i = 1_IK, 200_IK
                isPresentMethod = getUnifRand()
                lenArray = getUnifRand(0, 500)
                if (allocated(array)) deallocate(array)
#if             SK_ENABLED && D0_ENABLED
                allocate(character(lenArray,SKG) :: array)
                !call setUnifRand(array, repeat(SKG_"a", len(array)), repeat(SKG_"z", len(array)))
                call setUnifRand(array)
#else
                allocate(array(1 : lenArray))
                call setUnifRand(array)!, LOWER, UPPER) ! bounds are commented out due to a potential gfortran-13 bug. See commnets below for info.
#endif
#if             Ind_ENABLED
                call setResized(index, lenArray)
                if (present(isSortedElement)) then
                    call setSorted(array, index, isSortedElement)
                    assertion = assertion .and. isDescending(getRemapped(array, index))
                else
                    call setSorted(array, index)
                    assertion = assertion .and. isAscending(getRemapped(array, index))
                end if
#elif           Arr_ENABLED
                if (present(isSortedElement)) then
                    if (isPresentMethod) then
                        call setSorted(array, isSortedElement, method)
                    else
                        call setSorted(array, isSortedElement)
                    end if
                    assertion = assertion .and. isDescending(array)
                else
                    if (isPresentMethod) then
                        call setSorted(array, method)
                    else
                        call setSorted(array)
                    end if
                    assertion = assertion .and. isAscending(array)
                end if
#endif
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    call test%disp%skip()
                    call test%disp%show("present(isSortedElement)")
                    call test%disp%show( present(isSortedElement) )
#if                 Ind_ENABLED
                    call test%disp%show("index")
                    call test%disp%show( index )
#elif               Arr_ENABLED
                    call test%disp%show("isPresentMethod")
                    call test%disp%show( isPresentMethod )
#endif
                    call test%disp%show("array")
                    call test%disp%show( array )
                    call test%disp%skip()
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, SK_"setSorted() must correctly sort the input array or its index.", int(__LINE__, IK))
            end do
#if         Ind_ENABLED
            if (present(isSortedElement)) then
                if (isPresentMethod) then
                    call setSorted(empty, isSortedElement)
                else
                    call setSorted(empty, isSortedElement)
                end if
            else
                if (isPresentMethod) then
                    call setSorted(empty)
                else
                    call setSorted(empty)
                end if
            end if
            call test%assert(assertion, SK_"setSorted() must handle empty array sorting with present(isSortedElement) = "//getStr(present(isSortedElement)), int(__LINE__, IK))
#elif       Arr_ENABLED
            if (present(isSortedElement)) then
                if (isPresentMethod) then
                    call setSorted(empty, isSortedElement, method)
                else
                    call setSorted(empty, isSortedElement)
                end if
            else
                if (isPresentMethod) then
                    call setSorted(empty, method)
                else
                    call setSorted(empty)
                end if
            end if
            call test%assert(assertion, SK_"setSorted() must handle empty array sorting with present(isSortedElement) = "//getStr(present(isSortedElement)), int(__LINE__, IK))
#endif
        end subroutine runTestsWith

        pure function isSortedElement(a, b) result(sorted)
            logical(LK) :: sorted
#if         SK_ENABLED && D0_ENABLED
            character(1,SKG)        , intent(in) :: a, b
#elif       SK_ENABLED && D1_ENABLED
            character(*,SKG)        , intent(in) :: a, b
#elif       IK_ENABLED && D1_ENABLED
            integer(IKG)            , intent(in) :: a, b
#elif       LK_ENABLED && D1_ENABLED
            logical(LKG)            , intent(in) :: a, b
#elif       CK_ENABLED && D1_ENABLED
            complex(CKG)            , intent(in) :: a, b
#elif       RK_ENABLED && D1_ENABLED
            real(RKG)               , intent(in) :: a, b
#elif       PSSK_ENABLED && D1_ENABLED
            type(css_pdt(SKG))      , intent(in) :: a, b
#else
#error      "Unrecognized interface."
#endif
            sorted = a > b
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isSorted_ENABLED || isAscending_ENABLED || isDescending_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     isSorted_ENABLED
        use pm_arraySort, only: isOrdered => isSorted
#elif   isAscending_ENABLED
        use pm_arraySort, only: isOrdered => isAscending
#elif   isDescending_ENABLED
        use pm_arraySort, only: isOrdered => isDescending
#else
#error  "Unrecognized interface."
#endif
        integer(IK)     , parameter :: NDATA = 1000_IK
#if     D0_ENABLED && SK_ENABLED
        character(0,SKG)            :: Empty
        character(NDATA)            :: dataUnsorted
        call setUnifRand(dataUnsorted, repeat(SKG_"A", len(dataUnsorted,IK)), repeat(SKG_"Z", len(dataUnsorted,IK)))
#elif   D1_ENABLED && SK_ENABLED
        character(2,SKG)            :: Empty(0)
        character(2,SKG)            :: dataUnsorted(NDATA)
        call setUnifRand(dataUnsorted, SKG_"AA", SKG_"ZZ")
#elif   D1_ENABLED && IK_ENABLED
        integer(IKG)                :: Empty(0)
        integer(IKG)                :: dataUnsorted(NDATA)
        call setUnifRand(dataUnsorted, 1_IKG, huge(1_IKG))
#elif   D1_ENABLED && LK_ENABLED
        logical(LKG)                :: Empty(0)
        logical(LKG)                :: dataUnsorted(NDATA)
        call setUnifRand(dataUnsorted)
#elif   D1_ENABLED && CK_ENABLED
        complex(CKG)    , parameter :: LB = cmplx(0., -9., CKG), UB = cmplx(9., 0., CKG)
        complex(CKG)                :: dataUnsorted(NDATA)
        complex(CKG)                :: Empty(0)
        ! \bug
        ! gfortran-13 release mode heap-memory nocheck shared-lib passes NAN values to `setUnifRand()` for some of the input bounds.
        ! This caused infinite loops in `setUnifRand()`. Thus, the implementation of `setUnifRand()` was modified to handle NANs gracefully.
        ! The root cause of this remains unknown. For now, the bounds are excluded to allow testing to proceed.
        ! This issue could be related to the `elemental` attribute of `setUnifRand()` as similar problems
        ! have been also observed for other `elemental` routines.
        call setUnifRand(dataUnsorted)!, LB, UB)
#elif   D1_ENABLED && RK_ENABLED
        real(RKG)                   :: dataUnsorted(NDATA)
        real(RKG)                   :: Empty(0)
        call setUnifRand(dataUnsorted, 1._RKG, huge(1._RKG))
#elif   D1_ENABLED && PSSK_ENABLED
        integer(IK) :: i
        type(css_pdt(SKG))          :: Empty(0)
        type(css_pdt(SKG))          :: dataUnsorted(NDATA)
        do i = 1, NDATA
            allocate(character(SKG,2) :: dataUnsorted(i)%val)
            call setUnifRand(dataUnsorted(i)%val, SKG_"AA", SKG_"ZZ")
        end do
#else
#error  "Unrecognized Interface."
#endif
        ! Test for contiguous input arrays.
        ! The following tests may, in extremely rare conditions fail, for example, when the generated random array is fully sorted.

        !call random_seed()
        assertion = isOrdered(Empty)
        call test%assert(assertion, SK_"isOrdered() must return `.true.` for an input `contiguous` array of rank 1 of length 0.", int(__LINE__, IK))

        assertion = .not. isOrdered(dataUnsorted)
        call test%assert(assertion, SK_"isOrdered() must return `.false.` for an input contiguous unsorted array.", int(__LINE__, IK))

        call setSorted(dataUnsorted)
#if     isSorted_ENABLED || isAscending_ENABLED
        assertion = isOrdered(dataUnsorted)
#elif   isDescending_ENABLED
        assertion = .not. isOrdered(dataUnsorted)
#endif
        call test%assert(assertion, SK_"isOrdered() must return a valid result for an input contiguous ascending-sorted array.", int(__LINE__, IK))
        dataUnsorted = getReversed(dataUnsorted) ! This is called due to a GFortran bug as of GFortran version 10.3.
#if     isSorted_ENABLED || isDescending_ENABLED
        assertion = isOrdered(dataUnsorted)
#elif   isAscending_ENABLED
        assertion = .not. isOrdered(dataUnsorted)
#endif
        call test%assert(assertion, SK_"isOrdered() must return a valid result for an input contiguous descending-sorted array.", int(__LINE__, IK))

#if     D0_ENABLED && SK_ENABLED
        dataUnsorted(:) = repeat(dataUnsorted(1:1), len(dataUnsorted, kind = IK))
#else
        dataUnsorted(:) = dataUnsorted(1)
#endif
        assertion = isOrdered(dataUnsorted)
        call test%assert(assertion, SK_"isOrdered() must return `.true.` for an input contiguous identically-valued array.", int(__LINE__, IK))
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef METHOD_TYPE
#undef COMPONENT
#undef GET_INDEX
#undef GET_SIZE
#undef METHOD