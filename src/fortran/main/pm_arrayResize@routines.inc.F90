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
!>  This file contains the implementation details of the routines under the generic interfaces of [pm_arrayResize](@ref pm_arrayResize).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the procedure names.
#if     setResized_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setResized()"
#elif   setRefilled_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setRefilled()"
#elif   setRebound_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setRebound()"
#elif   setRebilled_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@setRebilled()"
#else
#error  "Unrecognized interface."
#endif
        integer :: stat
        ! Set the lower bound of the new `array` for fixed lower bound routines.
#if     setResized_ENABLED || setRefilled_ENABLED
#define lb lbold
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Set the dimensionality of `array` and the allocation dimension and define `array` bounds and copy slices.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     D0_ENABLED && (setResized_ENABLED || setRefilled_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define ALL(X) X
#define GET_SHAPE(X)len(X, IK)
#define ARRAY_SLICE array(lbcold : ubcold)
#define SET_DIM(OBJECT) character(ub,SKG) :: OBJECT
#define TEMP_SLICE temp(lbc : lbc - lbcold + ubcold)
        integer(IK), parameter :: lbold = 1_IK
        integer(IK) :: ubold, ub
#if     SLDD_ENABLED || SDDD_ENABLED || DDDD_ENABLED
        integer(IK) :: lbcold, ubcold
#if     SDDD_ENABLED || DDDD_ENABLED
        integer(IK) :: lbc
#endif
#elif   !SLLU_ENABLED
#error  "Unrecognized interface."
#endif
        ubold = len(array, IK)

        !%%%%%%%%%
#elif   D1_ENABLED
        !%%%%%%%%%

#define ALL(X) X
#define GET_SHAPE(X)shape(X, IK)
#define SET_DIM(OBJECT) OBJECT(lb : ub)
#define ARRAY_SLICE array(lbcold : ubcold)
#define TEMP_SLICE temp(lbc : lbc - lbcold + ubcold)
#if     setResized_ENABLED || setRefilled_ENABLED
        integer(IK) :: ub
#endif
        integer(IK) :: lbold, ubold
#if     SLDD_ENABLED || SDDD_ENABLED || DDDD_ENABLED
        integer(IK) :: lbcold, ubcold
#if     SDDD_ENABLED || DDDD_ENABLED
        integer(IK) :: lbc
#endif
#elif   !SLLU_ENABLED
#error  "Unrecognized interface."
#endif
#define GET_BOUND(BOUND, X) BOUND(X, 1, kind = IK)

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   D2_ENABLED || D3_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

#define GET_SHAPE(X)shape(X, IK)
#if     setResized_ENABLED || setRefilled_ENABLED
        integer(IK) :: ub(rank(array))
#endif
        integer(IK) :: lbold(rank(array)), ubold(rank(array))
#if     SLDD_ENABLED || SDDD_ENABLED || DDDD_ENABLED
        integer(IK) :: lbcold(rank(array)), ubcold(rank(array))
#if     SDDD_ENABLED || DDDD_ENABLED
        integer(IK) :: lbc(rank(array))
#endif
#elif   !SLLU_ENABLED
#error  "Unrecognized interface."
#endif
#define GET_BOUND(BOUND, X) BOUND(X, kind = IK)
#if     D2_ENABLED
#define SET_DIM(OBJECT) OBJECT(lb(1) : ub(1), lb(2) : ub(2))
#define ARRAY_SLICE array(lbcold(1) : ubcold(1), lbcold(2) : ubcold(2))
#define TEMP_SLICE temp(lbc(1) : lbc(1) - lbcold(1) + ubcold(1), lbc(2) : lbc(2) - lbcold(2) + ubcold(2))
#elif   D3_ENABLED
#define SET_DIM(OBJECT) OBJECT(lb(1) : ub(1), lb(2) : ub(2), lb(3) : ub(3))
#define ARRAY_SLICE array(lbcold(1) : ubcold(1), lbcold(2) : ubcold(2), lbcold(3) : ubcold(3))
#define TEMP_SLICE temp(lbc(1) : lbc(1) - lbcold(1) + ubcold(1), lbc(2) : lbc(2) - lbcold(2) + ubcold(2), lbc(3) : lbc(3) - lbcold(3) + ubcold(3))
#else
#error  "Unrecognized interface."
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Bypass the gfortran allocation statement error for objects of type `character` of non-zero rank.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     SK_ENABLED && !D0_ENABLED && __GFORTRAN__
#define TYPE_OF_ARRAY character(len(array,IK),SKG) ::
#else
#define TYPE_OF_ARRAY
#endif
        ! Check the consistency of the fill length with the string array elements length.
#if     SK_ENABLED && !D0_ENABLED && (setRefilled_ENABLED || setRebilled_ENABLED)
        CHECK_ASSERTION(__LINE__, len(fill, IK) <= len(array, IK), PROCEDURE_NAME//SK_": The condition `len(fill) <= len(array)` must hold. len(fill), len(array) = "//getStr([len(fill, IK), len(array, IK)]))
#endif
        ! Define the allocation statement.
!#define SET_ALLOCATION(OBJECT) \
!if (present(failed)) then; \
!allocate(SET_DIM(OBJECT), stat); \
!if (stat /= 0) return; \
!else; \
!allocate(SET_DIM(OBJECT)); \
!end if;
        !%%%%%%%%%%%
#if     SDDD_ENABLED
        !%%%%%%%%%%%

        ! Check the allocation status.
        if (.not. allocated(array)) then
#if         setResized_ENABLED || setRefilled_ENABLED
#if         !D0_ENABLED
            lb = 1_IK
#endif
            ub = size
#endif
            !SET_ALLOCATION(array)
            if (present(failed)) then
                if (present(errmsg)) then
                    allocate(TYPE_OF_ARRAY SET_DIM(array), stat = stat, errmsg = errmsg)
                else
                    allocate(TYPE_OF_ARRAY SET_DIM(array), stat = stat)
                end if
                failed = logical(stat /= 0, LK)
                if (failed) return ! LCOV_EXCL_LINE
            else
                allocate(TYPE_OF_ARRAY SET_DIM(array))
            end if
#if         setRefilled_ENABLED && D0_ENABLED
            block
                integer(IK) :: i
                do concurrent(i = 1 : len(array, IK))
                    array(i:i) = fill
                end do
            end block
#elif       setRefilled_ENABLED || setRebilled_ENABLED
            array = fill
#elif       !(setResized_ENABLED || setRebound_ENABLED)
#error      "Unrecognized interface."
#endif
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   DDDD_ENABLED || SLDD_ENABLED || SLLU_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, allocated(array), PROCEDURE_NAME//SK_": The condition `allocated(array)` must hold.")
#if     setRefilled_ENABLED && D0_ENABLED && SK_ENABLED
        CHECK_ASSERTION(__LINE__, len(fill) <= len(array), PROCEDURE_NAME//SK_": The condition `len(fill) <= len(array)` must hold. len(fill), len(array) = "//getStr([len(fill), len(array)]))
#endif
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

#if     !D0_ENABLED
        lbold = GET_BOUND(lbound, array)
        ubold = GET_BOUND(ubound, array)
#endif
        ! Set the new upper bound of `array` for `setResized` and `setRefilled`.
#if     DDDD_ENABLED && (setResized_ENABLED || setRefilled_ENABLED)
        CHECK_ASSERTION(__LINE__, all([0_IK < GET_SHAPE(array)]), PROCEDURE_NAME//SK_": The condition `all([0 < len/shape(array))` must hold when the input argument `size` is missing. len/shape(array) = "//getStr(GET_SHAPE(array)))
        ub = lbold - 1_IK + 2_IK * (ubold - lbold + 1_IK)
#elif   setResized_ENABLED || setRefilled_ENABLED
        ub = lbold - 1_IK + size
#endif
        ! Check or set the old contents bounds and new contents offset.
#if     SLLU_ENABLED
        ! Check the contents offset.
        ! \bug Bypass the Intel compiler bug in processing multiple `CHECK_ASSERTION`
        ! macros in a single routine in `debug` compile mode by merging all `CHECK_ASSERTION` macros.
        CHECK_ASSERTION(__LINE__, ALL(lbold <= lbcold .and. lbcold <= ubold), PROCEDURE_NAME//SK_": The condition `all(lbound(array) <= lbcold .and. lbcold <= ubound(array))` must hold. rank(array), lbound(array), lbcold, ubound(array) = "//getStr([int(rank(array), IK), lbold, lbcold, ubold]))
        CHECK_ASSERTION(__LINE__, ALL(lbold <= ubcold .and. ubcold <= ubold), PROCEDURE_NAME//SK_": The condition `all(lbound(array) <= ubcold .and. ubcold <= ubound(array))` must hold. rank(array), lbound(array), ubcold, ubound(array) = "//getStr([int(rank(array), IK), lbold, ubcold, ubold]))
#elif   SLDD_ENABLED
        lbcold = max(lbold, lb)
        ubcold = lbcold + min(ubold - lbcold, ub - lbc)
        ! Check the contents lower bound.
        ! \bug Bypass the Intel compiler bug in processing multiple `CHECK_ASSERTION`
        ! macros in a single routine in `debug` compile mode by merging all `CHECK_ASSERTION` macros.
        CHECK_ASSERTION(__LINE__, ALL(lb <= lbc), PROCEDURE_NAME//SK_": The condition `all(lb <= lbc)` must hold where `lb` is the lower bound of the output `array`. rank(array), lb, lbc = "//getStr([int(rank(array), IK), lb, lbc]))
        CHECK_ASSERTION(__LINE__, ALL(lbc - lbcold + ubcold <= ub), PROCEDURE_NAME//SK_": The condition `all(lbc - lbcold + ubcold <= ub)` must hold with `ub` as the output `array` ubound. rank(array), lbc, lbcold, ubcold, ub = "//getStr([int(rank(array), IK), lbc, lbcold, ubcold, ub]))
#else
        lbcold = max(lbold, lb)
        ubcold = min(ubold, ub)
        lbc = lbcold
#endif
        ! Check the output `array` size.
        CHECK_ASSERTION(__LINE__, ALL(0_IK <= ub - lb + 1_IK), PROCEDURE_NAME//SK_": The condition `all(0_IK <= ub - lb + 1_IK)` must hold where `lb, ub` are the lower and upper bounds of the output `array`. lb, ub = "//getStr([lb, ub]))
        !SET_ALLOCATION(temp)
        if (present(failed)) then
            if (present(errmsg)) then
                allocate(TYPE_OF_ARRAY SET_DIM(temp), stat = stat, errmsg = errmsg)
            else
                allocate(TYPE_OF_ARRAY SET_DIM(temp), stat = stat)
            end if
            failed = logical(stat /= 0, LK)
            if (failed) return ! LCOV_EXCL_LINE
        else
            allocate(TYPE_OF_ARRAY SET_DIM(temp))
        end if
        ! Copy contents.
#if     setResized_ENABLED || setRebound_ENABLED
        TEMP_SLICE = ARRAY_SLICE
#elif   setRefilled_ENABLED || setRebilled_ENABLED
        call setCoreHalo(temp, ARRAY_SLICE, fill, lbc - lb)
#else
#error  "Unrecognized interface."
#endif
        call move_alloc(from = temp, to = array)
#undef  SET_ALLOCATION
#undef  TYPE_OF_ARRAY
#undef  ARRAY_SLICE
#undef  TEMP_SLICE
#undef  GET_SHAPE
#undef  GET_BOUND
#undef  SET_DIM
#undef  lbcold
#undef  ubcold
#undef  SIZE
#undef  lbc
#undef  ALL
#undef  lb