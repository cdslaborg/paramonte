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
!>  This module contains implementations of the tests of the procedures under the generic interfaces of [pm_arrayRefill](@ref pm_arrayRefill).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Bypass gfortran bug 10-12.
#if     setRefilled_D1_SK_ENABLED || setRefilled_D2_SK_ENABLED || setRefilled_D3_SK_ENABLED
#define TYPE_KIND character(2,SKC) ::
#else
#define TYPE_KIND
#endif

#if     setRefilled_D1_LK_ENABLED || setRefilled_D2_LK_ENABLED || setRefilled_D3_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     setRefilled_D0_ENABLED || setRefilled_D1_ENABLED
#define SET_BND(X,lb,ub) X(lb : ub)
#define SET_DIM(X) X(:)
#define SET_SIZE(X) X
#elif   setRefilled_D2_ENABLED
#define SET_BND(X,lb,ub) X(lb(1) : ub(1), lb(2) : ub(2))
#define SET_SIZE(X) X(rank(array))
#define SET_DIM(X) X(:,:)
#elif   setRefilled_D3_ENABLED
#define SET_BND(X,lb,ub) X(lb(1) : ub(1), lb(2) : ub(2), lb(3) : ub(3))
#define SET_SIZE(X) X(rank(array))
#define SET_DIM(X) X(:,:,:)
#else
#error  "Unrecognized interface."
#endif

#if     setRefilled_D0_ENABLED
        integer(IK)     , parameter     :: lbmin = +1_IK, ubmax = +15_IK
        character(:,SKC), allocatable   :: array, arrayInit, array_ref, lower, upper
        character(1,SKC), parameter     :: fill = SKC_"-"
#define GET_UBOUND(X) len(X, kind = IK)
#define GET_LBOUND(X) 1
#define ALL
#else
#define GET_LBOUND(X) lbound(X, kind = IK)
#define GET_UBOUND(X) ubound(X, kind = IK)
        integer(IK)     , parameter     :: lbmin = -5_IK, ubmax = +10_IK
#if     setRefilled_D1_SK_ENABLED || setRefilled_D2_SK_ENABLED || setRefilled_D3_SK_ENABLED
        character(2,SKC), allocatable   :: SET_DIM(array), SET_DIM(arrayInit), SET_DIM(array_ref)
        character(2,SKC), parameter     :: fill = SKC_"--", lower = SKC_"aa", upper = SKC_"zz"
#elif   setRefilled_D1_IK_ENABLED || setRefilled_D2_IK_ENABLED || setRefilled_D3_IK_ENABLED
        integer(IKC)    , allocatable   :: SET_DIM(array), SET_DIM(arrayInit), SET_DIM(array_ref)
        integer(IKC)    , parameter     :: fill = huge(0_IKC), lower = -huge(0_IKC), upper = huge(0_IKC)
#elif   setRefilled_D1_LK_ENABLED || setRefilled_D2_LK_ENABLED || setRefilled_D3_LK_ENABLED
        logical(LKC)    , allocatable   :: SET_DIM(array), SET_DIM(arrayInit), SET_DIM(array_ref)
        logical(LKC)    , parameter     :: fill = .false._LKC, lower = .false._LKC, upper = .true._LKC
#elif   setRefilled_D1_CK_ENABLED || setRefilled_D2_CK_ENABLED || setRefilled_D3_CK_ENABLED
        complex(CKC)    , allocatable   :: SET_DIM(array), SET_DIM(arrayInit), SET_DIM(array_ref)
        complex(CKC)    , parameter     :: fill = cmplx(huge(0._CKC), huge(0._CKC), kind = CKC)
        complex(CKC)    , parameter     :: lower = -fill, upper = fill
#elif   setRefilled_D1_RK_ENABLED || setRefilled_D2_RK_ENABLED || setRefilled_D3_RK_ENABLED
        real(RKC)       , allocatable   :: SET_DIM(array), SET_DIM(arrayInit), SET_DIM(array_ref)
        real(RKC)       , parameter     :: fill = huge(0._RKC)
        real(RKC)       , parameter     :: lower = -fill, upper = fill
#else
#error  "Unrecognized interface."
#endif
#endif
        character(127, SK)  :: errmsg
        logical(LK)         :: failed
        integer(IK)         :: SET_SIZE(lb)
        integer(IK)         :: SET_SIZE(ub)
        integer(IK)         :: SET_SIZE(lbc)
        integer(IK)         :: SET_SIZE(lbold)
        integer(IK)         :: SET_SIZE(ubold)
        integer(IK)         :: SET_SIZE(lbcold)
        integer(IK)         :: SET_SIZE(ubcold)
        logical(LK)         :: assumedSize
        integer             :: itest
        type(display_type)  :: disp
        disp = display_type()
        assertion = .true._LK

        do itest = 1, 100
            call runTestsWith()
            call runTestsWith(failed)
            call runTestsWith(errmsg = errmsg)
            call runTestsWith(failed, errmsg)
        end do

        ! Test with unallocated input `array`.
        call runTestsWithUnalloc()
        call runTestsWithUnalloc(failed)
        call runTestsWithUnalloc(errmsg = errmsg)
        call runTestsWithUnalloc(failed, errmsg)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine checkFailure(line, failed, errmsg)
            integer, intent(in) :: line
            logical(LK), intent(in), optional :: failed
            character(*, SK), intent(in), optional :: errmsg
            if (present(failed)) then
                assertion = assertion .and. .not. failed
                call test%assert(assertion, SK_"The `array` resizing must not fail with present(failed), present(errmsg) = "//getStr([present(failed), present(errmsg)]), int(line, IK))
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWithUnalloc(failed, errmsg)
            logical(LK), intent(out), optional :: failed
            character(*, SK), intent(out), optional :: errmsg
            if (allocated(array_ref)) deallocate(array_ref)
            lb = 1_IK
            call setUnifRand(ub, lb - 1_IK, ubmax)
#if         setRefilled_D0_ENABLED
            array_ref = repeat(fill, ub)
#else
            allocate(TYPE_KIND SET_BND(array_ref, lb, ub))
            array_ref = fill
#endif
            if (allocated(array)) deallocate(array)
            call setRefilled(array, fill, ub - lb + 1_IK, failed, errmsg)
            call checkFailure(__LINE__, failed, errmsg)
            assertion = assertion .and. ALL(GET_LBOUND(array) == GET_LBOUND(array_ref))
            call test%assert(assertion, SK_"The lower bounds of the output `array` must be correctly set when the input `array` is unallocated with present(failed), present(errmsg), LBOUND(array), LBOUND(array_ref) = "// & ! LCOV_EXCL_LINE
            getStr([present(failed), present(errmsg)])//SK_", "//getStr([GET_LBOUND(array), GET_LBOUND(array_ref)]), int(__LINE__, IK))
            assertion = assertion .and. ALL(GET_UBOUND(array) == GET_UBOUND(array_ref))
            call test%assert(assertion, SK_"The upper bounds of the output `array` must be correctly set when the input `array` is unallocated with present(failed), present(errmsg), UBOUND(array), UBOUND(array_ref) = "// & ! LCOV_EXCL_LINE
            getStr([present(failed), present(errmsg)])//SK_", "//getStr([GET_UBOUND(array), GET_UBOUND(array_ref)]), int(__LINE__, IK))
            assertion = assertion .and. ALL(array IS_EQUAL array_ref)
#if         setRefilled_D0_ENABLED
            assertion = assertion .and. len_trim(array) == len_trim(array_ref)
#endif
            call test%assert(assertion, SK_"The contents of the output `array` must be correctly set when the input `array` is unallocated with present(failed), present(errmsg) = "//getStr([present(failed), present(errmsg)]), int(__LINE__, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith(failed, errmsg)

            logical(LK), intent(out), optional :: failed
            character(*, SK), intent(out), optional :: errmsg

            if (allocated(arrayInit)) deallocate(arrayInit)
            if (allocated(array_ref)) deallocate(array_ref)

            ! Test the missing `size` interface.
            assumedSize = logical(getUnifRand(0., 1.) < 0.1, LK)
#if         setRefilled_D0_ENABLED
            lb = 1_IK
            lbold = 1_IK
            call setUnifRand(ubold, lbold, ubmax)
            if (assumedSize) then
                ub = ubold * 2_IK
            else
                call setUnifRand(ub, lb, 2 * ubmax)
            end if
            allocate(character(ub,SKC) :: array_ref)
            allocate(character(ubold,SKC) :: arrayInit)
            lower = repeat(SKC_"a", len(arrayInit))
            upper = repeat(SKC_"z", len(arrayInit))
#else
            call setUnifRand(lbold, lbmin, ubmax)
            call setUnifRand(ubold, lbold, ubmax)
            lb = lbold ! call setUnifRand(lb, -15_IK, +5_IK)
            if (assumedSize) then
                ub = lb + (ubold - lbold + 1_IK) * 2_IK - 1_IK
            else
                call setUnifRand(ub, lb, 2 * ubmax)
            end if
            allocate(TYPE_KIND SET_BND(array_ref, lb, ub))
            allocate(TYPE_KIND SET_BND(arrayInit, lbold, ubold))
#endif
            call setUnifRand(arrayInit, lower, upper)
            if (all([lb] <= [lbold]) .and. all([ubold] <= [ub]) .and. getUnifRand()) then
                if (assumedSize) then
                    ! Test the expansion interface without `size`.
                    lbc = lbold
                    lbcold = lbold
                    ubcold = ubold
                    if (allocated(array)) deallocate(array)
                    allocate(array, source = arrayInit)
                    call setCoreHalo(array_ref, array, fill, lbc - lb)
                    !write(*,*) "lb, ub", lb, ub
                    call setRefilled(array, fill, failed, errmsg)
                    call report(__LINE__, failed, errmsg)
                end if
                ! Test the expansion interface.
                lbc = lbold
                lbcold = lbold
                ubcold = ubold
                if (allocated(array)) deallocate(array)
                allocate(array, source = arrayInit)
                call setCoreHalo(array_ref, array, fill, lbc - lb)
                !write(*,*) "lb, ub", lb, ub
                call setRefilled(array, fill, ub - lb + 1_IK, failed, errmsg)
                call report(__LINE__, failed, errmsg)
                ! Test the expansion + shift interface.
                call setUnifRand(lbc, lb, ub - (ubcold - lbcold))
                !write(*,*) "lb, ub, lbc, lbcold, ubcold", lb, ub, lbc, lbcold, ubcold
                if (allocated(array)) deallocate(array)
                allocate(array, source = arrayInit)
                call setCoreHalo(array_ref, array, fill, lbc - lb)
                call setRefilled(array, fill, ub - lb + 1_IK, lbc, failed, errmsg)
                call report(__LINE__, failed, errmsg)
            end if
            ! Test the expansion/contraction + shift + subset interface.
            if (allocated(array)) deallocate(array)
            allocate(array, source = arrayInit)
            call setUnifRand(lbcold, lbold, ubold)
            call setUnifRand(ubcold, lbcold, min(ubold, lbcold + min(ubold - lbold, ub - lb)))
            call setUnifRand(lbc, lb, ub - (ubcold - lbcold))
            !write(*,*) "GET_LBOUND(array), GET_UBOUND(array), lbcold, ubcold", GET_LBOUND(array), GET_UBOUND(array), lbcold, ubcold
            !write(*,*) "SET_BND(array, lbcold, ubcold)", SET_BND(array, lbcold, ubcold)
            !write(*,*) "array_ref", array_ref
            call setCoreHalo(array_ref, SET_BND(array, lbcold, ubcold), fill, lbc - lb)
            !write(*,*) "array_ref, lb, ub, lbc, lbcold, ubcold", array_ref, lb, ub, lbc, lbcold, ubcold
            call setRefilled(array, fill, ub - lb + 1_IK, lbc, lbcold, ubcold, failed, errmsg)
            call report(__LINE__, failed, errmsg)

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line, failed, errmsg)
            integer, intent(in) :: line
            logical(LK), intent(in), optional :: failed
            character(*, SK), intent(in), optional :: errmsg
            call checkFailure(line, failed, errmsg)
            assertion = assertion .and. ALL(GET_LBOUND(array) == GET_LBOUND(array_ref))
            call display()
            call test%assert(assertion, SK_"The lower bounds of the output `array` must be correctly set with present(failed), present(errmsg) = "//getStr([present(failed), present(errmsg)]), int(line, IK))
            assertion = assertion .and. ALL(GET_UBOUND(array) == GET_UBOUND(array_ref))
            call display()
            call test%assert(assertion, SK_"The upper bounds of the output `array` must be correctly set with present(failed), present(errmsg) = "//getStr([present(failed), present(errmsg)]), int(line, IK))
            assertion = assertion .and. ALL(array IS_EQUAL array_ref)
            call display()
            call test%assert(assertion, SK_"Call to setRefilled() must correctly rebind and refill `array` with present(failed), present(errmsg) = "//getStr([present(failed), present(errmsg)]), int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine display()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip()
                call disp%show("rank(array)")
                call disp%show( rank(array) )
                call disp%show("lbold")
                call disp%show( lbold )
                call disp%show("ubold")
                call disp%show( ubold )
                call disp%show("lb")
                call disp%show( lb )
                call disp%show("ub")
                call disp%show( ub )
                call disp%show("lbc")
                call disp%show( lbc )
                call disp%show("lbcold")
                call disp%show( lbcold )
                call disp%show("ubcold")
                call disp%show( ubcold )
                call disp%show("arrayInit")
                call disp%show( arrayInit )
                call disp%show("array_ref")
                call disp%show( array_ref )
                call disp%show("array")
                call disp%show( array )
                call disp%show("array == array_ref")
                call disp%show( array IS_EQUAL array_ref )
                call disp%skip()
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GET_LBOUND
#undef GET_UBOUND
#undef TYPE_KIND
#undef IS_EQUAL
#undef SET_SIZE
#undef SET_BND
#undef SET_DIM
#undef ALL
