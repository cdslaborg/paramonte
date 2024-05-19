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
!>  This file contains the implementations of the tests of module [pm_mathlinSpace](@ref pm_mathlinSpace).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPE_KIND complex(TKG)
#elif   RK_ENABLED
#define TYPE_KIND real(TKG)
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getLinSpace_ENABLED || setLinSpace_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: x1, x2
        real(TKG), parameter :: TOL = epsilon(0._TKG) * 10
        TYPE_KIND, allocatable :: linSpace(:)
        TYPE_KIND, allocatable :: linSpace_ref(:)
        real(TKG), allocatable :: diff(:)
        integer(IK) :: sign

        assertion = .true._LK

        sign = 1_IK
        call testWith()

        sign = -1_IK
        call testWith()

    contains

        subroutine testWith()

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (0._TKG, 10._TKG)
            x2 = (10._TKG, 0._TKG)
#elif       RK_ENABLED
            x1 = 0._TKG
            x2 = 10._TKG
#endif
            allocate(linSpace_ref(0))

            call report()
            call test%assert(assertion, SK_"getLinSpace() must return an empty `linSpace` when size of `linSpace` is zero.")

            call report(fopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an empty `linSpace` when size of `linSpace` is zero with `fopen = .false.`.")

            call report(lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an empty `linSpace` when size of `linSpace` is zero with `lopen = .false.`.")

            call report(fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an empty `linSpace` when size of `linSpace` is zero with `fopen = .false._LK, lopen = .false.`.")

            call report(fopen = .false._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an empty `linSpace` when size of `linSpace` is zero with `fopen = .false._LK, lopen = .true.`.")

            call report(fopen = .true._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an empty `linSpace` when size of `linSpace` is zero with `fopen = .true._LK, lopen = .false.`.")

            call report(fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an empty `linSpace` when size of `linSpace` is zero with `fopen = .true._LK, lopen = .true.`.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (0._TKG, 10._TKG)
            x2 = (10._TKG, 0._TKG)
            linSpace_ref = [(0._TKG, 10._TKG), (5._TKG, 5._TKG), (10._TKG, 0._TKG)]
#elif       RK_ENABLED
            x1 = 0._TKG
            x2 = 10._TKG
            linSpace_ref = [0._TKG, 5._TKG, 10._TKG]
#endif
            call report()
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(fopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `lopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .false._LK, lopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (0._TKG, 10._TKG)
            x2 = (10._TKG, 0._TKG)
            !linSpace_ref = [(2.5_TKG, 10_TKG), (5._TKG, 7.5_TKG), (7.5_TKG, 5._TKG), (10._TKG, 2.5_TKG)]
            linSpace_ref = [(2.5_TKG, 7.5_TKG), (5._TKG, 5._TKG), (7.5_TKG, 2.5_TKG), (10._TKG, 0._TKG)]
#elif       RK_ENABLED
            x1 = 0._TKG
            x2 = 10._TKG
            linSpace_ref = [2.5_TKG, 5._TKG, 7.5_TKG, 10._TKG]
#endif
            call report(fopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .true.` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(fopen = .true._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .true._LK, lopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (0._TKG, 10._TKG)
            x2 = (10._TKG, 0._TKG)
            linSpace_ref = [(0._TKG, 10._TKG), (2.5_TKG, 7.5_TKG), (5._TKG, 5._TKG), (7.5_TKG, 2.5_TKG)]
#elif       RK_ENABLED
            x1 = 0._TKG
            x2 = 10._TKG
            linSpace_ref = [0._TKG, 2.5_TKG, 5._TKG, 7.5_TKG]
#endif

            call report(lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `lopen = .true.` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(fopen = .false._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .false._LK, lopen = .true.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (0._TKG, 10._TKG)
            x2 = (10._TKG, 0._TKG)
            linSpace_ref = [(1.25_TKG, 8.75_TKG), (3.75_TKG, 6.25_TKG), (6.25_TKG, 3.75_TKG), (8.75_TKG, 1.25_TKG)]
#elif       RK_ENABLED
            x1 = 0._TKG
            x2 = 10._TKG
            linSpace_ref = [1.25_TKG, 3.75_TKG, 6.25_TKG, 8.75_TKG]
#endif
            call report(fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .true._LK, lopen = .true.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (-10._TKG, +10._TKG)
            x2 = (+10._TKG, -10._TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            x1 = -10._TKG
            x2 = +10._TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif
            call report(fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .true._LK, lopen = .true.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (-7.5_TKG, +7.5_TKG)
            x2 = (+7.5_TKG, -7.5_TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            x1 = -7.5_TKG
            x2 = +7.5_TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif
            call report()
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(fopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `lopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            call report(fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .false._LK, lopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (-7.5_TKG, +7.5_TKG)
            x2 = (+12.5_TKG, -12.5_TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            x1 = -7.5_TKG
            x2 = +12.5_TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif
            call report(fopen = .false._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .false._LK, lopen = .true.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (-12.5_TKG, +12.5_TKG)
            x2 = (+7.5_TKG, -7.5_TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            x1 = -12.5_TKG
            x2 = +7.5_TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif
            call report(fopen = .true._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .true._LK, lopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (-12.5_TKG, +12.5_TKG)
            x2 = (+7.5_TKG, -7.5_TKG)
            linSpace_ref = [x1]
#elif       RK_ENABLED
            x1 = -12.5_TKG
            x2 = +7.5_TKG
            linSpace_ref = [x1]
#endif
            call report(fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLinSpace() must return a `linSpace = x1` when `size(linSpace)==1` when x1, x2 = "//getStr([x1, x2]*sign))

            call report()
            call test%assert(assertion, SK_"getLinSpace() must return a `linSpace = x1` when `size(linSpace)==1` with `fopen = .false._LK, lopen = .false.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         CK_ENABLED
            x1 = (-10._TKG, +10._TKG)
            x2 = (+10._TKG, -10._TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            x1 = -10._TKG
            x2 = +10._TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif
            call report(fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLinSpace() must return an increasing `linSpace` with `fopen = .true._LK, lopen = .true.` when x1, x2 = "//getStr([x1, x2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        subroutine reset()
            if (allocated(linSpace)) deallocate(linSpace)
            if (allocated(linSpace_ref)) deallocate(linSpace_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(fopen, lopen)
            logical(LK), intent(in), optional :: fopen, lopen
#if         getLinSpace_ENABLED
            integer(IK) :: count
            count = size(linSpace_ref, kind = IK)
            linSpace = getLinSpace(x1 * sign, x2 * sign, count, fopen = fopen, lopen = lopen)
#elif       setLinSpace_ENABLED
            if (allocated(linSpace)) deallocate(linSpace)
            allocate(linSpace, mold = linSpace_ref)
            call setLinSpace(linSpace, x1 * sign, x2 * sign, fopen = fopen, lopen = lopen)
#else
#error      "Unrecognized interface."
#endif
            diff = abs(linSpace - linSpace_ref * sign)
            assertion = assertion .and. all(diff < tol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("x1")
                call test%disp%show( x1 )
                call test%disp%show("x2")
                call test%disp%show( x2 )
                call test%disp%show("linSpace_ref")
                call test%disp%show( linSpace_ref )
                call test%disp%show("linSpace")
                call test%disp%show( linSpace )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%show("TOL")
                call test%disp%show( TOL )
                call test%disp%show("sign")
                call test%disp%show( sign )
                if (present(fopen)) then
                    call test%disp%show("fopen")
                    call test%disp%show( fopen )
                end if
                if (present(lopen)) then
                    call test%disp%show("lopen")
                    call test%disp%show( lopen )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getLogSpace_ENABLED || setLogSpace_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_KIND :: logx1, logx2
        real(TKG), parameter :: TOL = epsilon(0._TKG) * 10 ! ipo requires the given eps precision.
        TYPE_KIND, allocatable :: logSpace(:)
        TYPE_KIND, allocatable :: linSpace_ref(:), logSpace_ref(:)
        TYPE_KIND, allocatable :: diff(:)
        integer(IK) :: sign

        assertion = .true._LK

        sign = 1_IK
        call testWith()

        sign = -1_IK
        call testWith()

        sign = 1_IK
        call testWith(base = 2._TKG)
        sign = -1_IK
        call testWith(base = 2._TKG)

    contains

        subroutine testWith(base)

            real(TKG), intent(in), optional :: base

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (0._TKG, 10._TKG)
            logx2 = (10._TKG, 0._TKG)
#elif       RK_ENABLED
            logx1 = 0._TKG
            logx2 = 10._TKG
#endif
            allocate(linSpace_ref(0))

            call report(base)
            call test%assert(assertion, SK_"getLogSpace() must return an empty `logSpace` when size of `logSpace` is zero.")

            call report(base, fopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an empty `logSpace` when size of `logSpace` is zero with `fopen = .false.`.")

            call report(base, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an empty `logSpace` when size of `logSpace` is zero with `lopen = .false.`.")

            call report(base, fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an empty `logSpace` when size of `logSpace` is zero with `fopen = .false._LK, lopen = .false.`.")

            call report(base, fopen = .false._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an empty `logSpace` when size of `logSpace` is zero with `fopen = .false._LK, lopen = .true.`.")

            call report(base, fopen = .true._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an empty `logSpace` when size of `logSpace` is zero with `fopen = .true._LK, lopen = .false.`.")

            call report(base, fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an empty `logSpace` when size of `logSpace` is zero with `fopen = .true._LK, lopen = .true.`.")

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (0._TKG, 10._TKG)
            logx2 = (10._TKG, 0._TKG)
            linSpace_ref = [(0._TKG, 10._TKG), (5._TKG, 5._TKG), (10._TKG, 0._TKG)]
#elif       RK_ENABLED
            logx1 = 0._TKG
            logx2 = 10._TKG
            linSpace_ref = [0._TKG, 5._TKG, 10._TKG]
#endif

            call report(base)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, fopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `lopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .false._LK, lopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (0._TKG, 10._TKG)
            logx2 = (10._TKG, 0._TKG)
            !linSpace_ref = [(2.5_TKG, 10_TKG), (5._TKG, 7.5_TKG), (7.5_TKG, 5._TKG), (10._TKG, 2.5_TKG)]
            linSpace_ref = [(2.5_TKG, 7.5_TKG), (5._TKG, 5._TKG), (7.5_TKG, 2.5_TKG), (10._TKG, 0._TKG)]
#elif       RK_ENABLED
            logx1 = 0._TKG
            logx2 = 10._TKG
            linSpace_ref = [2.5_TKG, 5._TKG, 7.5_TKG, 10._TKG]
#endif

            call report(base, fopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .true.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, fopen = .true._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .true._LK, lopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (0._TKG, 10._TKG)
            logx2 = (10._TKG, 0._TKG)
            linSpace_ref = [(0._TKG, 10._TKG), (2.5_TKG, 7.5_TKG), (5._TKG, 5._TKG), (7.5_TKG, 2.5_TKG)]
#elif       RK_ENABLED
            logx1 = 0._TKG
            logx2 = 10._TKG
            linSpace_ref = [0._TKG, 2.5_TKG, 5._TKG, 7.5_TKG]
#endif

            call report(base, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `lopen = .true.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, fopen = .false._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .false._LK, lopen = .true.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (0._TKG, 10._TKG)
            logx2 = (10._TKG, 0._TKG)
            linSpace_ref = [(1.25_TKG, 8.75_TKG), (3.75_TKG, 6.25_TKG), (6.25_TKG, 3.75_TKG), (8.75_TKG, 1.25_TKG)]
#elif       RK_ENABLED
            logx1 = 0._TKG
            logx2 = 10._TKG
            linSpace_ref = [1.25_TKG, 3.75_TKG, 6.25_TKG, 8.75_TKG]
#endif

            call report(base, fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .true._LK, lopen = .true.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (-10._TKG, +10._TKG)
            logx2 = (+10._TKG, -10._TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            logx1 = -10._TKG
            logx2 = +10._TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif

            call report(base, fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .true._LK, lopen = .true.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (-7.5_TKG, +7.5_TKG)
            logx2 = (+7.5_TKG, -7.5_TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            logx1 = -7.5_TKG
            logx2 = +7.5_TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif

            call report(base)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, fopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `lopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base, fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .false._LK, lopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (-7.5_TKG, +7.5_TKG)
            logx2 = (+12.5_TKG, -12.5_TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            logx1 = -7.5_TKG
            logx2 = +12.5_TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif

            call report(base, fopen = .false._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .false._LK, lopen = .true.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (-12.5_TKG, +12.5_TKG)
            logx2 = (+7.5_TKG, -7.5_TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            logx1 = -12.5_TKG
            logx2 = +7.5_TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif

            call report(base, fopen = .true._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .true._LK, lopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (-12.5_TKG, +12.5_TKG)
            logx2 = (+7.5_TKG, -7.5_TKG)
            linSpace_ref = [logx1]
#elif       RK_ENABLED
            logx1 = -12.5_TKG
            logx2 = +7.5_TKG
            linSpace_ref = [logx1]
#endif

            call report(base, fopen = .false._LK, lopen = .false._LK)
            call test%assert(assertion, SK_"getLogSpace() must return a `logSpace = logx1` when `size(logSpace)==1` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            call report(base)
            call test%assert(assertion, SK_"getLogSpace() must return a `logSpace = logx1` when `size(logSpace)==1` with `fopen = .false._LK, lopen = .false.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         CK_ENABLED
            logx1 = (-10._TKG, +10._TKG)
            logx2 = (+10._TKG, -10._TKG)
            linSpace_ref = [(-7.5_TKG, 7.5_TKG), (-2.5_TKG, 2.5_TKG), (2.5_TKG, -2.5_TKG), (7.5_TKG, -7.5_TKG)]
#elif       RK_ENABLED
            logx1 = -10._TKG
            logx2 = +10._TKG
            linSpace_ref = [-7.5_TKG, -2.5_TKG, 2.5_TKG, 7.5_TKG]
#endif

            call report(base, fopen = .true._LK, lopen = .true._LK)
            call test%assert(assertion, SK_"getLogSpace() must return an increasing `logSpace` with `fopen = .true._LK, lopen = .true.` when logx1, logx2 = "//getStr([logx1, logx2]*sign))

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(logSpace)) deallocate(logSpace)
            if (allocated(linSpace_ref)) deallocate(linSpace_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(base, fopen, lopen)
            logical(LK), intent(in), optional :: fopen, lopen
            real(TKG), intent(in), optional :: base
#if         getLogSpace_ENABLED
            integer(IK) :: count
            count = size(linSpace_ref, kind = IK)
            logSpace = getLogSpace(logx1 * sign, logx2 * sign, count, fopen = fopen, lopen = lopen, base = base)
#elif       setLogSpace_ENABLED
            if (allocated(logSpace)) deallocate(logSpace)
            allocate(logSpace, mold = linSpace_ref)
            call setLogSpace(logSpace, logx1 * sign, logx2 * sign, fopen = fopen, lopen = lopen, base = base)
#else
#error      "Unrecognized interface."
#endif
            if (present(base)) then
                logSpace_ref = base**(linSpace_ref * sign)
            else
                logSpace_ref = exp(linSpace_ref * sign)
            end if
            if (allocated(diff)) deallocate(diff)
            allocate(diff, mold = logSpace_ref)
#if         CK_ENABLED
            where (logSpace_ref%re > 0._TKG)
                diff%re = abs(logSpace%re - logSpace_ref%re) / logSpace_ref%re
            elsewhere
                diff%re = abs(logSpace%re - logSpace_ref%re)
            end where
            where (logSpace_ref%im > 0._TKG)
                diff%im = abs(logSpace%im - logSpace_ref%im) / logSpace_ref%im
            elsewhere
                diff%im = abs(logSpace%im - logSpace_ref%im)
            end where
            assertion = assertion .and. all(diff%re < TOL)
            assertion = assertion .and. all(diff%im < TOL)
#elif       RK_ENABLED
            where (logSpace_ref > 0._TKG)
                diff = abs(logSpace - logSpace_ref) / logSpace_ref
            elsewhere
                diff = abs(logSpace - logSpace_ref)
            end where
            assertion = assertion .and. all(diff < TOL)
#endif
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("logx1")
                call test%disp%show( logx1 )
                call test%disp%show("logx2")
                call test%disp%show( logx2 )
                call test%disp%show("linSpace_ref")
                call test%disp%show( linSpace_ref )
                call test%disp%show("logSpace_ref")
                call test%disp%show( logSpace_ref )
                call test%disp%show("logSpace")
                call test%disp%show( logSpace )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%show("TOL")
                call test%disp%show( TOL )
                call test%disp%show("sign")
                call test%disp%show( sign )
                if (present(base)) then
                    call test%disp%show("base")
                    call test%disp%show( base )
                end if
                if (present(fopen)) then
                    call test%disp%show("fopen")
                    call test%disp%show( fopen )
                end if
                if (present(lopen)) then
                    call test%disp%show("lopen")
                    call test%disp%show( lopen )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_KIND