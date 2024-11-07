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
!>  This include file contains procedure implementations of the tests of [val2complex_pmod](@ref pm_val2complex).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_getComplex_LK_ENABLED || test_getComplex32_LK_ENABLED || test_getComplex64_LK_ENABLED || test_getComplex128_LK_ENABLED
#define test_LK_ENABLED 1
        logical(LK)     , allocatable   :: val(:)
#elif   test_getComplex_SK_ENABLED || test_getComplex32_SK_ENABLED || test_getComplex64_SK_ENABLED || test_getComplex128_SK_ENABLED
#define test_SK_ENABLED 1
        character(20,SK), allocatable   :: val(:)
        integer(IK)     , allocatable   :: iostat(:)
#else
#error  "Unrecognized interface."
#endif
        complex(CK)     , allocatable   :: Conversion(:), Conversion_ref(:)

        integer(IK) :: i

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     test_SK_ENABLED
        val = [character(20,SK) ::]
        allocate(iostat(0))
#elif   test_LK_ENABLED
        val = [logical(LK) ::]
#endif
        allocate(Conversion_ref(0))

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion == Conversion_ref)
        call report()
        call test%assert(assertion, desc = "An empty val array must yield an empty Conversion array")

#if     test_SK_ENABLED
        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(Conversion == Conversion_ref)
        call report()
        call test%assert(assertion, desc = "An empty val array must yield an empty Conversion array with `iostat` present.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     test_SK_ENABLED
        val = [character(20,SK) :: "1.", "12.", "123."]
        Conversion_ref = [(1._CK,0._CK), (12._CK,0._CK), (123._CK,0._CK)]
        allocate(iostat(size(Conversion_ref)))
#elif   test_LK_ENABLED
        val = [logical(LK) :: .false., .true., .false._LK]
        Conversion_ref = [(0._CK, 0._CK), (1._CK, 1._CK), (0._CK, 0._CK)]
#endif
        allocate(Conversion(size(Conversion_ref)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i))
            assertion = assertion .and. (Conversion(i) == Conversion_ref(i))
            call report()
            call test%assert(assertion, desc = "An valid scalar value for val must yield a valid output value.")
#if         test_SK_ENABLED
            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) == 0_IK)
            call report()
            call test%assert(assertion, desc = "An valid scalar value for val must yield a valid output value with `iostat = 0`.")
#endif
        end do

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion == Conversion_ref)
        call report()
        call test%assert(assertion, desc = "An non-empty val array of legal values must yield a Conversion array of the same size without errors.")

#if     test_SK_ENABLED
        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat == 0_IK)
        call report()
        call test%assert(assertion, desc = "An non-empty val array of legal values must yield a Conversion array of the same size with `all(iostat == 0)`.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     test_SK_ENABLED
        val = [character(20,SK) :: "1., paramonte", "12., paramonte", "123., paramonte"]
        Conversion_ref = [(1._CK, 0._CK), (12._CK, 0._CK), (123._CK, 0._CK)]
        allocate(iostat(size(val)))
#elif   test_LK_ENABLED
        val = [logical(LK) :: .false., .true., .false., .true., .false._LK]
        Conversion_ref = [(0._CK, 0._CK), (1._CK, 1._CK), (0._CK, 0._CK), (1._CK, 1._CK), (0._CK, 0._CK)]
#endif

        allocate(Conversion(size(val)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i))
            assertion = assertion .and. (Conversion(i) == Conversion_ref(i))
            call report()
            call test%assert(assertion, desc = "An valid scalar value for val must yield a valid output value even in the presence subsequent non-sense values.")
#if         test_SK_ENABLED
            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) == 0_IK)
            call report()
            call test%assert(assertion, desc = "An valid scalar value for val must yield a valid output value with `iostat = 0` even in the presence subsequent non-sense values.")
#endif
        end do

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion == Conversion_ref)
        call report()
        call test%assert(assertion, desc = "An non-empty val array of legal values must yield a Conversion array of the same size without errors even in the presence subsequent non-sense values.")

#if     test_SK_ENABLED
        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat == 0_IK)
        call report()
        call test%assert(assertion, desc = "An non-empty val array of legal values must yield a Conversion array of the same size with `all(iostat == 0_IK)` even in the presence subsequent non-sense values.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_SK_ENABLED
        call reset()

        val = [character(20,SK) :: "x", "y"]
        allocate(iostat(size(val)), Conversion(size(val)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) /=  0_IK)
            call report()
            call test%assert(assertion, desc = "An invalid scalar value for val must yield an invalid output value with `iostat /= 0`.")
        end do

        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat /=  0_IK)
        call report()
        call test%assert(assertion, desc = "An invalid vector val must yield an invalid vector output with `iostat /= 0`.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if test_SK_ENABLED
        call reset()
        val = [character(20,SK) :: "1. -1.", "12. -12.", "123. -123."]
        Conversion_ref = cmplx([(1.,-1.), (12.,-12.), (123.,-123.)], kind = CK)
        allocate(iostat(size(val)), Conversion(size(val)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i))
            assertion = assertion .and. (Conversion(i)%re == Conversion_ref(i)%re)
            call report()
            call test%assert(assertion, desc = "An valid space-separated complex scalar value for val must yield a valid output value.")

            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) == 0_IK)
            call report()
            call test%assert(assertion, desc = "An valid space-separated complex scalar value for val must yield a valid output value with `iostat = 0`.")
        end do

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion%re == Conversion_ref%re)
        call report()
        call test%assert(assertion, desc = "An non-empty space-separated complex val array of legal values must yield a val array of the same size without errors.")

        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat == 0_IK)
        call report()
        call test%assert(assertion, desc = "An non-empty space-separated complex val array of legal values must yield a val array of the same size with `all(iostat == 0_IK)`.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if test_SK_ENABLED
        call reset()
        val = [character(20,SK) :: "1., -1.", "12., -12.", "123., -123."]
        Conversion_ref = cmplx([(1.,-1.), (12.,-12.), (123.,-123.)], kind = CK)
        allocate(iostat(size(val)), Conversion(size(val)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i))
            assertion = assertion .and. (Conversion(i)%re == Conversion_ref(i)%re)
            call report()
            call test%assert(assertion, desc = "An valid comma-separated complex scalar value for val must yield a valid output value.")

            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) == 0_IK)
            call report()
            call test%assert(assertion, desc = "An valid comma-separated complex scalar value for val must yield a valid output value with `iostat = 0`.")
        end do

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion%re == Conversion_ref%re)
        call report()
        call test%assert(assertion, desc = "An non-empty comma-separated complex val array of legal values must yield a val array of the same size without errors.")

        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat == 0_IK)
        call report()
        call test%assert(assertion, desc = "An non-empty comma-separated complex val array of legal values must yield a val array of the same size with `all(iostat == 0_IK)`.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if test_SK_ENABLED
        call reset()
        val = [character(20,SK) :: "1., x", "12., x", "123., x"]
        Conversion_ref = cmplx([(1.,-1.), (12.,-12.), (123.,-123.)], kind = CK)
        allocate(iostat(size(val)), Conversion(size(val)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i))
            assertion = assertion .and. (Conversion(i)%re == Conversion_ref(i)%re)
            call report()
            call test%assert(assertion, desc = "An valid comma-separated complex scalar value for val with invalid imaginary component must yield a valid output value.")

            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) == 0_IK)
            call report()
            call test%assert(assertion, desc = "An valid comma-separated complex scalar value for val with invalid imaginary component must yield a valid output value with `iostat = 0`.")
        end do

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion%re == Conversion_ref%re)
        call report()
        call test%assert(assertion, desc = "An non-empty comma-separated complex val array of legal values with invalid imaginary components must yield a val array of the same size without errors.")

        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat == 0_IK)
        call report()
        call test%assert(assertion, desc = "An non-empty comma-separated complex val array of legal values with invalid imaginary components must yield a val array of the same size with `all(iostat == 0_IK)`.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if test_SK_ENABLED
        call reset()
        val = [character(20,SK) :: "(1., -1.)", "(12., -12.)", "(123., -123.)"]
        Conversion_ref = cmplx([(1.,-1.), (12.,-12.), (123.,-123.)], kind = CK)
        allocate(iostat(size(val)), Conversion(size(val)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i))
            assertion = assertion .and. (Conversion(i) == Conversion_ref(i))
            call report()
            call test%assert(assertion, desc = "An valid parenthesized comma-separated complex scalar value for val must yield a valid output value.")

            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) == 0_IK)
            call report()
            call test%assert(assertion, desc = "An valid parenthesized comma-separated complex scalar value for val must yield a valid output value with `iostat = 0`.")
        end do

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion == Conversion_ref)
        call report()
        call test%assert(assertion, desc = "An non-empty parenthesized comma-separated complex val array of legal values must yield a val array of the same size without errors.")

        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat == 0_IK)
        call report()
        call test%assert(assertion, desc = "An non-empty parenthesized comma-separated complex val array of legal values must yield a val array of the same size with `all(iostat == 0_IK)`.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if test_SK_ENABLED
        call reset()
        val = [character(20,SK) :: "1. x", "12. x", "123. x"]
        Conversion_ref = cmplx([(1.,-1.), (12.,-12.), (123.,-123.)], kind = CK)
        allocate(iostat(size(val)), Conversion(size(val)))

        do i = 1, size(val)
            Conversion(i) = getComplex(val(i))
            assertion = assertion .and. (Conversion(i)%re == Conversion_ref(i)%re)
            call report()
            call test%assert(assertion, desc = "An valid space-separated complex scalar value for val with invalid imaginary component must yield a valid output value.")

            Conversion(i) = getComplex(val(i), iostat(i))
            assertion = assertion .and. (iostat(i) == 0_IK)
            call report()
            call test%assert(assertion, desc = "An valid space-separated complex scalar value for val with invalid imaginary component must yield a valid output value with `iostat = 0`.")
        end do

        Conversion = getComplex(val)
        assertion = assertion .and. all(Conversion%re == Conversion_ref%re)
        call report()
        call test%assert(assertion, desc = "An non-empty space-separated complex val array of legal values with invalid imaginary components must yield a val array of the same size without errors.")

        Conversion = getComplex(val, iostat)
        assertion = assertion .and. all(iostat == 0_IK)
        call report()
        call test%assert(assertion, desc = "An non-empty space-separated complex val array of legal values with invalid imaginary components must yield a val array of the same size with `all(iostat == 0_IK)`.")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(val)) deallocate(val)
            if (allocated(Conversion)) deallocate(Conversion)
            if (allocated(Conversion_ref)) deallocate(Conversion_ref)
#if         test_SK_ENABLED
            if (allocated(iostat)) deallocate(iostat)
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "val            ", val
                write(test%disp%unit,"(*(g0,:,', '))") "Conversion_ref ", Conversion_ref
                write(test%disp%unit,"(*(g0,:,', '))") "Conversion     ", Conversion
#if     test_SK_ENABLED
                write(test%disp%unit,"(*(g0,:,', '))") "iostat         ", iostat
#endif
                write(test%disp%unit,"(*(g0,:,', '))") "i              ", i
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  test_LK_ENABLED
#undef  test_SK_ENABLED