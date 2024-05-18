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
!>  This module contains implementations of the tests of the module [pm_mathFactorial](@ref pm_mathFactorial).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_kind, only: LK
        use pm_val2str, only: getStr

#if     getFactorial_ENABLED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IKC), allocatable   :: PosInt(:), Factorial(:), Factorial_ref(:), diff(:)

        assertion = .true._LK

        PosInt = [integer(IKC)::]
        Factorial_ref = [integer(IKC)::]
        Factorial = getFactorial(PosInt)

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield an empty `Factorial` with an empty input `n`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [0_IKC]
        Factorial_ref = [1_IKC]
        Factorial = [getFactorial(PosInt(1))]

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield 1 for `n = 0`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [1_IKC]
        Factorial_ref = [1_IKC]
        Factorial = [getFactorial(PosInt(1))]

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield 1 for `n = 1`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [1_IKC, 2_IKC]
        Factorial_ref = [1_IKC, 2_IKC]
        Factorial = getFactorial(PosInt)

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield `[1_IKC, 2_IKC]` for `n = [1_IKC, 2_IKC]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [5_IKC]
        Factorial_ref = [120_IKC]
        Factorial = getFactorial(PosInt)

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield `[120_IKC]` for `n = [5_IKC]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(Factorial - Factorial_ref)
            assertion = assertion .and. all(diff == 0_IKC)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "PosInt         ", PosInt
                write(test%disp%unit,"(*(g0,:,', '))") "Factorial      ", Factorial
                write(test%disp%unit,"(*(g0,:,', '))") "Factorial_ref  ", Factorial_ref
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#elif   getLogFactorial_ENABLED

        real(RK), parameter :: EPS = epsilon(0._RK) * 10
        real(RK), allocatable :: WholeNumber(:), Factorial(:), Factorial_ref(:), diff(:)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        WholeNumber = [real(RK)::]
        Factorial_ref = log_gamma([real(RK)::])
        Factorial = getLogFactorial(WholeNumber)

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield an empty `Factorial` with an empty input `n`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        WholeNumber = [0._RK]
        Factorial_ref = log_gamma(WholeNumber + 1._RK)
        Factorial = [getLogFactorial(WholeNumber(1))]

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = 0`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        WholeNumber = [1._RK]
        Factorial_ref = log_gamma([1._RK])
        Factorial = [getLogFactorial(WholeNumber(1))]

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = 1`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        WholeNumber = [1._RK, 2._RK]
        Factorial_ref = log_gamma(WholeNumber + 1._RK)
        Factorial = getLogFactorial(WholeNumber)

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = [1._RK, 2._RK]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        WholeNumber = [5._RK]
        Factorial_ref = log_gamma(WholeNumber + 1._RK)
        Factorial = getLogFactorial(WholeNumber)

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = [5._RK]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(Factorial - Factorial_ref)
            assertion = assertion .and. all(abs(diff) <= EPS)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "WholeNumber    ", WholeNumber
                write(test%disp%unit,"(*(g0,:,', '))") "Factorial      ", Factorial
                write(test%disp%unit,"(*(g0,:,', '))") "Factorial_ref  ", Factorial_ref
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
#error  "Unrecognized interface."
#endif
