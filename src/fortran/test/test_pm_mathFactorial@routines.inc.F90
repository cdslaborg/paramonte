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

        integer(IKG), allocatable   :: PosInt(:), factorial(:), factorial_ref(:), diff(:)

        assertion = .true._LK

        PosInt = [integer(IKG)::]
        factorial_ref = [integer(IKG)::]
        factorial = getFactorial(PosInt)

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield an empty `factorial` with an empty input `n`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [0_IKG]
        factorial_ref = [1_IKG]
        factorial = [getFactorial(PosInt(1))]

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield 1 for `n = 0`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [1_IKG]
        factorial_ref = [1_IKG]
        factorial = [getFactorial(PosInt(1))]

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield 1 for `n = 1`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [1_IKG, 2_IKG]
        factorial_ref = [1_IKG, 2_IKG]
        factorial = getFactorial(PosInt)

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield `[1_IKG, 2_IKG]` for `n = [1_IKG, 2_IKG]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PosInt = [5_IKG]
        factorial_ref = [120_IKG]
        factorial = getFactorial(PosInt)

        call report()
        call test%assert(assertion, SK_"getFactorial() must yield `[120_IKG]` for `n = [5_IKG]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(factorial - factorial_ref)
            assertion = assertion .and. all(diff == 0_IKG)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "PosInt         ", PosInt
                write(test%disp%unit,"(*(g0,:,', '))") "factorial      ", factorial
                write(test%disp%unit,"(*(g0,:,', '))") "factorial_ref  ", factorial_ref
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#elif   getLogFactorial_ENABLED

        real(RK), parameter :: EPS = epsilon(0._RK) * 10
        real(RK), allocatable :: wholeNumber(:), factorial(:), factorial_ref(:), diff(:)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wholeNumber = [real(RK)::]
        factorial_ref = log_gamma([real(RK)::])
        factorial = getLogFactorial(wholeNumber)

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield an empty `factorial` with an empty input `n`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wholeNumber = [0._RK]
        factorial_ref = log_gamma(wholeNumber + 1._RK)
        factorial = [getLogFactorial(wholeNumber(1))]

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = 0`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wholeNumber = [1._RK]
        factorial_ref = log_gamma([1._RK])
        factorial = [getLogFactorial(wholeNumber(1))]

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = 1`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wholeNumber = [1._RK, 2._RK]
        factorial_ref = log_gamma(wholeNumber + 1._RK)
        factorial = getLogFactorial(wholeNumber)

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = [1._RK, 2._RK]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        wholeNumber = [5._RK]
        factorial_ref = log_gamma(wholeNumber + 1._RK)
        factorial = getLogFactorial(wholeNumber)

        call report()
        call test%assert(assertion, SK_"getLogFactorial() must yield correct answer for `x = [5._RK]`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            diff = abs(factorial - factorial_ref)
            assertion = assertion .and. all(abs(diff) <= EPS)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "wholeNumber    ", wholeNumber
                write(test%disp%unit,"(*(g0,:,', '))") "factorial      ", factorial
                write(test%disp%unit,"(*(g0,:,', '))") "factorial_ref  ", factorial_ref
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
#error  "Unrecognized interface."
#endif
