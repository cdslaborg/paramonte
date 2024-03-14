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

!>  \brief This file contains the implementations of the tests of module [pm_logicalCompare](@ref pm_logicalCompare).
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_kind, only: IK, LK

        integer(IK) :: i
        logical(LK) :: array1(4), array2(4), Result(4), Result_def(4)

        assertion = .true._LK

        array1 = [.true._LK, .true._LK, .false._LK, .false._LK]
        array2 = [.true._LK, .false._LK, .true._LK, .false._LK]

#if     test_isless_LK_ENABLED
#define COMPARES_WITH <
        Result_def = [.false._LK, .false._LK, .true._LK, .false._LK]
#elif   test_isleq_LK_ENABLED
#define COMPARES_WITH <=
        Result_def = [.true._LK, .false._LK, .true._LK, .true._LK]
#elif   test_isneq_LK_ENABLED
#define COMPARES_WITH /=
        Result_def = [.false._LK, .true._LK, .true._LK, .false._LK]
#elif   test_iseq_LK_ENABLED
#define COMPARES_WITH ==
        Result_def = [.true._LK, .false._LK, .false._LK, .true._LK]
#elif   test_ismeq_LK_ENABLED
#define COMPARES_WITH >=
        Result_def = [.true._LK, .true._LK, .false._LK, .true._LK]
#elif   test_ismore_LK_ENABLED
#define COMPARES_WITH >
        Result_def = [.false._LK, .true._LK, .false._LK, .false._LK]
#else
#error  "Unrecognized interface."
#endif

        do i = 1, size(array1)
            Result(i) = array1(i) COMPARES_WITH array2(i)
        end do
        assertion = assertion .and. all(Result .eqv. Result_def)
        call report()
        call test%assert(assertion, SK_"Two scalar logical values must be compared correctly.")

        Result = array1 COMPARES_WITH array2
        assertion = assertion .and. all(Result .eqv. Result_def)
        call report()
        call test%assert(assertion, SK_"Two vector logical values must be compared correctly.")

    contains

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "array1     ", array1
                write(test%disp%unit,"(*(g0,:,', '))") "array2     ", array2
                write(test%disp%unit,"(*(g0,:,', '))") "Result     ", Result
                write(test%disp%unit,"(*(g0,:,', '))") "Result_def ", Result_def
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#undef  COMPARES_WITH