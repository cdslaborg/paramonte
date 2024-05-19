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
!>  This module contains implementations of the tests of the procedures under the generic interface [pm_option](@ref pm_option).
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if option_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

        use pm_val2str, only: getStr

#if option_CK_ENABLED
        complex(CK)     , allocatable   :: Optional(:,:), Default(:,:)
#elif option_RK_ENABLED
        real(RK)        , allocatable   :: Optional(:,:), Default(:,:)
#elif option_IK_ENABLED
        integer(IKG)    , allocatable   :: Optional(:,:), Default(:,:)
#elif option_LK_ENABLED
        logical(LK)     , allocatable   :: Optional(:,:), Default(:,:)
#elif option_SK_ENABLED
        character(3, SK), allocatable   :: Optional(:,:), Default(:,:)
#else
#error "Unrecognized interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

        allocate(Optional(0,0), Default(0,0))

        assertion = assertion .and. all(Optional IS_EQUAL getOption(Default, Optional))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of size zero when Optional is present and has size zero.")

        assertion = assertion .and. all(Default IS_EQUAL getOption(Default))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of size zero when Optional is missing and `Default` has size zero.")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if option_CK_ENABLED
        Default     = reshape([(1._CK, -1._CK)], shape = [1,1])
        Optional    = reshape([(-1._CK, 1._CK)], shape = [1,1])
#elif option_RK_ENABLED
        Default     = reshape([+1._RK], shape = [1,1])
        Optional    = reshape([-1._RK], shape = [1,1])
#elif option_IK_ENABLED
        Default     = reshape([+1_IKG], shape = [1,1])
        Optional    = reshape([-1_IKG], shape = [1,1])
#elif option_LK_ENABLED
        Default     = reshape([.true._LK], shape = [1,1])
        Optional    = reshape([.false._LK], shape = [1,1])
#elif option_SK_ENABLED
        Default     = reshape(["ABC"], shape = [1,1])
        Optional    = reshape(["CBA"], shape = [1,1])
#endif

        assertion = assertion .and. all(Optional IS_EQUAL getOption(Default, Optional))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of shape [1,1] one when Optional is present and has shape [1,1].")

        assertion = assertion .and. all(Optional(:,1) IS_EQUAL getOption(Default(:,1), Optional(:,1)))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of shape [1] when Optional is present and is given as [:,1].")

        assertion = assertion .and. all(Optional(1,:) IS_EQUAL getOption(Default(1,:), Optional(1,:)))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of shape [1] when Optional is present and is given as [1,:].")

        assertion = assertion .and. Optional(1,1) IS_EQUAL getOption(Default(1,1), Optional(1,1))
        call report()
        call test%assert(assertion, desc = "getOption() must return a scalar object when Optional is present and is given as [1,1].")

        assertion = assertion .and. all(Default IS_EQUAL getOption(Default))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of size one when Optional is missing and `Default` has size of one.")
        
        assertion = assertion .and. all(Default(:,1) IS_EQUAL getOption(Default(:,1)))
        call report()
        call test%assert(assertion, desc = "getOption() must return a 1D `Default` when Optional is missing and `Default` has size of one and given as `Default(:,1)`.")

        assertion = assertion .and. all(Default(1,:) IS_EQUAL getOption(Default(1,:)))
        call report()
        call test%assert(assertion, desc = "getOption() must return a 1D `Default` when Optional is missing and `Default` has size of one and given as a non-contiguous `Default(1,:)`.")

        assertion = assertion .and. Default(1,1) IS_EQUAL getOption(Default(1,1))
        call report()
        call test%assert(assertion, desc = "getOption() must return a 0D `Default(1,1)` when Optional is missing and `Default` has size of one and given as a non-contiguous `Default(1,1)`.")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if option_CK_ENABLED
        Default     = reshape([(+1._CK, -1._CK), (+2._CK, -2._CK), (+3._CK, -3._CK), (+4._CK, -4._CK)], shape = [2,2])
        Optional    = reshape([(+1._CK, -1._CK), (-2._CK, +2._CK), (-3._CK, +3._CK), (-4._CK, +4._CK)], shape = [2,2])
#elif option_RK_ENABLED
        Default     = reshape([+1._RK, +2._RK, +3._RK, +4._RK], shape = [2,2])
        Optional    = reshape([-1._RK, -2._RK, -3._RK, -4._RK], shape = [2,2])
#elif option_IK_ENABLED
        Default     = reshape([+1_IKG, +2_IKG, +3_IKG, +4_IKG], shape = [2,2])
        Optional    = reshape([-1_IKG, -2_IKG, -3_IKG, -4_IKG], shape = [2,2])
#elif option_LK_ENABLED
        Default     = reshape([.true._LK, .false._LK, .true._LK, .false._LK], shape = [2,2])
        Optional    = reshape([.false._LK, .true._LK, .false._LK, .true._LK], shape = [2,2])
#elif option_SK_ENABLED
        Default     = reshape(["ABC", "DEF", "GHI", "JKL"], shape = [2,2])
        Optional    = reshape(["CBA", "FED", "IHG", "LKJ"], shape = [2,2])
#endif

        assertion = assertion .and. all(Optional IS_EQUAL getOption(Default, Optional))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of shape [2,2] one when Optional is present and has shape [2,2].")

        assertion = assertion .and. all(Optional(:,1) IS_EQUAL getOption(Default(:,1), Optional(:,1)))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of shape [1] when Optional is present and is given as [:,1].")

        assertion = assertion .and. all(Optional(1,:) IS_EQUAL getOption(Default(1,:), Optional(1,:)))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of shape [1] when Optional is present and is given as [1,:].")

        assertion = assertion .and. Optional(2,2) IS_EQUAL getOption(Default(2,2), Optional(2,2))
        call report()
        call test%assert(assertion, desc = "getOption() must return a scalar object when Optional is present and is given as [2,2].")

        assertion = assertion .and. all(Default IS_EQUAL getOption(Default))
        call report()
        call test%assert(assertion, desc = "getOption() must return an object of size one when Optional is missing and `Default` has size of one.")
        
        assertion = assertion .and. all(Default(:,1) IS_EQUAL getOption(Default(:,1)))
        call report()
        call test%assert(assertion, desc = "getOption() must return a 1D `Default` when Optional is missing and `Default` has size of one and given as `Default(:,1)`.")

        assertion = assertion .and. all(Default(1,:) IS_EQUAL getOption(Default(1,:)))
        call report()
        call test%assert(assertion, desc = "getOption() must return a 1D `Default` when Optional is missing and `Default` has size of one and given as a non-contiguous `Default(1,:)`.")

        assertion = assertion .and. Default(2,2) IS_EQUAL getOption(Default(2,2))
        call report()
        call test%assert(assertion, desc = "getOption() must return a 0D `Default(2,2)` when Optional is missing and `Default` has size of one and given as a non-contiguous `Default(2,2)`.")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(Optional)) deallocate(Optional)
            if (allocated(Default)) deallocate(Default)
        end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Optional   ", Optional
                write(test%disp%unit,"(*(g0,:,', '))") "Default    ", Default
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef IS_EQUAL
