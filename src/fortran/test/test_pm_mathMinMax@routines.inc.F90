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
!>  This include file contains the implementations of the tests of procedures in [pm_mathMinMax](@ref pm_mathMinMax).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMinMax_SK_ENABLED || setMinMax_SK_ENABLED
        character(:, SK), allocatable   :: MinMax_ref(:)
        character(:, SK), allocatable   :: minMax(:)
        character(:, SK), allocatable   :: smaller
        character(:, SK), allocatable   :: larger
#elif   getMinMax_SSK_ENABLED || setMinMax_SSK_ENABLED
        use pm_container, only: css_pdt
        type(css_pdt)   :: MinMax_ref(2)
        type(css_pdt)   :: minMax(2)
        type(css_pdt)   :: smaller
        type(css_pdt)   :: larger
#elif   getMinMax_CK_ENABLED || setMinMax_CK_ENABLED
        complex(CK)     :: MinMax_ref(2)
        complex(CK)     :: minMax(2)
        complex(CK)     :: smaller
        complex(CK)     :: larger
#elif   getMinMax_IK_ENABLED || setMinMax_IK_ENABLED
        integer(IK)     :: MinMax_ref(2)
        integer(IK)     :: minMax(2)
        integer(IK)     :: smaller
        integer(IK)     :: larger
#elif   getMinMax_RK_ENABLED || setMinMax_RK_ENABLED
        real(RK)        :: MinMax_ref(2)
        real(RK)        :: minMax(2)
        real(RK)        :: smaller
        real(RK)        :: larger
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setValues()
#if     getMinMax_ENABLED
        minMax = getMinMax(smaller, larger)
#elif   setMinMax_ENABLED
        call setMinMax(smaller, larger)
        minMax = [smaller, larger]

#else
#error  "Unrecognized interface."
#endif
        call report()
        call test%assert(assertion, SK_"The procedure must correctly output min and max when `a < b`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setValues()
#if     getMinMax_ENABLED
        minMax = getMinMax(larger, smaller)
#elif   setMinMax_ENABLED
        call setMinMax(larger, smaller)
        minMax = [larger, smaller]
#else
#error  "Unrecognized interface."
#endif
        call report()
        call test%assert(assertion, SK_"The procedure must correctly output min and max when `a > b`.")

#if     getMinMax_ENABLED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setValues()
        minMax = getMinMax([smaller, larger])

        call report()
        call test%assert(assertion, SK_"The procedure must correctly output min and max when `Pair(1) < Pair(2)`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call setValues()
        minMax = getMinMax([larger, smaller])

        call report()
        call test%assert(assertion, SK_"The procedure must correctly output min and max when `Pair(1) > Pair(2)`.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setValues()
#if         getMinMax_SSK_ENABLED || setMinMax_SSK_ENABLED
            use pm_distUnif, only: getUnifRand
            smaller%val = getUnifRand(SK_"AAA", SK_"HHH")
            larger%val = getUnifRand(SK_"OOO", SK_"ZZZ")
#elif       getMinMax_SK_ENABLED || setMinMax_SK_ENABLED
            use pm_distUnif, only: getUnifRand
            smaller = getUnifRand(SK_"AAA", SK_"HHH")
            larger = getUnifRand(SK_"OOO", SK_"ZZZ")
#else
            use pm_distUnif, only: setUnifRand
            call setUnifRand(smaller)
            call setUnifRand(smaller)
            call setUnifRand(larger)
            call setUnifRand(larger)
            smaller = -smaller
#endif
            MinMax_ref = [smaller, larger]
        end subroutine setValues

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()

#if         getMinMax_SSK_ENABLED || setMinMax_SSK_ENABLED
            use pm_container, only: operator(==)
            integer(IK) :: i
#endif
            assertion = assertion .and. all(minMax == MinMax_ref)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
#if             getMinMax_SSK_ENABLED || setMinMax_SSK_ENABLED
                write(test%disp%unit,"(*(g0,:,', '))") "MinMax_ref ", (MinMax_ref(i)%val, i = 1, 2)
                write(test%disp%unit,"(*(g0,:,', '))") "minMax     ", (minMax(i)%val, i = 1, 2)
#else
                write(test%disp%unit,"(*(g0,:,', '))") "MinMax_ref ", MinMax_ref
                write(test%disp%unit,"(*(g0,:,', '))") "minMax     ", minMax
               !write(test%disp%unit,"(*(g0,:,', '))") "diff       ", abs(minMax - MinMax_ref)
#endif
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
