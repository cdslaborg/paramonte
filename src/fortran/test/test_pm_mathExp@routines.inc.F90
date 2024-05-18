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
!>  This file contains the implementations of the tests of module [pm_mathExp](@ref pm_mathExp).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK_ENABLED
        use pm_kind, only: RKC => RKH
        integer(IKC), parameter :: LOWER = 2_IKC, UPPER = 100_IKC
        integer(IKC), parameter :: ZERO = 0_IKC
        integer(IKC)            :: absx(100)
        integer(IKC)            :: exponent
#elif   RK_ENABLED
        integer(IK)             :: exponent
        real(RKC)   , parameter :: LOWER = 1.01_RKC, UPPER = 100._RKC
        real(RKC)   , parameter :: ZERO = 0._RKC
        real(RKC)               :: absx(100)
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        call runTestsWith()
        call runTestsWith(base = getUnifRand(LOWER, UPPER))
        call runTestsWith(base = getUnifRand(LOWER, UPPER, s1 = 4_IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        impure elemental subroutine runTestsWith(base)
            integer(IK) :: i
#if         IK_ENABLED
            integer(IKC), intent(in), optional :: base
            absx = getUnifRand(1_IKC, nint(sqrt(real(huge(0_IKC), RKC)), kind = IKC), size(absx, 1, IK))
#elif       RK_ENABLED
            real(RKC), intent(in), optional :: base
            absx = getUnifRand(2 * epsilon(0._RKC), sqrt(huge(0._RKC)), size(absx, 1, IK))
#endif
            do i = 1, size(absx)
                exponent = getExpNext(absx(i), base)
                call report(absx(i), base)
            end do
            !exponent = getExpNext(ZERO, base)
            !call report(ZERO, base)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        impure elemental subroutine report(absx, base)
#if         IK_ENABLED
            integer(IKC), intent(in)            :: absx
            integer(IKC), intent(in), optional  :: base
            integer(IKC)                        :: base_def
            base_def = getOption(2_IKC, base)
#elif       RK_ENABLED
            real(RKC)   , intent(in)            :: absx
            real(RKC)   , intent(in), optional  :: base
            real(RKC)                           :: base_def
            base_def = getOption(2._RKC, base)
#endif
            assertion = assertion .and. logical(absx <= base_def**exponent, LK)
            assertion = assertion .and. logical(absx >= base_def**(max(0, int(exponent - 1))) .or. absx == ZERO, LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("exponent")
                call test%disp%show( exponent )
                call test%disp%show("present(base)")
                call test%disp%show( present(base) )
                call test%disp%show("base_def")
                call test%disp%show( base_def )
#if             getExpNext_ENABLED
                call test%disp%show("[real(RKC) :: base_def**(exponent - 1), absx, base_def**exponent]")
                call test%disp%show( [real(RKC) :: base_def**(exponent - 1), absx, base_def**exponent] )
#elif           getExpPrev_ENABLED
#else
#error          "Unrecognized interface."
#endif
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The test condition `absx <= base_def**exponent` must hold.", int(__LINE__, IK))
            call test%assert(assertion, SK_"The test condition `absx >= base_def**(exponent - 1)` must hold.", int(__LINE__, IK))
        end subroutine
