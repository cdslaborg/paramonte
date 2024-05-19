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
        use pm_kind, only: RKG => RKH
        integer(IKG), parameter :: LOWER = 2_IKG, UPPER = 100_IKG
        integer(IKG), parameter :: ZERO = 0_IKG
        integer(IKG)            :: absx(100)
        integer(IKG)            :: exponent
#elif   RK_ENABLED
        integer(IK)             :: exponent
        real(RKG)   , parameter :: LOWER = 1.01_RKG, UPPER = 100._RKG
        real(RKG)   , parameter :: ZERO = 0._RKG
        real(RKG)               :: absx(100)
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
            integer(IKG), intent(in), optional :: base
            absx = getUnifRand(1_IKG, nint(sqrt(real(huge(0_IKG), RKG)), kind = IKG), size(absx, 1, IK))
#elif       RK_ENABLED
            real(RKG), intent(in), optional :: base
            absx = getUnifRand(2 * epsilon(0._RKG), sqrt(huge(0._RKG)), size(absx, 1, IK))
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
            integer(IKG), intent(in)            :: absx
            integer(IKG), intent(in), optional  :: base
            integer(IKG)                        :: base_def
            base_def = getOption(2_IKG, base)
#elif       RK_ENABLED
            real(RKG)   , intent(in)            :: absx
            real(RKG)   , intent(in), optional  :: base
            real(RKG)                           :: base_def
            base_def = getOption(2._RKG, base)
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
                call test%disp%show("[real(RKG) :: base_def**(exponent - 1), absx, base_def**exponent]")
                call test%disp%show( [real(RKG) :: base_def**(exponent - 1), absx, base_def**exponent] )
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
