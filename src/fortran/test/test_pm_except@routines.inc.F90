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
!>  This include file contains the implementations of the tests of procedures with generic interfaces of [pm_except](@ref pm_except).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi, Sunday 4:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getInfPos_ENABLED || setInfPos_ENABLED || getInfNeg_ENABLED || setInfNeg_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getInfPos_ENABLED
#define test_getInf_ENABLED 1
#define GEN_INF getInfPos
#define IS_INF isInfPos
#elif   setInfPos_ENABLED
#define test_setInf_ENABLED 1
#define GET_INF setInfPos
#define IS_INF isInfPos
#elif   getInfNeg_ENABLED
#define test_getInf_ENABLED 1
#define GEN_INF getInfNeg
#define IS_INF isInfNeg
#elif   setInfNeg_ENABLED
#define test_setInf_ENABLED 1
#define GET_INF setInfNeg
#define IS_INF isInfNeg
#else
#error  "Unrecognized interface."
#endif

#if     CK_ENABLED
        complex(CKC), allocatable   :: Inf(:)
        real(CKC)   , allocatable   :: Dummy(:) ! \bug bypass Intel 2021.4 bug.
#elif   RK_ENABLED
        real(RKC)   , allocatable   :: Inf(:)
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK
        allocate(Inf(3))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_getInf_ENABLED
        Inf(1) = GEN_INF(Inf(1))
#elif   test_setInf_ENABLED
        call GET_INF(Inf(1))
#endif

        assertion = assertion .and. IS_INF(Inf(1))
        call report()
        call test%assert(assertion, SK_"GEN_INF()/GET_INF() must return must return a scalar `Inf`.", int(__LINE__, IK))

        assertion = assertion .and. isInf(Inf(1))
        call report()
        call test%assert(assertion, SK_"GEN_INF()/GET_INF() must return must return a scalar `Inf`.", int(__LINE__, IK))

#if     test_getInf_ENABLED
        Inf = GEN_INF(Inf)
#elif   test_setInf_ENABLED
        call GET_INF(Inf)
#endif

        assertion = assertion .and. all(IS_INF(Inf))
        call report()
        call test%assert(assertion, SK_"GEN_INF()/GET_INF() must return must return a vector `Inf`.", int(__LINE__, IK))

        assertion = assertion .and. all(isInf(Inf))
        call report()
        call test%assert(assertion, SK_"GEN_INF()/GET_INF() must return must return a vector `Inf`.", int(__LINE__, IK))

#if     test_getInf_ENABLED
        Inf(1) = GEN_INF(Inf(1))
        Inf(3) = GEN_INF(Inf(3))
#elif   test_setInf_ENABLED
        call GET_INF(Inf(1))
        call GET_INF(Inf(3))
#endif
        call setUnifRand(Inf(2))

        assertion = assertion .and. IS_INF(Inf(1)) .and. IS_INF(Inf(3)) .and. .not. IS_INF(Inf(2))
        call report()
        call test%assert(assertion, SK_"IS_INF() must properly recognize two `Inf` values in a vector of 3 values.", int(__LINE__, IK))

        assertion = assertion .and. isInf(Inf(1)) .and. isInf(Inf(3)) .and. .not. isInf(Inf(2))
        call report()
        call test%assert(assertion, SK_"isInf() must properly recognize two `Inf` values in a vector of 3 values.", int(__LINE__, IK))

#if     CK_ENABLED
        allocate(Dummy(size(Inf)))
#if     test_getInf_ENABLED
        Dummy(1) = GEN_INF(Dummy(1))
#elif   test_setInf_ENABLED
        call GET_INF(Dummy(1))
#endif
        Inf(1)%re = Dummy(1)
        call setUnifRand(Dummy(1))
        Inf(1)%im = Dummy(1)

        assertion = assertion .and. IS_INF(Inf(1)%re) .and. .not. IS_INF(Inf(1)%im)
        call report()
        call test%assert(assertion, SK_"IS_INF() must properly recognize a scalar `Inf` real component and a scalar `Inf` imaginary component.", int(__LINE__, IK))

        assertion = assertion .and. isInf(Inf(1)%re) .and. .not. isInf(Inf(1)%im)
        call report()
        call test%assert(assertion, SK_"isInf() must properly recognize a scalar `Inf` real component and a scalar `Inf` imaginary component.", int(__LINE__, IK))
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "Inf", Inf
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getNAN_ENABLED || setNAN_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#if     CK_ENABLED
        complex(CKC), allocatable   :: NAN(:)
        real(CKC)   , allocatable   :: Dummy(:) ! \bug bypass Intel 2021.4 bug.
#elif   RK_ENABLED
        real(RKC)   , allocatable   :: NAN(:)
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        allocate(NAN(3))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getNAN_ENABLED
        NAN(1) = getNAN(NAN(1))
#elif   setNAN_ENABLED
        call setNAN(NAN(1))
#endif
        assertion = assertion .and. isNAN(NAN(1))
        call report()
        call test%assert(assertion, SK_"getNAN() must return must return a scalar `NAN`.", int(__LINE__, IK))

#if     getNAN_ENABLED
        NAN = getNAN(NAN)
#elif   setNAN_ENABLED
        call setNAN(NAN)
#endif
        assertion = assertion .and. all(isNAN(NAN))
        call report()
        call test%assert(assertion, SK_"getNAN() must return must return a vector `NAN`.", int(__LINE__, IK))

#if     getNAN_ENABLED
        NAN(1) = getNAN(NAN(1))
        NAN(3) = getNAN(NAN(3))
#elif   setNAN_ENABLED
        call setNAN(NAN(1))
        call setNAN(NAN(3))
#endif
        call setUnifRand(NAN(2))
        assertion = assertion .and. isNAN(NAN(1)) .and. isNAN(NAN(3)) .and. .not. isNAN(NAN(2))
        call report()
        call test%assert(assertion, SK_"isNAN() must properly recognize two `NAN` values in a vector of 3 values.", int(__LINE__, IK))

#if     CK_ENABLED
        allocate(Dummy(size(NAN)))
#if     getNAN_ENABLED
        Dummy(1) = getNAN(Dummy(1))
#elif   setNAN_ENABLED
        call setNAN(Dummy(1))
#endif
        NAN(1)%re = Dummy(1)
        call setUnifRand(Dummy(1))
        NAN(1)%im = Dummy(1)
        assertion = assertion .and. isNAN(NAN(1)%re) .and. .not. isNAN(NAN(1)%im)
        call report()
        call test%assert(assertion, SK_"isNAN() must properly recognize a scalar `NAN` real component and a scalar `NAN` imaginary component.", int(__LINE__, IK))
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "NAN", NAN
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef test_getInf_ENABLED
#undef test_setInf_ENABLED
#undef GEN_INF
#undef GET_INF
#undef IS_INF