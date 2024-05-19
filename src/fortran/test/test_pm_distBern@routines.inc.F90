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
!>  This include file contains procedure implementations of the tests of [pm_distBern](@ref pm_distBern).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%
#if     isHead_ENABLED
        !%%%%%%%%%%%%%

        integer(IK)     :: i
        integer(IK)     , parameter :: NSIM = 20000_IK
        logical(LK)     :: rand(NSIM)

        assertion = .true._LK

        do i = 1_IK, NSIM
            rand(i) = isHead()
        end do
        assertion = assertion .and. logical(abs(NSIM / 2_IK - count(rand, kind = IK)) <= 5 * sqrt(real(NSIM / 2_IK)), LK)
        call test%assert(assertion, SK_"The procedure `isHead()` must be unbiased. There is a 5\sgima chance this test fails. Rerunning the test may resolve the failure.", int(__LINE__, IK))

        rand = isHead(size = NSIM)
        assertion = assertion .and. logical(abs(NSIM / 2_IK - count(rand, kind = IK)) <= 5 * sqrt(real(NSIM / 2_IK)), LK)
        call test%assert(assertion, SK_"The procedure `isHead(size = NSIM)` must be unbiased. There is a 5\sgima chance this test fails. Rerunning the test may resolve the failure.", int(__LINE__, IK))

        do i = 1_IK, NSIM
            rand(i) = isHead(p = 1._RKG)
        end do
        assertion = assertion .and. logical(all(rand), LK)
        call test%assert(assertion, SK_"The procedure `isHead(p = 1._RKG)` must always yield `.true.`.", int(__LINE__, IK))

        do i = 1_IK, NSIM
            rand(i) = isHead(p = 0._RKG)
        end do
        assertion = assertion .and. logical(.not. any(rand), LK)
        call test%assert(assertion, SK_"The procedure `isHead(p = 0._RKG)` must always yield `.false.`.", int(__LINE__, IK))

        rand = isHead(p = 1._RKG, size = NSIM)
        assertion = assertion .and. logical(all(rand), LK)
        call test%assert(assertion, SK_"The procedure `isHead(p = 1._RKG, size = NSIM)` must always yield `.true.`.", int(__LINE__, IK))

        rand = isHead(p = 0._RKG, size = NSIM)
        assertion = assertion .and. logical(.not. any(rand), LK)
        call test%assert(assertion, SK_"The procedure `isHead(p = 0._RKG, size = NSIM)` must always yield `.false.`.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%
#elif   getBernRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        integer(IK)     :: i
        integer(IK)     , parameter :: NSIM = 20000_IK
        integer(IK)     :: rand(NSIM)

        assertion = .true._LK

        do i = 1_IK, NSIM
            rand(i) = getBernRand(p = .5_RKG)
        end do
        assertion = assertion .and. logical(abs(NSIM / 2_IK - count(rand == 1_IK, kind = IK)) <= 5 * sqrt(real(NSIM / 2_IK)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p)` must be unbiased. There is a 5\sgima chance this test fails. Rerunning the test may resolve the failure.", int(__LINE__, IK))

        rand = getBernRand(p = .5_RKG, size = NSIM)
        assertion = assertion .and. logical(abs(NSIM / 2_IK - count(rand == 1_IK, kind = IK)) <= 5 * sqrt(real(NSIM / 2_IK)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p, size = NSIM)` must be unbiased. There is a 5\sgima chance this test fails. Rerunning the test may resolve the failure.", int(__LINE__, IK))

        do i = 1_IK, NSIM
            rand(i) = getBernRand(p = 1._RKG)
        end do
        assertion = assertion .and. logical(all(rand == 1_IK), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 1._RKG)` must always yield `1`.", int(__LINE__, IK))

        do i = 1_IK, NSIM
            rand(i) = getBernRand(p = 0._RKG)
        end do
        assertion = assertion .and. logical(.not. any(rand == 1_IK), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 0._RKG)` must always yield `0`.", int(__LINE__, IK))

        rand = getBernRand(p = 1._RKG, size = NSIM)
        assertion = assertion .and. logical(all(rand == 1_IK), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 1._RKG, size = NSIM)` must always yield `1`.", int(__LINE__, IK))

        rand = getBernRand(p = 0._RKG, size = NSIM)
        assertion = assertion .and. logical(.not. any(rand == 1_IK), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 0._RKG, size = NSIM)` must always yield `0`.", int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%
#elif   setBernRand_ENABLED
        !%%%%%%%%%%%%%%%%%%

        use pm_distUnif, only: getUnifRand
        integer(IK)     :: i
        integer(IK)     , parameter :: NSIM = 20000_IK
#if     IK_ENABLED
#define IS_TRUE(x) x == 1_IKG
        integer(IKG)    :: rand(NSIM)
#elif   LK_ENABLED
#define IS_TRUE(x) x
        logical(LKG)    :: rand(NSIM)
#elif   RK_ENABLED
#define IS_TRUE(x) x == 1._RKG
        real(RKG)       :: rand(NSIM)
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        do i = 1_IK, NSIM
            call setBernRand(rand(i), getUnifRand(0._RKG, 1._RKG), p = .5_RKG)
        end do
        assertion = assertion .and. logical(abs(NSIM / 2_IK - count(IS_TRUE(rand), kind = IK)) <= 5 * sqrt(real(NSIM / 2_IK)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p)` must be unbiased. There is a 5\sgima chance this test fails. Rerunning the test may resolve the failure.", int(__LINE__, IK))

        call setBernRand(rand, getUnifRand(0._RKG, 1._RKG, size(rand, 1, IK)), p = .5_RKG)
        assertion = assertion .and. logical(abs(NSIM / 2_IK - count(IS_TRUE(rand), kind = IK)) <= 5 * sqrt(real(NSIM / 2_IK)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p, size = NSIM)` must be unbiased. There is a 5\sgima chance this test fails. Rerunning the test may resolve the failure.", int(__LINE__, IK))

        do i = 1_IK, NSIM
            call setBernRand(rand(i), getUnifRand(0._RKG, 1._RKG), p = 1._RKG)
        end do
        assertion = assertion .and. logical(all(IS_TRUE(rand)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 1._RKG)` must always yield `1`.", int(__LINE__, IK))

        do i = 1_IK, NSIM
            call setBernRand(rand(i), getUnifRand(0._RKG, 1._RKG), p = 0._RKG)
        end do
        assertion = assertion .and. logical(.not. any(IS_TRUE(rand)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 0._RKG)` must always yield `0`.", int(__LINE__, IK))

        call setBernRand(rand, getUnifRand(0._RKG, 1._RKG, size(rand, 1, IK)), p = 1._RKG)
        assertion = assertion .and. logical(all(IS_TRUE(rand)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 1._RKG, size = NSIM)` must always yield `1`.", int(__LINE__, IK))

        call setBernRand(rand, getUnifRand(0._RKG, 1._RKG, size(rand, 1, IK)), p = 0._RKG)
        assertion = assertion .and. logical(.not. any(IS_TRUE(rand)), LK)
        call test%assert(assertion, SK_"The procedure `getBernRand(p = 0._RKG, size = NSIM)` must always yield `0`.", int(__LINE__, IK))

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  IS_TRUE