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
!>  This include file contains the implementations of procedures of [test_pm_complexDiv](@ref test_pm_complexDiv).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_val2str, only: getStr
        use pm_kind, only: SK, IK, LK
        use pm_distUnif, only: setUnifRand
        use pm_complexCompareAll, only: operator(<=)
        use pm_complexAbs, only: abs

#if     getDiv_ENABLED

        integer(IK)     , parameter :: NP = 5_IK
        complex(CKC)    , parameter :: SQRT_HUGE = sqrt(huge(0._CKC)), EPS = 100 * cmplx(epsilon(0._CKC), epsilon(0._CKC), CKC)
        complex(CKC)    , parameter :: LOWER = -SQRT_HUGE, UPPER = +SQRT_HUGE, ZERO = (0._CKC, 0._CKC)
        complex(CKC)    :: dividend(NP,NP), divisor(NP,NP), Quotient(NP,NP), Quotient_ref(NP,NP)
        integer(IK)     :: i, j, itest

        assertion = .true._LK

        do itest = 1, 10
            call runTests()
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTests()

            call setUnifRand(dividend, LOWER, UPPER)
            do
                call setUnifRand(divisor, LOWER, UPPER)
                if (all(divisor /= ZERO)) exit
            end do
            Quotient_ref = getDiv_ref(dividend, divisor)

            do i = 1, NP
                do j = 1, NP
                    Quotient(i,j) = getDiv(dividend(i,j), divisor(i,j))
                end do
            end do
            call report(int(__LINE__, IK))

            do j = 1, NP
                Quotient(:,j) = getDiv(dividend(:,j), divisor(:,j))
            end do
            call report(int(__LINE__, IK))

            Quotient(:,:) = getDiv(dividend(:,:), divisor(:,:))
            call report(int(__LINE__, IK))

            ! Finite-range tests.

            call setUnifRand(dividend)
            do
                call setUnifRand(divisor)
                if (all(divisor /= ZERO)) exit
            end do
            Quotient_ref = dividend / divisor
            Quotient(:,:) = getDiv(dividend(:,:), divisor(:,:))
            call report(int(__LINE__, IK))

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pure elemental function getDiv_ref(dividend, divisor) result(quotient)
            complex(CKC), intent(in)    :: dividend, divisor
            complex(CKC)                :: quotient
            real(CKC)                   :: r, d
            if (abs(divisor%re) < abs(divisor%im)) then
                r = divisor%re / divisor%im
                d = divisor%im + r * divisor%re
                quotient%re = (dividend%re * r + dividend%im) / d
                quotient%im = (dividend%im * r - dividend%re) / d
            else
                r = divisor%im / divisor%re
                d = divisor%re + r * divisor%im
                quotient%re = (dividend%re + dividend%im * r) / d
                quotient%im = (dividend%im - dividend%re * r) / d
            end if
        end function getDiv_ref

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line)
            integer(IK), intent(in) :: line
            integer(IK)             :: i, j
            complex(CKC)            :: diff
            do i = 1, NP
                do j = 1, NP
                    diff = abs(Quotient(i,j) - Quotient_ref(i,j))
                    assertion = assertion .and. diff <= EPS
                    if (test%traceable .and. .not. assertion) then
                        ! LCOV_EXCL_START
                        write(test%disp%unit,"(*(g0,:,', '))")
                        write(test%disp%unit,"(*(g0,:,', '))") "i, j", i, j
                        write(test%disp%unit,"(*(g0,:,', '))") "dividend(i,j)", dividend(i,j)
                        write(test%disp%unit,"(*(g0,:,', '))") "divisor(i,j)", divisor(i,j)
                        write(test%disp%unit,"(*(g0,:,', '))") "Quotient(i,j)", Quotient(i,j)
                        write(test%disp%unit,"(*(g0,:,', '))") "Quotient_ref(i,j)", Quotient_ref(i,j)
                        write(test%disp%unit,"(*(g0,:,', '))") "diff", diff
                        write(test%disp%unit,"(*(g0,:,', '))")
                        ! LCOV_EXCL_STOP
                    end if
                    call test%assert(assertion, SK_"The ratio of the complex values must be correctly computed without overflow.", line)
                end do
            end do
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else

#error  "Unrecognized interface."

#endif

#undef COMPARE
