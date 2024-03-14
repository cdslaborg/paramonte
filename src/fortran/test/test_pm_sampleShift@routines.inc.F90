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
!>  This file contains procedure implementations of tests of [pm_sampleMean](@ref pm_sampleMean).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the conjugation rule.
#if     CK_ENABLED
#define GET_CONJG(X)conjg(X)
#define TYPE_OF_SAMPLE complex(TKC)
        complex(TKC), parameter :: ZERO = 0._TKC, ONE = (1._TKC, 1._TKC), tol = cmplx(epsilon(1._TKC), epsilon(1._TKC), TKC) * 10
#elif   RK_ENABLED
#define GET_CONJG(X)X
#define TYPE_OF_SAMPLE real(TKC)
        real(TKC), parameter :: ZERO = 0._TKC, ONE = 1._TKC, tol = epsilon(1._TKC) * 10
#else
#error  "Unrecognized interface."
#endif
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), sampleShifted(:,:), diff(:,:), amount(:)
        integer(IK) :: itry, nsam, ndim, dim
        logical(LK) :: isPresentDim ! is dim present or not.
        logical(LK) :: isTransHerm
        logical(LK) :: isd1
        assertion = .true.
        do itry = 1, 100
            nsam = getUnifRand(1_IK, 5_IK)
            ndim = getUnifRand(1_IK, 5_IK)
            dim = merge(1, getChoice([1, 2]), isd1)
            isPresentDim = getUnifRand()
            isTransHerm = getUnifRand()
            if (dim == 2) then
                sample = getUnifRand(-ONE, ONE, ndim, nsam)
            else
                sample = getUnifRand(-ONE, ONE, nsam, ndim)
            end if
            amount = getUnifRand(0._TKC, 1._TKC, ndim)
            isd1 = ndim == 1 .and. dim == 1 .and. getUnifRand()
            sampleShifted = sample
#if         getShifted_ENABLED
            if (isd1) then
                if (isPresentDim) then
                    sampleShifted(:,1) = getShifted(getShifted(sampleShifted(:,1), dim, amount(1)), dim, -amount(1))
                else
                    sampleShifted(:,1) = getShifted(getShifted(sampleShifted(:,1), amount(1)), -amount(1))
                end if
            else
                if (isPresentDim) then
                    if (isTransHerm) then
                        sample = GET_CONJG(transpose(getShifted(sample, dim, -amount)))
                        sampleShifted = getShifted(sampleShifted, dim, -amount, transHerm)
                    else
                        sampleShifted = getShifted(getShifted(sampleShifted, dim, amount), dim, -amount)
                    end if
                else
                    if (isTransHerm) then
                        sample = GET_CONJG(transpose(getShifted(sample, -amount(1))))
                        sampleShifted = getShifted(sampleShifted, -amount(1), transHerm)
                    else
                        sampleShifted = getShifted(getShifted(sampleShifted, amount(1)), -amount(1))
                    end if
                end if
            end if
#elif       setShifted_ENABLED
            if (isd1) then
                call setShifted(sampleShifted(:,1), +amount(1))
                call setShifted(sampleShifted(:,1), -amount(1))
            else
                call setShifted(sampleShifted, dim, +amount)
                call setShifted(sampleShifted, dim, -amount)
            end if
#endif
            diff = abs(sample - sampleShifted)
            assertion = assertion .and. all(diff < tol)
            call report(__LINE__)
        end do

    contains

        subroutine report(line)
            integer, intent(in) :: line
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("isd1")
                call test%disp%show( isd1 )
                call test%disp%show("isTransHerm")
                call test%disp%show( isTransHerm )
                call test%disp%show("isPresentDim")
                call test%disp%show( isPresentDim )
                call test%disp%show("[ndim, nsam, dim]")
                call test%disp%show( [ndim, nsam, dim] )
                call test%disp%show("sample")
                call test%disp%show( sample )
                call test%disp%show("sampleShifted")
                call test%disp%show( sampleShifted )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The sample must be shifted correctly.", int(line, IK))
        end subroutine
#undef  TYPE_OF_SAMPLE
#undef  GET_CONJG