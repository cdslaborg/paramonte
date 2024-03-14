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
!>  This file contains procedure implementations of tests of [pm_sampleCor](@ref pm_sampleCor).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the comparison precision and tolerance.
#if     CK_ENABLED && (getCor_ENABLED || setCor_ENABLED)
        real(TKC), parameter :: rtol = epsilon(1._TKC) * 100
        complex(TKC), parameter :: ctol = (rtol, rtol)
#define GET_CONJG(X)conjg(X)
#elif   RK_ENABLED && (getCor_ENABLED || setCor_ENABLED)
        real(TKC), parameter :: rtol = epsilon(1._TKC) * 100
        real(TKC), parameter :: ctol = rtol
#define GET_CONJG(X)X
#elif   getCor_ENABLED || setCor_ENABLED
#error  "Unrecognized interface."
#else
        real(RK), parameter :: rtol = epsilon(1._RK) * 100
#endif
        ! Define the sample type.
#if     SK_ENABLED && D0_ENABLED
#define GET_SHAPE len
#define TYPE_OF_SAMPLE character(:,TKC)
        character(1,TKC), parameter :: lb = TKC_"a", ub = TKC_"z"
#else
#define GET_SHAPE shape
#if     SK_ENABLED
#define TYPE_OF_SAMPLE character(2,TKC)
        TYPE_OF_SAMPLE, parameter :: lb = TKC_"aa", ub = TKC_"zz"
#elif   IK_ENABLED
#define TYPE_OF_SAMPLE integer(TKC)
        TYPE_OF_SAMPLE, parameter :: lb = -9_TKC, ub = +9_TKC
#elif   CK_ENABLED
#define TYPE_OF_SAMPLE complex(TKC)
        TYPE_OF_SAMPLE, parameter :: lb = (2._TKC, -3._TKC), ub = (3._TKC, -2._TKC), ZERO = (0._TKC, 0._TKC), ONE = (1._TKC, 0._TKC), ONES = (1._TKC, 1._TKC), TWOS = ONES + ONES
#elif   RK_ENABLED
#define TYPE_OF_SAMPLE real(TKC)
        TYPE_OF_SAMPLE, parameter :: lb = 2._TKC, ub = 3._TKC, ZERO = 0._TKC, ONE = 1._TKC, ONES = 1._TKC, TWOS = ONES + ONES
#elif   PSSK_ENABLED
#define TYPE_OF_SAMPLE type(css_pdt(TKC))
        character(1,TKC), parameter :: lb = TKC_"a", ub = TKC_"z"
#elif   BSSK_ENABLED
#define TYPE_OF_SAMPLE type(css_type)
        character(1,TKC), parameter :: lb = TKC_"a", ub = TKC_"z"
#elif   !setCordance_ENABLED
#error  "Unrecognized interface."
#endif
#endif
        !%%%%%%%%%%%%%
#if     getCor_ENABLED
        !%%%%%%%%%%%%%

        real(TKC) :: rweisum
        integer(IK) :: iweisum
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam, ndim, dim
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:)
        TYPE_OF_SAMPLE, allocatable :: cor(:,:), cov(:,:), covupp(:,:), covlow(:,:), cor_ref(:,:), cdiff(:,:)
        real(TKC), allocatable :: std(:), stdinv(:)
        assertion = .true._LK

        do itry = 1, 50

            ! Test CFC: The strategy is to construct a cot mat, convert to cov mat using `setCov`, and recover the original cor mat using `setCor`.

            block

                ndim = getUnifRand(1_IK, 5_IK)
                cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                call setMatCopy(cor_ref, rdpack, cor_ref, rdpack, upp, transHerm)
                call setMatInit(cor_ref, dia, ONE)
                std = getUnifRand(1._TKC, 2._TKC, ndim)
                stdinv = 1._TKC / std
                cov = getCov(cor_ref, uppDia, std)
                covupp = cov; call setMatInit(covupp, low, ZERO)
                covlow = cov; call setMatInit(covlow, upp, ZERO)

                cor = getCor(covupp, uppDia)
                call reportCFC(__LINE__)

                cor = getCor(covlow, lowDia)
                call reportCFC(__LINE__)

                call setMatInit(covupp, dia, ZERO)
                call setMatInit(covlow, dia, ZERO)

                cor = getCor(covupp, upp, stdinv)
                call reportCFC(__LINE__, stdinv)

                cor = getCor(covlow, low, stdinv)
                call reportCFC(__LINE__, stdinv)

            end block

            ! Test Prs: The strategy is to construct a sample, compute its correlation and compare it against the reference computed here.

            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample interface.

            block

                type(uppDia_type), parameter :: subsetr = uppDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
                    sample = getUnifRand(ONES, TWOS, ndim, nsam)
                else
                    sample = getUnifRand(ONES, TWOS, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                cor_ref = getFilled(ZERO, ndim, ndim)
                cor = getFilled(ZERO, ndim, ndim)

                ! integer weighted

                mean = getMean(sample, dim, iweight)

                cor = getCor(sample, dim, iweight)
                call setCor(cor_ref, subsetr, mean, sample, dim, iweight, iweisum)
                call setMatCopy(cor_ref, rdpack, cor_ref, rdpack, subsetr, transHerm)
                call report(__LINE__, SK_"integer-weighted")

                ! real weighted

                mean = getMean(sample, dim, rweight)

                cor = getCor(sample, dim, rweight)
                call setCor(cor_ref, subsetr, mean, sample, dim, rweight, rweisum)
                call setMatCopy(cor_ref, rdpack, cor_ref, rdpack, subsetr, transHerm)
                call report(__LINE__, SK_"real-weighted")

                ! unweighted

                mean = getMean(sample, dim)

                cor = getCor(sample, dim)
                call setCor(cor_ref, subsetr, mean, sample, dim)
                call setMatCopy(cor_ref, rdpack, cor_ref, rdpack, subsetr, transHerm)
                call report(__LINE__, SK_"unweighted")

            end block

            ! test sample XY interface.

            block

                dim = 1
                ndim = 2
                sample = getUnifRand(ONES, TWOS, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                cor_ref = getFilled(ZERO, ndim, ndim)
                cor = getFilled(ZERO, ndim, ndim)

                ! integer weighted

                mean = getMean(sample, dim, iweight)

                cor(1,2) = getCor(sample(:,1), sample(:,2), iweight)
                call setCor(cor_ref(1,2), mean, sample(:,1), sample(:,2), iweight, iweisum)
                call report(__LINE__, SK_"integer-weighted")

                ! real weighted

                mean = getMean(sample, dim, rweight)

                cor(1,2) = getCor(sample(:,1), sample(:,2), rweight)
                call setCor(cor_ref(1,2), mean, sample(:,1), sample(:,2), rweight, rweisum)
                call report(__LINE__, SK_"real-weighted")

                ! unweighted

                mean = getMean(sample, dim)

                cor(1,2) = getCor(sample(:,1), sample(:,2))
                call setCor(cor_ref(1,2), mean, sample(:,1), sample(:,2))
                call report(__LINE__, SK_"unweighted")

            end block

        end do

    contains

        subroutine reportCFC(line, stdinv)
            real(TKC), intent(in), optional :: stdinv(:)
            integer, intent(in) :: line
            cdiff = abs(cor - cor_ref)
            assertion = assertion .and. all(cdiff < ctol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("ndim")
                call test%disp%show( ndim )
                call test%disp%show("present(stdinv)")
                call test%disp%show( present(stdinv) )
                if (present(stdinv)) then
                    call test%disp%show("stdinv")
                    call test%disp%show( stdinv )
                end if
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cor")
                call test%disp%show( cor )
                call test%disp%show("cor_ref")
                call test%disp%show( cor_ref )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cor` must be correctly computed from the specified `cor` and `stdinv`.", int(line, IK))
        end subroutine

        subroutine report(line, this)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: this
            cdiff = abs(cor - cor_ref)
            assertion = assertion .and. all(cdiff < ctol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[ndim, nsam, dim, shape(sample, IK)]")
                call test%disp%show( [ndim, nsam, dim, shape(sample, IK)] )
                call test%disp%show("sample")
                call test%disp%show( sample )
                call test%disp%show("iweight")
                call test%disp%show( iweight )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("cor_ref")
                call test%disp%show( cor_ref )
                call test%disp%show("cor")
                call test%disp%show( cor )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cor` must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%
#elif   setCor_ENABLED
        !%%%%%%%%%%%%%

        real(TKC) :: rweisum
        integer(IK) :: iweisum
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam, ndim, dim
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:)
        TYPE_OF_SAMPLE, allocatable :: cor(:,:), cov(:,:), cor_ref(:,:), cdiff(:,:)
        real(TKC), allocatable :: std(:), stdinv(:)
        assertion = .true._LK

        do itry = 1, 50

            ! Test CFC: The strategy is to construct a cot mat, convert to cov mat using `setCov`, and recover the original cor mat using `setCor`.

            block

                block

                    type(low_type) :: subsetr
                    type(low_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubComp(subsetr), ZERO, ZERO)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block

                block

                    type(low_type) :: subsetr
                    type(upp_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubComp(subsetr), ZERO, ZERO)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block

                block

                    type(upp_type) :: subsetr
                    type(low_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubComp(subsetr), ZERO, ZERO)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block

                block

                    type(upp_type) :: subsetr
                    type(upp_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubComp(subsetr), ZERO, ZERO)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block

                block

                    type(lowDia_type) :: subsetr
                    type(low_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubSymm(subsetr), ZERO, ONE)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block

                block

                    type(lowDia_type) :: subsetr
                    type(upp_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubSymm(subsetr), ZERO, ONE)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block

                block

                    type(uppDia_type) :: subsetr
                    type(low_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubSymm(subsetr), ZERO, ONE)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block

                block

                    type(uppDia_type) :: subsetr
                    type(upp_type) :: subsetv
                    ndim = getUnifRand(1_IK, 5_IK)
                    cor_ref = getUnifRand(-ONES, ONES, ndim, ndim)
                    call setMatInit(cor_ref, getSubSymm(subsetr), ZERO, ONE)
                    std = getUnifRand(1._TKC, 2._TKC, ndim)
                    stdinv = 1._TKC / std
                    cov = getFilled(ZERO, ndim, ndim)
                    call setCov(cov, getSubUnion(subsetv, dia), cor_ref, getSubUnion(subsetr, dia), std)

                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, getSubUnion(subsetv, dia))
                    call reportCFC(__LINE__)

                    call setMatInit(cov, dia, ZERO)
                    cor = getFilled(ZERO, ndim, ndim)
                    call setCor(cor, subsetr, cov, subsetv, stdinv)
                    call reportCFC(__LINE__, stdinv)

                end block
            end block

            ! Test Prs: The strategy is to construct a sample, compute its correlation and compare it against the reference computed here.

            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample lowDia interface.

            block

                type(lowDia_type), parameter :: subsetr = lowDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
                    sample = getUnifRand(ONES, TWOS, ndim, nsam)
                else
                    sample = getUnifRand(ONES, TWOS, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                cor_ref = getCorRef(mean, sample, dim, real(iweight, TKC))

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, mean, sample, dim, iweight, iweisum)
                call report(__LINE__, subsetr, SK_"integer-weighted")

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, getShifted(sample, dim, -mean), dim, iweight, iweisum)
                call report(__LINE__, subsetr, SK_"integer-weighted shifted")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                cor_ref = getCorRef(mean, sample, dim, rweight)

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, mean, sample, dim, rweight, rweisum)
                call report(__LINE__, subsetr, SK_"real-weighted")

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, getShifted(sample, dim, -mean), dim, rweight, rweisum)
                call report(__LINE__, subsetr, SK_"real-weighted shifted")

                ! unweighted

                mean = getMean(sample, dim)
                cor_ref = getCorRef(mean, sample, dim)

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, mean, sample, dim)
                call report(__LINE__, subsetr, SK_"unweighted")

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, getShifted(sample, dim, -mean), dim)
                call report(__LINE__, subsetr, SK_"unweighted shifted")

            end block

            ! test sample uppDia interface.

            block

                type(uppDia_type), parameter :: subsetr = uppDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
                    sample = getUnifRand(ONES, TWOS, ndim, nsam)
                else
                    sample = getUnifRand(ONES, TWOS, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                cor_ref = getCorRef(mean, sample, dim, real(iweight, TKC))

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, mean, sample, dim, iweight, iweisum)
                call report(__LINE__, subsetr, SK_"integer-weighted")

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, getShifted(sample, dim, -mean), dim, iweight, iweisum)
                call report(__LINE__, subsetr, SK_"integer-weighted shifted")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                cor_ref = getCorRef(mean, sample, dim, rweight)

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, mean, sample, dim, rweight, rweisum)
                call report(__LINE__, subsetr, SK_"real-weighted")

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, getShifted(sample, dim, -mean), dim, rweight, rweisum)
                call report(__LINE__, subsetr, SK_"real-weighted shifted")

                ! unweighted

                mean = getMean(sample, dim)
                cor_ref = getCorRef(mean, sample, dim)

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, mean, sample, dim)
                call report(__LINE__, subsetr, SK_"unweighted")

                cor = getFilled(ZERO, ndim, ndim)
                call setCor(cor, subsetr, getShifted(sample, dim, -mean), dim)
                call report(__LINE__, subsetr, SK_"unweighted shifted")

            end block

            ! test sample XY interface.

            block

                type(lowDia_type), parameter :: subsetr = lowDia_type()
                dim = 1_IK
                ndim = 2_IK
                mean = getFilled(ZERO, ndim)
                cor_ref = getFilled(ZERO, ndim, ndim)
                sample = getUnifRand(ONES, TWOS, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                cor_ref = getCorRef(mean, sample, dim, real(iweight, TKC))

                cor = getMatInit([ndim, ndim], uppLowDia, ZERO, ZERO, ONE)
                call setCor(cor(1,2), mean, sample(:,1), sample(:,2), iweight, iweisum)
                call report(__LINE__, uppLowDia, SK_"integer-weighted")

                cor = getMatInit([ndim, ndim], uppLowDia, ZERO, ZERO, ONE)
                call setCor(cor(1,2), sample(:,1) - mean(1), sample(:,2) - mean(2), iweight, iweisum)
                call report(__LINE__, uppLowDia, SK_"integer-weighted")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                cor_ref = getCorRef(mean, sample, dim, rweight)

                cor = getMatInit([ndim, ndim], uppLowDia, ZERO, ZERO, ONE)
                call setCor(cor(1,2), mean, sample(:,1), sample(:,2), rweight, rweisum)
                call report(__LINE__, uppLowDia, SK_"real-weighted")

                cor = getMatInit([ndim, ndim], uppLowDia, ZERO, ZERO, ONE)
                call setCor(cor(1,2), sample(:,1) - mean(1), sample(:,2) - mean(2), rweight, rweisum)
                call report(__LINE__, uppLowDia, SK_"real-weighted")

                ! unweighted

                mean = getMean(sample, dim)
                cor_ref = getCorRef(mean, sample, dim)

                cor = getMatInit([ndim, ndim], uppLowDia, ZERO, ZERO, ONE)
                call setCor(cor(1,2), mean, sample(:,1), sample(:,2))
                call report(__LINE__, uppLowDia, SK_"unweighted")

                cor = getMatInit([ndim, ndim], uppLowDia, ZERO, ZERO, ONE)
                call setCor(cor(1,2), sample(:,1) - mean(1), sample(:,2) - mean(2))
                call report(__LINE__, uppLowDia, SK_"unweighted")

            end block

        end do

    contains

        subroutine reportCFC(line, stdinv)
            real(TKC), intent(in), optional :: stdinv(:)
            integer, intent(in) :: line
            cdiff = abs(cor - cor_ref)
            assertion = assertion .and. all(cdiff < ctol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("ndim")
                call test%disp%show( ndim )
                call test%disp%show("present(stdinv)")
                call test%disp%show( present(stdinv) )
                if (present(stdinv)) then
                    call test%disp%show("stdinv")
                    call test%disp%show( stdinv )
                end if
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cor")
                call test%disp%show( cor )
                call test%disp%show("cor_ref")
                call test%disp%show( cor_ref )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cor` must be correctly computed from the specified `cor` and `stdinv`.", int(line, IK))
        end subroutine

        PURE function getCorRef(mean, sample, dim, weight) result(cor)
            ! returns the full correlation matrix of the sample.
            integer(IK), intent(in) :: dim
            real(TKC), intent(in), optional :: weight(:)
            TYPE_OF_SAMPLE, intent(in), contiguous :: sample(:,:), mean(:)
            TYPE_OF_SAMPLE :: cor(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
            TYPE_OF_SAMPLE :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
            integer(IK) :: idim, jdim, ndim
            real(TKC) :: stdinv(size(sample, 3 - dim, IK))
            ndim = size(sample, 3 - dim, IK)
            sampleShifted = getShifted(sample, dim, -mean)
            if (present(weight)) then
                call setCov(cor, uppDia, sampleShifted, dim, weight, sum(weight))
            else
                call setCov(cor, uppDia, sampleShifted, dim)
            end if
            do jdim = 1, ndim
                stdinv(jdim) = 1._TKC / sqrt(real(cor(jdim, jdim), TKC))
                cor(jdim, jdim) = 1._TKC
                do idim = jdim - 1, 1, -1
                    cor(idim, jdim) = cor(idim, jdim) * stdinv(idim) * stdinv(jdim)
                    cor(jdim, idim) = GET_CONJG(cor(idim, jdim))
                end do
            end do
        end function

        subroutine report(line, subsetr, this)
            integer, intent(in) :: line
            class(*), intent(in) :: subsetr
            character(*, SK), intent(in) :: this
            if (same_type_as(subsetr, uppDia)) then
                call setMatCopy(cor, rdpack, cor, rdpack, uppDia, transHerm)
            elseif (same_type_as(subsetr, lowDia)) then
                call setMatCopy(cor, rdpack, cor, rdpack, lowDia, transHerm)
            elseif (same_type_as(subsetr, uppLowDia)) then
                cor(2,1) = GET_CONJG(cor(1,2))
            else
                error stop "Unrecognized subset."
            end if
            cdiff = abs(cor - cor_ref)
            assertion = assertion .and. all(cdiff < ctol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[ndim, nsam, dim, shape(sample, IK)]")
                call test%disp%show( [ndim, nsam, dim, shape(sample, IK)] )
                call test%disp%show("sample")
                call test%disp%show( sample )
                call test%disp%show("iweight")
                call test%disp%show( iweight )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("cor_ref")
                call test%disp%show( cor_ref )
                call test%disp%show("cor")
                call test%disp%show( cor )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cor` must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getRho_ENABLED && (XY_ENABLED || D0_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: itry, nsam
        real(RK), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        real(RK), allocatable :: frankx(:), franky(:)
#if     D0_ENABLED
        TYPE_OF_SAMPLE, allocatable :: x, y
#elif   XY_ENABLED
        TYPE_OF_SAMPLE, allocatable :: x(:), y(:)
#else
#error  "Unrecognized interface."
#endif
        real(RK) :: rho, rho_ref, rdiff
        assertion = .true._LK
        do itry = 1, 50

            nsam = getUnifRand(20_IK, 30_IK)
#if         D0_ENABLED && SK_ENABLED
            x = getUnifRand(repeat(lb, nsam), repeat(ub, nsam))
            y = getUnifRand(repeat(lb, nsam), repeat(ub, nsam))
#elif       XY_ENABLED && (BSSK_ENABLED || PSSK_ENABLED)
            if (allocated(x)) deallocate(x); allocate(x(nsam))
            if (allocated(y)) deallocate(y); allocate(y(nsam))
            block
                integer(IK) :: isam
                do isam = 1, nsam
                    x(isam)%val = getUnifRand(lb, ub)
                    y(isam)%val = getUnifRand(lb, ub)
                end do
            end block
#elif       XY_ENABLED
            x = getUnifRand(lb, ub, nsam)
            y = getUnifRand(lb, ub, nsam)
#else
#error      "Unrecognized interface."
#endif
            iweight = getUnifRand(1_IK, 9_IK, nsam)
            rweight = getUnifRand(1._RK, 9._RK, nsam)

            ! integer weighted

            rho_ref = getRhoCoef(real(iweight, RK))

            rho = getRho(x, y, iweight)
            call report(__LINE__, SK_"integer-weighted")

            ! real weighted

            rho_ref = getRhoCoef(rweight)

            rho = getRho(x, y, rweight)
            call report(__LINE__, SK_"real-weighted")

            ! unweighted

            rho_ref = getRhoCoef()

            rho = getRho(x, y)
            call report(__LINE__, SK_"unweighted")

        end do

    contains

        function getRhoCoef(weight) result(rho)
            ! returns the full correlation matrix of the sample.
            real(RK), intent(in), optional :: weight(:)
            real(RK) :: rho
            frankx = getRankFractional(x)
            franky = getRankFractional(y)
            if (present(weight)) then
                rho = getCor(frankx, franky, weight)
            else
                rho = getCor(frankx, franky)
            end if
        end function

        subroutine report(line, this)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: this
            rdiff = abs(rho - rho_ref)
            assertion = assertion .and. rdiff < rtol
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("nsam")
                call test%disp%show( nsam )
                call test%disp%show("x")
                call test%disp%show( x , deliml = SK_"""" )
                call test%disp%show("y" )
                call test%disp%show( y , deliml = SK_"""" )
                call test%disp%show("iweight")
                call test%disp%show( iweight )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("rho_ref")
                call test%disp%show( rho_ref )
                call test%disp%show("rho")
                call test%disp%show( rho )
                call test%disp%show("rdiff")
                call test%disp%show( rdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `rho` and fractional ranks must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getRho_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RK), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam
        real(RK), allocatable :: rho(:,:), rho_ref(:,:), rdiff(:,:)
        TYPE_OF_SAMPLE, allocatable :: sample(:,:)
        integer(IK) :: ndim, dim
        assertion = .true._LK

        do itry = 1, 50

            ! Test Rho: The strategy is to construct a sample, compute its correlation and compare it against the reference computed here.

            nsam = getUnifRand(20_IK, 30_IK)

            ! test sample lowDia interface.

            block

                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
#if                 PSSK_ENABLED
                    sample = css_pdt(getUnifRand(lb, ub, ndim, nsam))
#elif               BSSK_ENABLED
                    sample = css_type(getUnifRand(lb, ub, ndim, nsam))
#else
                    sample = getUnifRand(lb, ub, ndim, nsam)
#endif
                else
#if                 PSSK_ENABLED
                    sample = css_pdt(getUnifRand(lb, ub, nsam, ndim))
#elif               BSSK_ENABLED
                    sample = css_type(getUnifRand(lb, ub, nsam, ndim))
#else
                    sample = getUnifRand(lb, ub, nsam, ndim)
#endif
                end if
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                rweight = getUnifRand(1., 9., nsam)

                ! integer weighted

                rho_ref = getRhoRef(real(iweight, RK))
                rho = getRho(sample, dim, iweight)
                call report(__LINE__, SK_"integer-weighted")

                ! real weighted

                rho_ref = getRhoRef(rweight)
                rho = getRho(sample, dim, rweight)
                call report(__LINE__, SK_"real-weighted")

                ! unweighted

                rho_ref = getRhoRef()
                rho = getRho(sample, dim)
                call report(__LINE__, SK_"unweighted")

            end block

        end do

    contains

        function getRhoRef(weight) result(rho)
            ! returns the full correlation matrix of the sample.
            real(RK), intent(in), optional :: weight(:)
            real(RK) :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
            real(RK) :: frank(size(sample, 1, IK), size(sample, 2, IK))
            integer(IK) :: idim, ndim, nsam
            ndim = size(sample, 3 - dim, IK)
            nsam = size(sample, dim, IK)
            if (dim == 2_IK) then
                do idim = 1, ndim
                    frank(idim, :) = getRankFractional(sample(idim, :))
                end do
            else
                do idim = 1, ndim
                    frank(:, idim) = getRankFractional(sample(:, idim))
                end do
            end if
            if (present(weight)) then
                rho = getCor(frank, dim, weight)
            else
                rho = getCor(frank, dim)
            end if
        end function

        subroutine report(line, this)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: this
            rdiff = abs(rho - rho_ref)
            assertion = assertion .and. all(rdiff < rtol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[ndim, nsam, dim, GET_SHAPE(sample, IK)]")
                call test%disp%show( [ndim, nsam, dim, GET_SHAPE(sample, IK)] )
                call test%disp%show("sample")
                call test%disp%show( sample )
                call test%disp%show("iweight")
                call test%disp%show( iweight )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("rho_ref")
                call test%disp%show( rho_ref )
                call test%disp%show("rho")
                call test%disp%show( rho )
                call test%disp%show("rdiff")
                call test%disp%show( rdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `rho` and fractional ranks must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRho_ENABLED && (XY_ENABLED || D0_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: itry, nsam
        real(RK), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        real(RK), allocatable :: frankx(:), franky(:), frankx_ref(:), franky_ref(:)
#if     D0_ENABLED
        TYPE_OF_SAMPLE, allocatable :: x, y
#elif   XY_ENABLED
        TYPE_OF_SAMPLE, allocatable :: x(:), y(:)
#else
#error  "Unrecognized interface."
#endif
        real(RK) :: rho, rho_ref, rdiff
        assertion = .true._LK
        do itry = 1, 50

            nsam = getUnifRand(20_IK, 30_IK)
#if         D0_ENABLED && SK_ENABLED
            x = getUnifRand(repeat(lb, nsam), repeat(ub, nsam))
            y = getUnifRand(repeat(lb, nsam), repeat(ub, nsam))
#elif       XY_ENABLED && (BSSK_ENABLED || PSSK_ENABLED)
            if (allocated(x)) deallocate(x); allocate(x(nsam))
            if (allocated(y)) deallocate(y); allocate(y(nsam))
            block
                integer(IK) :: isam
                do isam = 1, nsam
                    x(isam)%val = getUnifRand(lb, ub)
                    y(isam)%val = getUnifRand(lb, ub)
                end do
            end block
#elif       XY_ENABLED
            x = getUnifRand(lb, ub, nsam)
            y = getUnifRand(lb, ub, nsam)
#else
#error      "Unrecognized interface."
#endif
            frankx = getFilled(0._RK, nsam)
            franky = getFilled(0._RK, nsam)
            iweight = getUnifRand(1_IK, 9_IK, nsam)
            rweight = getUnifRand(1._RK, 9._RK, nsam)

            ! integer weighted

            rho_ref = getRhoCoef(real(iweight, RK))

            call setRho(rho, frankx, franky, x, y, iweight)
            call report(__LINE__, SK_"integer-weighted")

            ! real weighted

            rho_ref = getRhoCoef(rweight)

            call setRho(rho, frankx, franky, x, y, rweight)
            call report(__LINE__, SK_"real-weighted")

            ! unweighted

            rho_ref = getRhoCoef()

            call setRho(rho, frankx, franky, x, y)
            call report(__LINE__, SK_"unweighted")

        end do

    contains

        function getRhoCoef(weight) result(rho)
            ! returns the full correlation matrix of the sample.
            real(RK), intent(in), optional :: weight(:)
            real(RK) :: rho
            frankx_ref = getRankFractional(x)
            franky_ref = getRankFractional(y)
            if (present(weight)) then
                rho = getCor(frankx_ref, franky_ref, weight)
            else
                rho = getCor(frankx_ref, franky_ref)
            end if
        end function

        subroutine report(line, this)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: this
            rdiff = abs(rho - rho_ref)
            assertion = assertion .and. rdiff < rtol
            assertion = assertion .and. all(abs(frankx - frankx_ref) < rtol)
            assertion = assertion .and. all(abs(franky - franky_ref) < rtol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("nsam")
                call test%disp%show( nsam )
                call test%disp%show("x")
                call test%disp%show( x , deliml = SK_"""" )
                call test%disp%show("y" )
                call test%disp%show( y , deliml = SK_"""" )
                call test%disp%show("iweight")
                call test%disp%show( iweight )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("frankx")
                call test%disp%show( frankx )
                call test%disp%show("frankx_ref")
                call test%disp%show( frankx_ref )
                call test%disp%show("franky")
                call test%disp%show( franky )
                call test%disp%show("franky_ref")
                call test%disp%show( franky_ref )
                call test%disp%show("rho_ref")
                call test%disp%show( rho_ref )
                call test%disp%show("rho")
                call test%disp%show( rho )
                call test%disp%show("rdiff")
                call test%disp%show( rdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            ! Three tests, so three calls.
            call test%assert(assertion, SK_"The `rho` and fractional ranks must be computed correctly for "//this//SK_" sample.", int(line, IK))
            call test%assert(assertion, SK_"The `rho` and fractional ranks must be computed correctly for "//this//SK_" sample.", int(line, IK))
            call test%assert(assertion, SK_"The `rho` and fractional ranks must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setRho_ENABLED && D2_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RK), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam
        real(RK), allocatable :: frank(:,:), frank_ref(:,:), rho(:,:), rho_ref(:,:), rdiff(:,:)
        TYPE_OF_SAMPLE, allocatable :: sample(:,:)
        integer(IK) :: ndim, dim
        assertion = .true._LK

        do itry = 1, 50

            ! Test Rho: The strategy is to construct a sample, compute its correlation and compare it against the reference computed here.

            nsam = getUnifRand(20_IK, 30_IK)

            ! test sample lowDia interface.

            block

                type(lowDia_type), parameter :: subsetr = lowDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
#if                 PSSK_ENABLED
                    sample = css_pdt(getUnifRand(lb, ub, ndim, nsam))
#elif               BSSK_ENABLED
                    sample = css_type(getUnifRand(lb, ub, ndim, nsam))
#else
                    sample = getUnifRand(lb, ub, ndim, nsam)
#endif
                else
#if                 PSSK_ENABLED
                    sample = css_pdt(getUnifRand(lb, ub, nsam, ndim))
#elif               BSSK_ENABLED
                    sample = css_type(getUnifRand(lb, ub, nsam, ndim))
#else
                    sample = getUnifRand(lb, ub, nsam, ndim)
#endif
                end if
                frank = getFilled(0._RK, size(sample, 1, IK), size(sample, 2, IK))
                frank_ref = frank
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                rweight = getUnifRand(1., 9., nsam)

                ! integer weighted

                rho_ref = getRhoRef(real(iweight, RK))

                rho = getFilled(0._RK, ndim, ndim)
                call setRho(rho, subsetr, frank, sample, dim, iweight)
                call report(__LINE__, subsetr, SK_"integer-weighted")

                ! real weighted

                rho_ref = getRhoRef(rweight)

                rho = getFilled(0._RK, ndim, ndim)
                call setRho(rho, subsetr, frank, sample, dim, rweight)
                call report(__LINE__, subsetr, SK_"real-weighted")

                ! unweighted

                rho_ref = getRhoRef()

                rho = getFilled(0._RK, ndim, ndim)
                call setRho(rho, subsetr, frank, sample, dim)
                call report(__LINE__, subsetr, SK_"unweighted")

            end block

            ! test sample uppDia interface.

            block

                type(uppDia_type), parameter :: subsetr = uppDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
#if                 PSSK_ENABLED
                    sample = css_pdt(getUnifRand(lb, ub, ndim, nsam))
#elif               BSSK_ENABLED
                    sample = css_type(getUnifRand(lb, ub, ndim, nsam))
#else
                    sample = getUnifRand(lb, ub, ndim, nsam)
#endif
                else
#if                 PSSK_ENABLED
                    sample = css_pdt(getUnifRand(lb, ub, nsam, ndim))
#elif               BSSK_ENABLED
                    sample = css_type(getUnifRand(lb, ub, nsam, ndim))
#else
                    sample = getUnifRand(lb, ub, nsam, ndim)
#endif
                end if
                frank = getFilled(0._RK, size(sample, 1, IK), size(sample, 2, IK))
                frank_ref = frank
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                rweight = getUnifRand(1., 9., nsam)

                ! integer weighted

                rho_ref = getRhoRef(real(iweight, RK))

                rho = getFilled(0._RK, ndim, ndim)
                call setRho(rho, subsetr, frank, sample, dim, iweight)
                call report(__LINE__, subsetr, SK_"integer-weighted")

                ! real weighted

                rho_ref = getRhoRef(rweight)

                rho = getFilled(0._RK, ndim, ndim)
                call setRho(rho, subsetr, frank, sample, dim, rweight)
                call report(__LINE__, subsetr, SK_"real-weighted")

                ! unweighted

                rho_ref = getRhoRef()

                rho = getFilled(0._RK, ndim, ndim)
                call setRho(rho, subsetr, frank, sample, dim)
                call report(__LINE__, subsetr, SK_"unweighted")

            end block

        end do

    contains

        function getRhoRef(weight) result(rho)
            ! returns the full correlation matrix of the sample.
            real(RK), intent(in), optional :: weight(:)
            real(RK) :: rho(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
            integer(IK) :: idim, ndim, nsam
            ndim = size(sample, 3 - dim, IK)
            nsam = size(sample, dim, IK)
            if (dim == 2_IK) then
                do idim = 1, ndim
                    frank_ref(idim, :) = getRankFractional(sample(idim, :))
                end do
            else
                do idim = 1, ndim
                    frank_ref(:, idim) = getRankFractional(sample(:, idim))
                end do
            end if
            if (present(weight)) then
                rho = getCor(frank_ref, dim, weight)
            else
                rho = getCor(frank_ref, dim)
            end if
        end function

        subroutine report(line, subsetr, this)
            integer, intent(in) :: line
            class(*), intent(in) :: subsetr
            character(*, SK), intent(in) :: this
            if (same_type_as(subsetr, uppDia)) then
                call setMatCopy(rho, rdpack, rho, rdpack, uppDia, transHerm)
            elseif (same_type_as(subsetr, lowDia)) then
                call setMatCopy(rho, rdpack, rho, rdpack, lowDia, transHerm)
            elseif (same_type_as(subsetr, uppLowDia)) then
                rho(2,1) = rho(1,2)
            else
                error stop "Unrecognized subset."
            end if
            rdiff = abs(rho - rho_ref)
            assertion = assertion .and. all(rdiff < rtol)
            assertion = assertion .and. all(abs(frank - frank_ref) < rtol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[ndim, nsam, dim, GET_SHAPE(sample, IK)]")
                call test%disp%show( [ndim, nsam, dim, GET_SHAPE(sample, IK)] )
                call test%disp%show("sample")
                call test%disp%show( sample )
                call test%disp%show("iweight")
                call test%disp%show( iweight )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("frank_ref")
                call test%disp%show( frank_ref )
                call test%disp%show("frank")
                call test%disp%show( frank )
                call test%disp%show("rho_ref")
                call test%disp%show( rho_ref )
                call test%disp%show("rho")
                call test%disp%show( rho )
                call test%disp%show("rdiff")
                call test%disp%show( rdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            ! two tests, two reports.
            call test%assert(assertion, SK_"The `rho` and fractional ranks must be computed correctly for "//this//SK_" sample.", int(line, IK))
            call test%assert(assertion, SK_"The `rho` and fractional ranks must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%
#elif   setCordance_ENABLED
        !%%%%%%%%%%%%%%%%%%

        assertion = .true.

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_SAMPLE
#undef  GET_CONJG
#undef  GET_SHAPE