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
!>  This file contains procedure implementations of tests of [pm_sampleCov](@ref pm_sampleCov).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC), parameter :: rtol = epsilon(1._TKC) * 100
        ! Define the sample type.
#if     CK_ENABLED
#define GET_CONJG(X)conjg(X)
#define TYPE_OF_SAMPLE complex(TKC)
        complex(TKC), parameter :: ZERO = (0._TKC, 0._TKC), ONE = (1._TKC, 1._TKC), TWO = 2 * (1._TKC, 1._TKC), ctol = (rtol, rtol), RONE = (1._TKC, 0._TKC)
#elif   RK_ENABLED
#define GET_CONJG(X)X
#define TYPE_OF_SAMPLE real(TKC)
        real(TKC), parameter :: ZERO = 0._TKC, ONE = 1._TKC, RONE = 1._TKC, TWO = 2._TKC, ctol = rtol
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%
#if     getCov_ENABLED
        !%%%%%%%%%%%%%

        real(TKC) :: rweisqs
        real(TKC) :: rweisum
        integer(IK) :: iweisum
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam, ndim, dim
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:), meang(:)
        TYPE_OF_SAMPLE, allocatable :: cov(:,:), cor(:,:), cov_ref(:,:), cdiff(:,:)
        real(TKC), allocatable :: std(:)
        assertion = .true._LK

        do itry = 1, 50

            ! test CFC subsetr = lowDia.

            block

                type(lowDia_type), parameter :: subsetr = lowDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                cov_ref = getFilled(ZERO, ndim, ndim)
                std = getUnifRand(1._TKC, 2._TKC, ndim)
                cor = getUnifRand(-ONE, ONE, ndim, ndim)
                call setMatInit(cor, getSubSymm(subsetr), ZERO, RONE)
                call setCov(cov_ref, subsetr, cor, subsetr, std)
                call setMatCopy(cov_ref, rdpack, cov_ref, rdpack, subsetr, transHerm)

                cov = getCov(cor, subsetr, std)
                call reportCFC(__LINE__)

            end block

            ! test CFC subsetr = uppDia.

            block

                type(uppDia_type), parameter :: subsetr = uppDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                cov_ref = getFilled(ZERO, ndim, ndim)
                std = getUnifRand(1._TKC, 2._TKC, ndim)
                cor = getUnifRand(-ONE, ONE, ndim, ndim)
                call setMatInit(cor, getSubSymm(subsetr), ZERO, RONE)
                call setCov(cov_ref, subsetr, cor, subsetr, std)
                call setMatCopy(cov_ref, rdpack, cov_ref, rdpack, subsetr, transHerm)

                cov = getCov(cor, subsetr, std)
                call reportCFC(__LINE__)

            end block

            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample ULD interface.

            block

                type(lowDia_type), parameter :: subset = lowDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                cov_ref = getFilled(ZERO, ndim, ndim)
                mean = getFilled(ZERO, ndim)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                    meang = sample(:,1)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                    meang = sample(1,:)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                rweisqs = sum(rweight**2)

                ! integer weighted

                cov = getCov(sample, dim, iweight)
                call setCovMean(cov_ref, subset, mean, sample, dim, iweight, iweisum, meang)
                call setMatCopy(cov_ref, rdpack, cov_ref, rdpack, subset, transHerm)
                call report(__LINE__, SK_"integer-weighted")

                cov = getCov(sample, dim, iweight, fweight_type())
                cov_ref = cov_ref * getVarCorrection(real(iweisum, TKC))
                call report(__LINE__, SK_"integer-weighted")

                ! real weighted

                cov = getCov(sample, dim, rweight)
                call setCovMean(cov_ref, subset, mean, sample, dim, rweight, rweisum, meang)
                call setMatCopy(cov_ref, rdpack, cov_ref, rdpack, subset, transHerm)
                call report(__LINE__, SK_"real-weighted")

                cov = getCov(sample, dim, rweight, rweight_type())
                cov_ref = cov_ref * getVarCorrection(rweisum, rweisqs)
                call report(__LINE__, SK_"real-weighted")

                ! unweighted

                cov = getCov(sample, dim)
                call setCovMean(cov_ref, subset, mean, sample, dim, meang)
                call setMatCopy(cov_ref, rdpack, cov_ref, rdpack, subset, transHerm)
                call report(__LINE__, SK_"unweighted")

                cov = getCov(sample, dim, rweight_type())
                cov_ref = cov_ref * getVarCorrection(real(nsam, TKC))
                call report(__LINE__, SK_"unweighted")

            end block

            ! test sample XY interface.

            block

                type(lowDia_type), parameter :: subset = lowDia_type()
                ndim = 2_IK
                mean = getFilled(ZERO, ndim)
                cov_ref = getFilled(ZERO, ndim, ndim)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                rweisqs = sum(rweight**2)
                meang = sample(1,:)

                ! integer weighted

                cov = getCov(sample(:,1), sample(:,2), iweight)
                call setCovMean(cov_ref, mean, sample(:,1), sample(:,2), iweight, iweisum, meang)
                call report(__LINE__, SK_"integer-weighted")

                cov = getCov(sample(:,1), sample(:,2), iweight, fweight_type())
                cov_ref = cov_ref * getVarCorrection(real(iweisum, TKC))
                call report(__LINE__, SK_"integer-weighted")

                ! real weighted

                cov = getCov(sample(:,1), sample(:,2), rweight)
                call setCovMean(cov_ref, mean, sample(:,1), sample(:,2), rweight, rweisum, meang)
                call report(__LINE__, SK_"real-weighted")

                cov = getCov(sample(:,1), sample(:,2), rweight, rweight_type())
                cov_ref = cov_ref * getVarCorrection(rweisum, rweisqs)
                call report(__LINE__, SK_"real-weighted")

                ! unweighted

                cov = getCov(sample(:,1), sample(:,2))
                call setCovMean(cov_ref, mean, sample(:,1), sample(:,2), meang)
                call report(__LINE__, SK_"unweighted")

                cov = getCov(sample(:,1), sample(:,2), rweight_type())
                cov_ref = cov_ref * getVarCorrection(real(nsam, TKC))
                call report(__LINE__, SK_"unweighted")

            end block

        end do

    contains

        subroutine setAssertedCov(cov)
            TYPE_OF_SAMPLE, intent(inout) :: cov(:,:)
            cdiff = abs(cov - cov_ref)
            assertion = assertion .and. all(cdiff < ctol)
        end subroutine

        subroutine reportCFC(line)
            integer, intent(in) :: line
            call setAssertedCov(cov)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("ndim")
                call test%disp%show( ndim )
                call test%disp%show("std")
                call test%disp%show( std )
                call test%disp%show("cor")
                call test%disp%show( cor )
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cov_ref")
                call test%disp%show( cov_ref )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cov` must be correctly computed from the specified `cor` and `std`.", int(line, IK))
        end subroutine

        subroutine report(line, this)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: this
            call setAssertedCov(cov)
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
                call test%disp%show("cov_ref")
                call test%disp%show( cov_ref )
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cov` must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%
#elif   setCov_ENABLED
        !%%%%%%%%%%%%%

        real(TKC) :: rweisum
        integer(IK) :: iweisum
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam, ndim, dim
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:)
        TYPE_OF_SAMPLE, allocatable :: cov(:,:), cor(:,:), cov_ref(:,:), cdiff(:,:)
        real(TKC), allocatable :: std(:)
        assertion = .true._LK

        do itry = 1, 50

            ! test CFC subsetr = lowDia.

            block

                type(lowDia_type), parameter :: subsetr = lowDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                std = getUnifRand(1._TKC, 2._TKC, ndim)
                cor = getUnifRand(-ONE, ONE, ndim, ndim)
                call setMatInit(cor, getSubSymm(subsetr), ZERO, RONE)

                cov_ref = getCFC(cor, subsetr, std)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subsetr, cor, subsetr, std)
                call reportCFC(__LINE__, subsetr)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, getSubSymm(subsetr), cor, subsetr, std)
                call reportCFC(__LINE__, getSubSymm(subsetr))

            end block

            ! test CFC subsetr = uppDia.

            block

                type(uppDia_type), parameter :: subsetr = uppDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                std = getUnifRand(1._TKC, 2._TKC, ndim)
                cor = getUnifRand(-ONE, ONE, ndim, ndim)
                call setMatInit(cor, getSubSymm(subsetr), ZERO, RONE)

                cov_ref = getCFC(cor, subsetr, std)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subsetr, cor, subsetr, std)
                call reportCFC(__LINE__, subsetr)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, getSubSymm(subsetr), cor, subsetr, std)
                call reportCFC(__LINE__, getSubSymm(subsetr))

            end block

            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample lowDia interface.

            block

                type(lowDia_type), parameter :: subset = lowDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                cov_ref = getCovRef(mean, sample, dim, real(iweight, TKC))

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, mean, sample, dim, iweight, iweisum)
                call report(__LINE__, subset, SK_"integer-weighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, getShifted(sample, dim, -mean), dim, iweight, iweisum)
                call report(__LINE__, subset, SK_"integer-weighted shifted")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                cov_ref = getCovRef(mean, sample, dim, rweight)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, mean, sample, dim, rweight, rweisum)
                call report(__LINE__, subset, SK_"real-weighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, getShifted(sample, dim, -mean), dim, rweight, rweisum)
                call report(__LINE__, subset, SK_"real-weighted shifted")

                ! unweighted

                mean = getMean(sample, dim)
                cov_ref = getCovRef(mean, sample, dim)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, mean, sample, dim)
                call report(__LINE__, subset, SK_"unweighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, getShifted(sample, dim, -mean), dim)
                call report(__LINE__, subset, SK_"unweighted shifted")

            end block

            ! test sample uppDia interface.

            block

                type(uppDia_type), parameter :: subset = uppDia_type()
                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                cov_ref = getCovRef(mean, sample, dim, real(iweight, TKC))

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, mean, sample, dim, iweight, iweisum)
                call report(__LINE__, subset, SK_"integer-weighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, getShifted(sample, dim, -mean), dim, iweight, iweisum)
                call report(__LINE__, subset, SK_"integer-weighted shifted")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                cov_ref = getCovRef(mean, sample, dim, rweight)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, mean, sample, dim, rweight, rweisum)
                call report(__LINE__, subset, SK_"real-weighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, getShifted(sample, dim, -mean), dim, rweight, rweisum)
                call report(__LINE__, subset, SK_"real-weighted shifted")

                ! unweighted

                mean = getMean(sample, dim)
                cov_ref = getCovRef(mean, sample, dim)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, mean, sample, dim)
                call report(__LINE__, subset, SK_"unweighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, subset, getShifted(sample, dim, -mean), dim)
                call report(__LINE__, subset, SK_"unweighted shifted")

            end block

            ! test sample XY interface.

            block

                type(lowDia_type), parameter :: subset = lowDia_type()
                dim = 1_IK
                ndim = 2_IK
                mean = getFilled(ZERO, ndim)
                cov_ref = getFilled(ZERO, ndim, ndim)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                cov_ref = getCovRef(mean, sample, dim, real(iweight, TKC))

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, mean, sample(:,1), sample(:,2), iweight, iweisum)
                call report(__LINE__, uppLowDia, SK_"integer-weighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, sample(:,1) - mean(1), sample(:,2) - mean(2), iweight, iweisum)
                call report(__LINE__, uppLowDia, SK_"integer-weighted")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                cov_ref = getCovRef(mean, sample, dim, rweight)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, mean, sample(:,1), sample(:,2), rweight, rweisum)
                call report(__LINE__, uppLowDia, SK_"real-weighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, sample(:,1) - mean(1), sample(:,2) - mean(2), rweight, rweisum)
                call report(__LINE__, uppLowDia, SK_"real-weighted")

                ! unweighted

                mean = getMean(sample, dim)
                cov_ref = getCovRef(mean, sample, dim)

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, mean, sample(:,1), sample(:,2))
                call report(__LINE__, uppLowDia, SK_"unweighted")

                cov = getFilled(ZERO, ndim, ndim)
                call setCov(cov, sample(:,1) - mean(1), sample(:,2) - mean(2))
                call report(__LINE__, uppLowDia, SK_"unweighted")

            end block

        end do

    contains

        subroutine setAssertedCov(cov, subset)
            TYPE_OF_SAMPLE, intent(inout) :: cov(:,:)
            class(*), intent(in) :: subset
            select type(subset)
            type is (lowDia_type)
                call setMatCopy(cov(1:,2:), rdpack, cov(2:,1:), rdpack, subset, transHerm)
            type is (uppDia_type)
                call setMatCopy(cov(2:,1:), rdpack, cov(1:,2:), rdpack, subset, transHerm)
            type is (uppLowDia_type)
                continue
            class default
                error stop "Unrecognized subset." ! LCOV_EXCL_LINE
            end select
            cdiff = abs(cov - cov_ref)
            assertion = assertion .and. all(cdiff < ctol)
        end subroutine

        pure function getCFC(cor, subsetr, std) result(cov)
            class(*), intent(in) :: subsetr
            real(TKC), intent(in) :: std(:)
            TYPE_OF_SAMPLE, intent(in) :: cor(:,:)
            TYPE_OF_SAMPLE :: cov(size(cor, 1, IK), size(cor, 2, IK))
            integer(IK) :: idim, jdim, ndim
            ndim = size(std, 1, IK)
            do jdim = 1, ndim
                do idim = jdim, ndim
                    if (same_type_as(subsetr, lowDia)) then
                        cov(idim, jdim) = cor(idim, jdim) * std(idim) * std(jdim)
                        cov(jdim, idim) = GET_CONJG(cov(idim, jdim))
                    elseif (same_type_as(subsetr, uppDia)) then
                        cov(jdim, idim) = cor(jdim, idim) * std(idim) * std(jdim)
                        cov(idim, jdim) = GET_CONJG(cov(jdim, idim))
                    else
                        error stop "Unrecognized subsetr." ! LCOV_EXCL_LINE
                    end if
                end do
            end do
        end function

        subroutine reportCFC(line, subset)
            class(*), intent(in) :: subset
            integer, intent(in) :: line
            call setAssertedCov(cov, subset)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("ndim")
                call test%disp%show( ndim )
                call test%disp%show("std")
                call test%disp%show( std )
                call test%disp%show("cor")
                call test%disp%show( cor )
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cov_ref")
                call test%disp%show( cov_ref )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cov` must be correctly computed from the specified `cor` and `std`.", int(line, IK))
        end subroutine

        PURE function getCovRef(mean, sample, dim, weight) result(cov)
            integer(IK), intent(in) :: dim
            real(TKC), intent(in), optional :: weight(:)
            TYPE_OF_SAMPLE, intent(in), contiguous :: sample(:,:), mean(:)
            TYPE_OF_SAMPLE :: cov(size(sample, 3 - dim, IK), size(sample, 3 - dim, IK))
            TYPE_OF_SAMPLE :: sampleShifted(size(sample, 1, IK), size(sample, 2, IK))
            integer(IK) :: idim, jdim, ndim
            real(TKC) :: normfac
            ndim = size(sample, 3 - dim, IK)
            sampleShifted = getShifted(sample, dim, -mean)
            if (present(weight)) then
                normfac = 1._TKC / sum(weight)
                if (dim == 2) then
                    do jdim = 1, ndim
                        do idim = 1, ndim
                            cov(idim, jdim) = dot_product(sampleShifted(jdim,:), sampleShifted(idim,:) * weight) * normfac
                        end do
                    end do
                else
                    do jdim = 1, ndim
                        do idim = 1, ndim
                            cov(idim, jdim) = dot_product(sampleShifted(:,jdim), sampleShifted(:,idim) * weight) * normfac
                        end do
                    end do
                end if
            else
                normfac = 1._TKC / size(sample, dim, IK)
                if (dim == 2) then
                    do jdim = 1, ndim
                        do idim = 1, ndim
                            cov(idim, jdim) = dot_product(sampleShifted(jdim,:), sampleShifted(idim,:)) * normfac
                        end do
                    end do
                else
                    do jdim = 1, ndim
                        do idim = 1, ndim
                            cov(idim, jdim) = dot_product(sampleShifted(:,jdim), sampleShifted(:,idim)) * normfac
                        end do
                    end do
                end if
            end if
        end function

        subroutine report(line, subset, this)
            integer, intent(in) :: line
            class(*), intent(in) :: subset
            character(*, SK), intent(in) :: this
            call setAssertedCov(cov, subset)
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
                call test%disp%show("cov_ref")
                call test%disp%show( cov_ref )
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `cov` must be computed correctly for "//this//SK_" sample.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%
#elif   setCovMean_ENABLED
        !%%%%%%%%%%%%%%%%%

        real(TKC) :: rweisum, rweisum_ref
        integer(IK) :: iweisum, iweisum_ref
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam, ndim, dim
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:), meang(:), mean_ref(:), mdiff(:)
        TYPE_OF_SAMPLE, allocatable :: cov(:,:), cov_ref(:,:), cdiff(:,:)
        assertion = .true._LK

        do itry = 1, 50

            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample lowDia interface.

            block

                type(lowDia_type), parameter :: subset = lowDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                dim = getChoice([1, 2])
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                    meang = sample(:,1)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                    meang = sample(1,:)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                mean_ref = getFilled(ZERO, ndim)
                cov_ref = getFilled(ZERO, ndim, ndim)

                ! integer weighted

                mean_ref = getMean(sample, dim, iweight)
                call setCov(cov_ref, subset, mean_ref, sample, dim, iweight, iweisum_ref)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, subset, mean, sample, dim, iweight, iweisum, meang)
                call setAssertedSum(__LINE__, SK_"integer-weighted", real(iweisum, TKC), real(iweisum_ref, TKC))
                call setAssertedCov(__LINE__, SK_"integer-weighted", subset)
                call setAssertedAvg(__LINE__, SK_"integer-weighted")

                ! real weighted

                mean_ref = getMean(sample, dim, rweight)
                call setCov(cov_ref, subset, mean_ref, sample, dim, rweight, rweisum_ref)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, subset, mean, sample, dim, rweight, rweisum, meang)
                call setAssertedSum(__LINE__, SK_"real-weighted", rweisum, rweisum_ref)
                call setAssertedCov(__LINE__, SK_"real-weighted", subset)
                call setAssertedAvg(__LINE__, SK_"real-weighted")

                ! unweighted

                mean_ref = getMean(sample, dim)
                call setCov(cov_ref, subset, mean_ref, sample, dim)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, subset, mean, sample, dim, meang)
                call setAssertedCov(__LINE__, SK_"unweighted", subset)
                call setAssertedAvg(__LINE__, SK_"unweighted")

            end block

            ! test sample uppDia interface.

            block

                type(uppDia_type), parameter :: subset = uppDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                dim = getChoice([1, 2])
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                    meang = sample(:,1)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                    meang = sample(1,:)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                mean_ref = getFilled(ZERO, ndim)
                cov_ref = getFilled(ZERO, ndim, ndim)

                ! integer weighted

                mean_ref = getMean(sample, dim, iweight)
                call setCov(cov_ref, subset, mean_ref, sample, dim, iweight, iweisum_ref)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, subset, mean, sample, dim, iweight, iweisum, meang)
                call setAssertedSum(__LINE__, SK_"integer-weighted", real(iweisum, TKC), real(iweisum_ref, TKC))
                call setAssertedCov(__LINE__, SK_"integer-weighted", subset)
                call setAssertedAvg(__LINE__, SK_"integer-weighted")

                ! real weighted

                mean_ref = getMean(sample, dim, rweight)
                call setCov(cov_ref, subset, mean_ref, sample, dim, rweight, rweisum_ref)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, subset, mean, sample, dim, rweight, rweisum, meang)
                call setAssertedSum(__LINE__, SK_"real-weighted", rweisum, rweisum_ref)
                call setAssertedCov(__LINE__, SK_"real-weighted", subset)
                call setAssertedAvg(__LINE__, SK_"real-weighted")

                ! unweighted

                mean_ref = getMean(sample, dim)
                call setCov(cov_ref, subset, mean_ref, sample, dim)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, subset, mean, sample, dim, meang)
                call setAssertedCov(__LINE__, SK_"unweighted", subset)
                call setAssertedAvg(__LINE__, SK_"unweighted")

            end block

            ! test sample XY interface.

            block

                type(uppLowDia_type), parameter :: subset = uppLowDia_type()
                dim = 1_IK
                ndim = 2_IK
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                meang = sample(1,:)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                mean_ref = getFilled(ZERO, ndim)
                cov_ref = getFilled(ZERO, ndim, ndim)

                ! integer weighted

                mean_ref = getMean(sample, dim, iweight)
                call setCov(cov_ref, mean_ref, sample(:,1), sample(:,2), iweight, iweisum_ref)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, mean, sample(:,1), sample(:,2), iweight, iweisum, meang)
                call setAssertedSum(__LINE__, SK_"integer-weighted", real(iweisum, TKC), real(iweisum_ref, TKC))
                call setAssertedCov(__LINE__, SK_"integer-weighted", subset)
                call setAssertedAvg(__LINE__, SK_"integer-weighted")

                ! real weighted

                mean_ref = getMean(sample, dim, rweight)
                call setCov(cov_ref, mean_ref, sample(:,1), sample(:,2), rweight, rweisum_ref)

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, mean, sample(:,1), sample(:,2), rweight, rweisum, meang)
                call setAssertedSum(__LINE__, SK_"real-weighted", rweisum, rweisum_ref)
                call setAssertedCov(__LINE__, SK_"real-weighted", subset)
                call setAssertedAvg(__LINE__, SK_"real-weighted")

                ! unweighted

                mean_ref = getMean(sample, dim)
                call setCov(cov_ref, mean_ref, sample(:,1), sample(:,2))

                mean = getFilled(ZERO, ndim)
                cov = getFilled(ZERO, ndim, ndim)
                call setCovMean(cov, mean, sample(:,1), sample(:,2), meang)
                call setAssertedCov(__LINE__, SK_"unweighted", subset)
                call setAssertedAvg(__LINE__, SK_"unweighted")

            end block

        end do

    contains

        subroutine setAssertedSum(line, this, weisum, weisum_ref)
            real(TKC), intent(in) :: weisum, weisum_ref
            character(*, SK), intent(in) :: this
            integer, intent(in) :: line
            assertion = assertion .and. abs(weisum - weisum_ref) < rtol * weisum_ref
            call report(line, SK_"The `weisum` must be computed correctly for "//this//SK_" sample.")
        end subroutine

        subroutine setAssertedAvg(line, this)
            character(*, SK), intent(in) :: this
            integer, intent(in) :: line
            mdiff = abs(mean - mean_ref)
            assertion = assertion .and. all(mdiff < ctol)
            call report(line, SK_"The `mean` must be computed correctly for "//this//SK_" sample.")
        end subroutine

        subroutine setAssertedCov(line, this, subset)
            character(*, SK), intent(in) :: this
            class(*), intent(in) :: subset
            integer, intent(in) :: line
            select type(subset)
            type is (lowDia_type)
                call setMatCopy(cov(1:,2:), rdpack, cov(2:,1:), rdpack, subset, transHerm)
                call setMatCopy(cov_ref(1:,2:), rdpack, cov_ref(2:,1:), rdpack, subset, transHerm)
            type is (uppDia_type)
                call setMatCopy(cov(2:,1:), rdpack, cov(1:,2:), rdpack, subset, transHerm)
                call setMatCopy(cov_ref(2:,1:), rdpack, cov_ref(1:,2:), rdpack, subset, transHerm)
            type is (uppLowDia_type)
                continue
            class default
                error stop "Unrecognized subset." ! LCOV_EXCL_LINE
            end select
            cdiff = abs(cov - cov_ref)
            assertion = assertion .and. all(cdiff < ctol)
            call report(line, SK_"The `cov` must be computed correctly for "//this//SK_" sample.")
        end subroutine

        subroutine report(line, msg)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: msg
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
                call test%disp%show("cov_ref")
                call test%disp%show( cov_ref )
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%show("mean_ref")
                call test%disp%show( mean_ref )
                call test%disp%show("mean")
                call test%disp%show( mean )
                call test%disp%show("mdiff")
                call test%disp%show( mdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCovMerged_ENABLED || setCovMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: fracA
        integer(IK) :: itry, ndim, dim
        integer(IK) :: nsam, nsamA, nsamB
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), meanDiff(:)
        TYPE_OF_SAMPLE, allocatable :: covA(:,:), covB(:,:), cov(:,:), cov_ref(:,:), cdiff(:,:)
        assertion = .true._LK
        dim = 2_IK

        do itry = 1, 50

            nsamA = getUnifRand(2_IK, 5_IK)
            nsamB = getUnifRand(2_IK, 5_IK)
            nsam = nsamA + nsamB
            fracA = real(nsamA, TKC) / real(nsam, TKC)

            ! test D2 interface.

            block

                type(uppDia_type), parameter :: subset = uppDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                cov_ref = getCov(sample, dim)
                covA = getCov(sample(:,1:nsamA), dim)
                covB = getCov(sample(:,nsamA+1:), dim)
                meanDiff = getMean(sample(:,1:nsamA), dim) - getMean(sample(:,nsamA+1:), dim)

#if             getCovMerged_ENABLED
                cov = getFilled(ZERO, ndim, ndim)
                cov = getCovMerged(covB, covA, meanDiff, fracA)
                cdiff = abs(cov - cov_ref)
                assertion = assertion .and. all(cdiff < ctol)
                call report(__LINE__, SK_"The new `cov` must be computed correctly.")
#elif           setCovMerged_ENABLED

                ! new

                cov = getFilled(ZERO, ndim, ndim)
                call setCovMerged(cov, covB, covA, meanDiff, fracA, subset)
                call setAssertedCov(__LINE__, subset, SK_"new")

                cov = getFilled(ZERO, ndim, ndim)
                call setCovMerged(cov, covB, covA, meanDiff, fracA, getSubSymm(subset))
                call setAssertedCov(__LINE__, getSubSymm(subset), SK_"new")

                ! old

                cov = covB
                call setCovMerged(cov, covA, meanDiff, fracA, subset)
                call setAssertedCov(__LINE__, subset, SK_"in-place")

                cov = covB
                call setCovMerged(cov, covA, meanDiff, fracA, getSubSymm(subset))
                call setAssertedCov(__LINE__, getSubSymm(subset), SK_"in-place")

#endif

            end block

        end do

    contains

#if     setCovMerged_ENABLED
        subroutine setAssertedCov(line, subset, specific)
            character(*, SK), intent(in) :: specific
            class(*), intent(in) :: subset
            integer, intent(in) :: line
            select type (subset)
            type is (lowDia_type)
                call setMatCopy(cov, rdpack, cov, rdpack, subset, transHerm)
            type is (uppDia_type)
                call setMatCopy(cov, rdpack, cov, rdpack, subset, transHerm)
            class default
                error stop "Unrecognized subset." ! LCOV_EXCL_LINE
            end select
            cdiff = abs(cov - cov_ref)
            assertion = assertion .and. all(cdiff < ctol)
            call report(line, SK_"The "//specific//SK_" `cov` must be computed correctly.")
        end subroutine
#endif

        subroutine report(line, msg)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: msg
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[ndim, nsamA, nsamB, nsam, dim, shape(sample, IK)]")
                call test%disp%show( [ndim, nsamA, nsamB, nsam, dim, shape(sample, IK)] )
                call test%disp%show("sample")
                call test%disp%show( sample )
                call test%disp%show("meanDiff")
                call test%disp%show( meanDiff )
                call test%disp%show("covA")
                call test%disp%show( covA )
                call test%disp%show("covB")
                call test%disp%show( covB )
                call test%disp%show("cov_ref")
                call test%disp%show( cov_ref )
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCovMeanMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: fracA
        integer(IK) :: itry, ndim, dim
        integer(IK) :: nsam, nsamA, nsamB
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:), mean_ref(:), meanA(:), meanB(:), mdiff(:)
        TYPE_OF_SAMPLE, allocatable :: covA_ref(:,:), covB_ref(:,:), covA(:,:), covB(:,:), cov(:,:), cov_ref(:,:), cdiff(:,:)
        assertion = .true._LK
        dim = 2_IK

        do itry = 1, 50

            nsamA = getUnifRand(2_IK, 5_IK)
            nsamB = getUnifRand(2_IK, 5_IK)
            nsam = nsamA + nsamB
            fracA = real(nsamA, TKC) / real(nsam, TKC)

            ! test D2 interface.

            block

                type(uppDia_type), parameter :: subset = uppDia_type()
                ndim = getUnifRand(1_IK, 5_IK)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                cov_ref = getCov(sample, dim)
                mean_ref = getMean(sample, dim)
                covA_ref = getCov(sample(:,1:nsamA), dim)
                covB_ref = getCov(sample(:,nsamA+1:), dim)
                meanA = getMean(sample(:,1:nsamA), dim)
                meanB = getMean(sample(:,nsamA+1:), dim)

                ! new subset

                covA = getFilled(ZERO, ndim, ndim); call setMatCopy(covA, rdpack, covA_ref, rdpack, subset)
                covB = getFilled(ZERO, ndim, ndim); call setMatCopy(covB, rdpack, covB_ref, rdpack, subset)
                cov = getFilled(ZERO, ndim, ndim)
                mean = getFilled(ZERO, ndim)

                call setCovMeanMerged(cov, mean, covB, meanB, covA, meanA, fracA, subset)
                call setAssertedCov(__LINE__, subset, SK_"new")
                call setAssertedAvg(__LINE__, SK_"new")

                ! new subset opposite

                covA = getFilled(ZERO, ndim, ndim); call setMatCopy(covA, rdpack, covA_ref, rdpack, getSubSymm(subset))
                covB = getFilled(ZERO, ndim, ndim); call setMatCopy(covB, rdpack, covB_ref, rdpack, getSubSymm(subset))
                cov = getFilled(ZERO, ndim, ndim)
                mean = getFilled(ZERO, ndim)

                call setCovMeanMerged(cov, mean, covB, meanB, covA, meanA, fracA, getSubSymm(subset))
                call setAssertedCov(__LINE__, getSubSymm(subset), SK_"new")
                call setAssertedAvg(__LINE__, SK_"new")

                ! old subset

                covA = getFilled(ZERO, ndim, ndim); call setMatCopy(covA, rdpack, covA_ref, rdpack, subset)
                covB = getFilled(ZERO, ndim, ndim); call setMatCopy(covB, rdpack, covB_ref, rdpack, subset)
                mean = meanB
                cov = covB

                call setCovMeanMerged(cov, mean, covA, meanA, fracA, subset)
                call setAssertedCov(__LINE__, subset, SK_"in-place")
                call setAssertedAvg(__LINE__, SK_"in-place")

                ! old subset opposite

                covA = getFilled(ZERO, ndim, ndim); call setMatCopy(covA, rdpack, covA_ref, rdpack, getSubSymm(subset))
                covB = getFilled(ZERO, ndim, ndim); call setMatCopy(covB, rdpack, covB_ref, rdpack, getSubSymm(subset))
                mean = meanB
                cov = covB

                call setCovMeanMerged(cov, mean, covA, meanA, fracA, getSubSymm(subset))
                call setAssertedCov(__LINE__, getSubSymm(subset), SK_"in-place")
                call setAssertedAvg(__LINE__, SK_"in-place")

            end block

        end do

    contains

        subroutine setAssertedAvg(line, specific)
            character(*, SK), intent(in) :: specific
            integer, intent(in) :: line
            mdiff = abs(mean - mean_ref)
            assertion = assertion .and. all(mdiff < ctol)
            call report(line, SK_"The "//specific//SK_" `meanMerged` must be computed correctly.")
        end subroutine

        subroutine setAssertedCov(line, subset, specific)
            character(*, SK), intent(in) :: specific
            class(*), intent(in) :: subset
            integer, intent(in) :: line
            select type (subset)
            type is (lowDia_type)
                call setMatCopy(cov, rdpack, cov, rdpack, subset, transHerm)
            type is (uppDia_type)
                call setMatCopy(cov, rdpack, cov, rdpack, subset, transHerm)
            class default
                error stop "Unrecognized subset." ! LCOV_EXCL_LINE
            end select
            cdiff = abs(cov - cov_ref)
            assertion = assertion .and. all(cdiff < ctol)
            call report(line, SK_"The "//specific//SK_" `cov` must be computed correctly.")
        end subroutine

        subroutine report(line, msg)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: msg
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[ndim, nsamA, nsamB, nsam, dim, shape(sample, IK)]")
                call test%disp%show( [ndim, nsamA, nsamB, nsam, dim, shape(sample, IK)] )
                call test%disp%show("sample")
                call test%disp%show( sample )
                call test%disp%show("meanA")
                call test%disp%show( meanA )
                call test%disp%show("meanB")
                call test%disp%show( meanB )
                call test%disp%show("mean_ref")
                call test%disp%show( mean_ref )
                call test%disp%show("mean")
                call test%disp%show( mean )
                call test%disp%show("mdiff")
                call test%disp%show( mdiff )
                call test%disp%show("covA")
                call test%disp%show( covA )
                call test%disp%show("covB")
                call test%disp%show( covB )
                call test%disp%show("cov_ref")
                call test%disp%show( cov_ref )
                call test%disp%show("cov")
                call test%disp%show( cov )
                call test%disp%show("cdiff")
                call test%disp%show( cdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_SAMPLE
#undef  GET_CONJG