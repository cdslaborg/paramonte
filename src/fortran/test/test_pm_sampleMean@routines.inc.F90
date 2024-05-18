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
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC), parameter :: rtol = epsilon(1._TKC) * 100
        ! Define the sample type.
#if     CK_ENABLED
#define TYPE_OF_SAMPLE complex(TKC)
        complex(TKC), parameter :: ONE = (1._TKC, 1._TKC), TWO = 2 * (1._TKC, 1._TKC), ctol = (rtol, rtol)
#elif   RK_ENABLED
#define TYPE_OF_SAMPLE real(TKC)
        real(TKC), parameter :: ONE = 1._TKC, TWO = 2._TKC, ctol = rtol
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%
#if     getMean_ENABLED
        !%%%%%%%%%%%%%%

        real(TKC) :: rweisum
        integer(IK) :: iweisum
        integer(IK) :: itry, nsam, ndim, dim
        real(TKC), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:), mean_ref(:), diff(:)
        assertion = .true._LK

        do itry = 1, 50

            dim = 0
            nsam = getUnifRand(1_IK, 5_IK)

            ! test sample D2 DIM interface.

            block

                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(mean_ref, ndim)
                call setResized(mean, ndim)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                call setMean(mean_ref, sample, dim, iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                call setMean(mean_ref, sample, dim, rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean = getMean(sample, dim)
                call setMean(mean_ref, sample, dim)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

            ! test sample D2 ALL interface.

            block

                dim = 0
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(mean_ref, 1_IK)
                call setResized(mean, 1_IK)
                if (getUnifRand()) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam * ndim)
                iweight = getUnifRand(1_IK, 9_IK, nsam * ndim)

                ! integer weighted

                mean(1) = getMean(sample, iweight)
                call setMean(mean_ref(1), sample, iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean(1) = getMean(sample, rweight)
                call setMean(mean_ref(1), sample, rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean(1) = getMean(sample)
                call setMean(mean_ref(1), sample)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

            ! test sample D1 DIM interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(mean, 1_IK)
                call setResized(mean_ref, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)

                ! integer weighted

                mean(1) = getMean(sample(:,1), dim, iweight)
                call setMean(mean_ref(1), sample(:,1), dim, iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean(1) = getMean(sample(:,1), dim, rweight)
                call setMean(mean_ref(1), sample(:,1), dim, rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean(1) = getMean(sample(:,1), dim)
                call setMean(mean_ref(1), sample(:,1), dim)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

            ! test sample D1 ALL interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(mean, 1_IK)
                call setResized(mean_ref, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)

                ! integer weighted

                mean(1) = getMean(sample(:,1), iweight)
                call setMean(mean_ref(1), sample(:,1), iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean(1) = getMean(sample(:,1), rweight)
                call setMean(mean_ref(1), sample(:,1), rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean(1) = getMean(sample(:,1))
                call setMean(mean_ref(1), sample(:,1))

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

        end do

    contains

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
                call test%disp%show("mean_ref")
                call test%disp%show( mean_ref )
                call test%disp%show("mean")
                call test%disp%show( mean )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%
#elif   setMean_ENABLED
        !%%%%%%%%%%%%%%

        integer(IK) :: itry, nsam, ndim, dim
        real(TKC) :: rweisum, rweisum_ref
        real(TKC), allocatable :: rweight(:)
        integer(IK) :: iweisum, iweisum_ref
        integer(IK), allocatable :: iweight(:)
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:), mean_ref(:), diff(:)
        assertion = .true._LK

        do itry = 1, 50

            dim = 0
            nsam = getUnifRand(1_IK, 5_IK)

            ! test XY interface.

            block

                ndim = 2_IK
                call setResized(mean, ndim)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                mean_ref = getMeanD2(sample, 1_IK, real(iweight, TKC))
                call setMean(mean, sample(:,1), sample(:,2), iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                assertion = assertion .and. iweisum == iweisum_ref
                call report(__LINE__, SK_"The `weisum` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean_ref = getMeanD2(sample, 1_IK, rweight)
                call setMean(mean, sample(:,1), sample(:,2), rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                assertion = assertion .and. abs(rweisum - rweisum_ref) <= rtol * abs(rweisum_ref)
                call report(__LINE__, SK_"The `weisum` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean_ref = getMeanD2(sample, 1_IK)
                call setMean(mean, sample(:,1), sample(:,2))

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

            ! test sample D2 DIM interface.

            block

                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(mean, ndim)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                mean_ref = getMeanD2(sample, dim, real(iweight, TKC))
                call setMean(mean, sample, dim, iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                assertion = assertion .and. iweisum == iweisum_ref
                call report(__LINE__, SK_"The `weisum` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean_ref = getMeanD2(sample, dim, rweight)
                call setMean(mean, sample, dim, rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                assertion = assertion .and. abs(rweisum - rweisum_ref) <= rtol * abs(rweisum_ref)
                call report(__LINE__, SK_"The `weisum` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean_ref = getMeanD2(sample, dim)
                call setMean(mean, sample, dim)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

            ! test sample D2 ALL interface.

            block

                dim = 0
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(mean, 1_IK)
                if (getUnifRand()) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKC, 9._TKC, nsam * ndim)
                iweight = getUnifRand(1_IK, 9_IK, nsam * ndim)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                mean_ref = getMeanD2(sample, weight = real(iweight, TKC))
                call setMean(mean(1), sample, iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                assertion = assertion .and. iweisum == iweisum_ref
                call report(__LINE__, SK_"The `weisum` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean_ref = getMeanD2(sample, weight = rweight)
                call setMean(mean(1), sample, rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                assertion = assertion .and. abs(rweisum - rweisum_ref) <= rtol * abs(rweisum_ref)
                call report(__LINE__, SK_"The `weisum` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean_ref = getMeanD2(sample)
                call setMean(mean(1), sample)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

            ! test sample D1 DIM interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(mean, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                mean_ref = getMeanD1(sample(:,1), real(iweight, TKC))
                call setMean(mean(1), sample(:,1), dim, iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                assertion = assertion .and. iweisum == iweisum_ref
                call report(__LINE__, SK_"The `weisum` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean_ref = getMeanD1(sample(:,1), weight = rweight)
                call setMean(mean(1), sample(:,1), dim, rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                assertion = assertion .and. abs(rweisum - rweisum_ref) <= rtol * abs(rweisum_ref)
                call report(__LINE__, SK_"The `weisum` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean_ref = getMeanD1(sample(:,1))
                call setMean(mean(1), sample(:,1), dim)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

            ! test sample D1 ALL interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(mean, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKC, 9._TKC, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                mean_ref = getMeanD1(sample(:,1), real(iweight, TKC))
                call setMean(mean(1), sample(:,1), iweight, iweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for integer-weighted sample.")

                assertion = assertion .and. iweisum == iweisum_ref
                call report(__LINE__, SK_"The `weisum` must be computed correctly for integer-weighted sample.")

                ! real weighted

                mean_ref = getMeanD1(sample(:,1), weight = rweight)
                call setMean(mean(1), sample(:,1), rweight, rweisum)

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for real-weighted sample.")

                assertion = assertion .and. abs(rweisum - rweisum_ref) <= rtol * abs(rweisum_ref)
                call report(__LINE__, SK_"The `weisum` must be computed correctly for real-weighted sample.")

                ! unweighted

                mean_ref = getMeanD1(sample(:,1))
                call setMean(mean(1), sample(:,1))

                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The `mean` must be computed correctly for unweighted sample.")

            end block

        end do

    contains

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
                call test%disp%show("iweisum")
                call test%disp%show( iweisum )
                call test%disp%show("iweisum_ref")
                call test%disp%show( iweisum_ref )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("rweisum")
                call test%disp%show( rweisum )
                call test%disp%show("rweisum_ref")
                call test%disp%show( rweisum_ref )
                call test%disp%show("mean_ref")
                call test%disp%show( mean_ref )
                call test%disp%show("mean")
                call test%disp%show( mean )
                call test%disp%show("diff")
                call test%disp%show( diff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        pure function getMeanD1(sample, weight) result(mean)
            real(TKC), intent(in), optional :: weight(:)
            TYPE_OF_SAMPLE, intent(in) :: sample(:)
            TYPE_OF_SAMPLE :: mean(1)
            if (present(weight)) then
                mean = sum(sample * weight) / sum(weight)
            else
                mean = sum(sample) / size(sample, 1, IK)
            end if
        end function

        pure function getMeanD2(sample, dim, weight) result(mean)
            real(TKC), intent(in), optional :: weight(:)
            TYPE_OF_SAMPLE, intent(in) :: sample(:,:)
            integer(IK), intent(in), optional :: dim
            TYPE_OF_SAMPLE, allocatable :: mean(:)
            integer(IK) :: idim
            if (present(dim)) then
                allocate(mean(size(sample, 3 - dim, IK)))
                do idim = 1, size(sample, 3 - dim, IK)
                    if (dim == 2) then
                        mean(idim:idim) = getMeanD1(sample(idim,:), weight)
                    else
                        mean(idim:idim) = getMeanD1(sample(:,idim), weight)
                    end if
                end do
            else
                if (present(weight)) then
                    mean = [sum(reshape(sample, [product(shape(sample))]) * weight) / sum(weight)]
                else
                    mean = [sum(sample) / size(sample, kind = IK)]
                end if
            end if
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getMeanMerged_ENABLED || setMeanMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: fracA
        integer(IK) :: itry, ndim, dim
        integer(IK) :: nsam, nsamA, nsamB
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), meanA(:), meanB(:), mean(:), mean_ref(:), diff(:)
        assertion = .true._LK
        dim = 2_IK

        do itry = 1, 50

            nsamA = getUnifRand(1_IK, 5_IK)
            nsamB = getUnifRand(1_IK, 5_IK)
            nsam = nsamA + nsamB
            fracA = real(nsamA, TKC) / real(nsam, TKC)

            ! test D2 interface.

            block

                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(mean, ndim)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                mean_ref = getMean(sample, dim)
                meanA = getMean(sample(:,1:nsamA), dim)
                meanB = getMean(sample(:,nsamA+1:), dim)

#if             getMeanMerged_ENABLED
                mean = getMeanMerged(meanB, meanA, fracA)
                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The new `meanMerged` must be computed correctly.")
#elif           setMeanMerged_ENABLED
                ! new
                call setMeanMerged(mean, meanB, meanA, fracA)
                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The new `meanMerged` must be computed correctly.")
                ! old
                mean = meanB
                call setMeanMerged(mean, meanA, fracA)
                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The in-place `meanMerged` must be computed correctly.")
#endif

            end block

            ! test D1 interface.

            block

                ndim = 1_IK
                call setResized(mean, ndim)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                mean_ref = getMean(sample, dim)
                meanA = getMean(sample(:,1:nsamA), dim)
                meanB = getMean(sample(:,nsamA+1:), dim)

#if             getMeanMerged_ENABLED
                mean(1) = getMeanMerged(meanB(1), meanA(1), fracA)
                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The new `meanMerged` must be computed correctly.")
#elif           setMeanMerged_ENABLED
                ! new
                call setMeanMerged(mean(1), meanB(1), meanA(1), fracA)
                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The new `meanMerged` must be computed correctly.")
                ! old
                mean = meanB
                call setMeanMerged(mean(1), meanA(1), fracA)
                diff = abs(mean - mean_ref)
                assertion = assertion .and. all(diff < ctol)
                call report(__LINE__, SK_"The in-place `meanMerged` must be computed correctly.")
#endif
            end block

        end do

    contains

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
                call test%disp%show("diff")
                call test%disp%show( diff )
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