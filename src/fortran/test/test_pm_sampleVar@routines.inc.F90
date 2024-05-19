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
!>  This file contains procedure implementations of tests of [pm_sampleVar](@ref pm_sampleVar).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG), parameter :: rtol = epsilon(1._TKG) * 100
        ! Define the sample type.
#if     CK_ENABLED
#define TYPE_OF_SAMPLE complex(TKG)
        complex(TKG), parameter :: ONE = (1._TKG, 1._TKG), TWO = 2 * (1._TKG, 1._TKG), ctol = (rtol, rtol)
#elif   RK_ENABLED
#define TYPE_OF_SAMPLE real(TKG)
        real(TKG), parameter :: ONE = 1._TKG, TWO = 2._TKG, ctol = rtol
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%
#if     getVarCorrection_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: itry
        real(TKG) :: weisum, weisqs, correction, correction_ref
        do itry = 1, 10

            ! Test reliability weight.

            weisum = getUnifRand(2._TKG, 3._TKG)
            weisqs = getUnifRand(3._TKG, 4._TKG)
            correction = getVarCorrection(weisum, weisqs)
            correction_ref = weisum**2 / (weisum**2 - weisqs)
            assertion = logical(correction == correction_ref, LK)
            call report(__LINE__, SK_"The Bessel Reliability correction must be correctly computed.")

            weisum = getUnifRand(2._TKG, 3._TKG)
            correction = getVarCorrection(weisum)
            correction_ref = weisum / (weisum - 1._TKG)
            assertion = logical(correction == correction_ref, LK)
            call report(__LINE__, SK_"The Bessel Frequency correction must be correctly computed.")

        end do

    contains

        subroutine report(line, msg)
            integer, intent(in) :: line
            character(*, SK), intent(in) :: msg
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[weisum, weisqs, correction, correction_ref]")
                call test%disp%show( [weisum, weisqs, correction, correction_ref] )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%
#elif   getVar_ENABLED
        !%%%%%%%%%%%%%

        real(TKG) :: rweisum
        integer(IK) :: iweisum
        integer(IK) :: itry, nsam, ndim, dim
        real(TKG), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:)
        real(TKG), allocatable :: var(:), var_ref(:), vdiff(:)
        assertion = .true._LK

        do itry = 1, 50

            dim = 0
            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample D2 DIM interface.

            block

                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(var_ref, ndim)
                call setResized(var, ndim)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                var = getVar(sample, dim, iweight)
                mean = getMean(sample, dim, iweight)
                call setVar(var_ref, mean, sample, dim, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                var = getVar(sample, dim, iweight, fweight_type())
                var_ref = var_ref * getVarCorrection(real(iweisum, TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! real weighted

                var = getVar(sample, dim, rweight)
                mean = getMean(sample, dim, rweight)
                call setVar(var_ref, mean, sample, dim, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                var = getVar(sample, dim, rweight, rweight_type())
                var_ref = var_ref * getVarCorrection(rweisum, sum(rweight**2))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! unweighted

                var = getVar(sample, dim)
                mean = getMean(sample, dim)
                call setVar(var_ref, mean, sample, dim)
                call setAssertedVar(__LINE__, SK_"unweighted")

                var = getVar(sample, dim, fweight_type())
                var_ref = var_ref * getVarCorrection(real(size(sample, dim), TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

            end block

            ! test sample D2 ALL interface.

            block

                dim = 0
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(var_ref, 1_IK)
                call setResized(mean, 1_IK)
                call setResized(var, 1_IK)
                if (getUnifRand()) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKG, 9._TKG, nsam * ndim)
                iweight = getUnifRand(1_IK, 9_IK, nsam * ndim)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                var(1) = getVar(sample, iweight)
                mean(1) = getMean(sample, iweight)
                call setVar(var_ref(1), mean(1), sample, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                var(1) = getVar(sample, iweight, fweight_type())
                var_ref = var_ref * getVarCorrection(real(iweisum, TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! real weighted

                var(1) = getVar(sample, rweight)
                mean(1) = getMean(sample, rweight)
                call setVar(var_ref(1), mean(1), sample, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                var(1) = getVar(sample, rweight, rweight_type())
                var_ref = var_ref * getVarCorrection(rweisum, sum(rweight**2))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! unweighted

                var(1) = getVar(sample)
                mean(1) = getMean(sample)
                call setVar(var_ref(1), mean(1), sample)
                call setAssertedVar(__LINE__, SK_"unweighted")

                var(1) = getVar(sample)
                var_ref = var_ref * getVarCorrection(real(size(sample), TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

            end block

            ! test sample D1 DIM interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(var, 1_IK)
                call setResized(var_ref, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                var(1) = getVar(sample(:,1), dim, iweight)
                mean(1) = getMean(sample(:,1), dim, iweight)
                call setVar(var_ref(1), mean(1), sample(:,1), dim, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                var(1) = getVar(sample(:,1), dim, iweight, fweight_type())
                var_ref = var_ref * getVarCorrection(real(iweisum, TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! real weighted

                var(1) = getVar(sample(:,1), dim, rweight)
                mean(1) = getMean(sample(:,1), dim, rweight)
                call setVar(var_ref(1), mean(1), sample(:,1), dim, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                var(1) = getVar(sample(:,1), dim, rweight, rweight_type())
                var_ref = var_ref * getVarCorrection(rweisum, sum(rweight**2))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! unweighted

                var(1) = getVar(sample(:,1), dim)
                mean(1) = getMean(sample(:,1), dim)
                call setVar(var_ref(1), mean(1), sample(:,1), dim)
                call setAssertedVar(__LINE__, SK_"unweighted")

                var(1) = getVar(sample(:,1), dim)
                var_ref = var_ref * getVarCorrection(real(size(sample), TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

            end block

            ! test sample D1 ALL interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(var, 1_IK)
                call setResized(var_ref, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                var(1) = getVar(sample(:,1), iweight)
                mean(1) = getMean(sample(:,1), iweight)
                call setVar(var_ref(1), mean(1), sample(:,1), iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                var(1) = getVar(sample(:,1), iweight, fweight_type())
                var_ref = var_ref * getVarCorrection(real(iweisum, TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! real weighted

                var(1) = getVar(sample(:,1), rweight)
                mean(1) = getMean(sample(:,1), rweight)
                call setVar(var_ref(1), mean(1), sample(:,1), rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                var(1) = getVar(sample(:,1), rweight, rweight_type())
                var_ref = var_ref * getVarCorrection(rweisum, sum(rweight**2))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

                ! unweighted

                var(1) = getVar(sample(:,1))
                mean(1) = getMean(sample(:,1))
                call setVar(var_ref(1), mean(1), sample(:,1))
                call setAssertedVar(__LINE__, SK_"unweighted")

                var(1) = getVar(sample(:,1))
                var_ref = var_ref * getVarCorrection(real(size(sample), TKG))
                call report(__LINE__, SK_"The unbiased `var` must be computed correctly.")

            end block

        end do

    contains

        subroutine setAssertedVar(line, this)
            character(*, SK), intent(in) :: this
            integer, intent(in) :: line
            vdiff = abs(var - var_ref)
            assertion = assertion .and. all(vdiff < rtol)
            call report(line, SK_"The `var` must be computed correctly for "//this//SK_" sample.")
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
                call test%disp%show("var_ref")
                call test%disp%show( var_ref )
                call test%disp%show("var")
                call test%disp%show( var )
                call test%disp%show("vdiff")
                call test%disp%show( vdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%
#elif   setVar_ENABLED
        !%%%%%%%%%%%%%

        real(TKG) :: rweisum
        integer(IK) :: iweisum
        real(TKG), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam, ndim, dim
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:)
        real(TKG), allocatable :: var(:), var_ref(:), vdiff(:)
        assertion = .true._LK

        do itry = 1, 50

            dim = 0
            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample D2 DIM interface.

            block

                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(var, ndim)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, dim, iweight)
                var_ref = getVarD2(mean, sample, dim, real(iweight, TKG))

                call setVar(var, mean, sample, dim, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                call setVar(var, getShifted(sample, dim, -mean), dim, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted shifted")

                ! real weighted

                mean = getMean(sample, dim, rweight)
                var_ref = getVarD2(mean, sample, dim, rweight)

                call setVar(var, mean, sample, dim, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                call setVar(var, getShifted(sample, dim, -mean), dim, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted shifted")

                ! unweighted

                mean = getMean(sample, dim)
                var_ref = getVarD2(mean, sample, dim)

                call setVar(var, mean, sample, dim)
                call setAssertedVar(__LINE__, SK_"unweighted")

                call setVar(var, getShifted(sample, dim, -mean), dim)
                call setAssertedVar(__LINE__, SK_"unweighted shifted")

            end block

            ! test sample D2 ALL interface.

            block

                dim = 0
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(var_ref, 1_IK)
                call setResized(mean, 1_IK)
                call setResized(var, 1_IK)
                if (getUnifRand()) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                end if
                rweight = getUnifRand(1._TKG, 9._TKG, nsam * ndim)
                iweight = getUnifRand(1_IK, 9_IK, nsam * ndim)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample, iweight)
                var_ref = getVarD2(mean, sample, weight = real(iweight, TKG))

                call setVar(var(1), mean(1), sample, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                call setVar(var(1), getShifted(sample, -mean(1)), iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted shifted")

                ! real weighted

                mean(1) = getMean(sample, rweight)
                var_ref = getVarD2(mean, sample, weight = rweight)

                call setVar(var(1), mean(1), sample, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                call setVar(var(1), getShifted(sample, -mean(1)), rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted shifted")

                ! unweighted

                mean(1) = getMean(sample)
                var_ref = getVarD2(mean, sample)

                call setVar(var(1), mean(1), sample)
                call setAssertedVar(__LINE__, SK_"unweighted")

                call setVar(var(1), getShifted(sample, -mean(1)))
                call setAssertedVar(__LINE__, SK_"unweighted shifted")

            end block

            ! test sample D1 DIM interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(var, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample(:,1), dim, iweight)
                var_ref = getVarD1(mean(1), sample(:,1), real(iweight, TKG))

                call setVar(var(1), mean(1), sample(:,1), dim, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                call setVar(var(1), getShifted(sample(:,1), dim, -mean(1)), dim, iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted shifted")

                ! real weighted

                mean = getMean(sample(:,1), rweight)
                var_ref = getVarD1(mean(1), sample(:,1), weight = rweight)

                call setVar(var(1), mean(1), sample(:,1), dim, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                call setVar(var(1), getShifted(sample(:,1), dim, -mean(1)), dim, rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted shifted")

                ! unweighted

                mean = getMean(sample(:,1))
                var_ref = getVarD1(mean(1), sample(:,1))

                call setVar(var(1), mean(1), sample(:,1), dim)
                call setAssertedVar(__LINE__, SK_"unweighted")

                call setVar(var(1), getShifted(sample(:,1), dim, -mean(1)), dim)
                call setAssertedVar(__LINE__, SK_"unweighted shifted")

            end block

            ! test sample D1 ALL interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(var, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum = sum(iweight)
                rweisum = sum(rweight)

                ! integer weighted

                mean = getMean(sample(:,1), iweight)
                var_ref = getVarD1(mean(1), sample(:,1), real(iweight, TKG))

                call setVar(var(1), mean(1), sample(:,1), iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                call setVar(var(1), getShifted(sample(:,1), -mean(1)), iweight, iweisum)
                call setAssertedVar(__LINE__, SK_"integer-weighted shifted")

                ! real weighted

                mean = getMean(sample(:,1), rweight)
                var_ref = getVarD1(mean(1), sample(:,1), weight = rweight)

                call setVar(var(1), mean(1), sample(:,1), rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted")

                call setVar(var(1), getShifted(sample(:,1), -mean(1)), rweight, rweisum)
                call setAssertedVar(__LINE__, SK_"real-weighted shifted")

                ! unweighted

                mean = getMean(sample(:,1))
                var_ref = getVarD1(mean(1), sample(:,1))
                call setVar(var(1), mean(1), sample(:,1))
                call setAssertedVar(__LINE__, SK_"unweighted")

                call setVar(var(1), getShifted(sample(:,1), -mean(1)))
                call setAssertedVar(__LINE__, SK_"unweighted shifted")

            end block

        end do

    contains

        subroutine setAssertedVar(line, this)
            character(*, SK), intent(in) :: this
            integer, intent(in) :: line
            vdiff = abs(var - var_ref)
            assertion = assertion .and. all(vdiff < rtol)
            call report(line, SK_"The `var` must be computed correctly for "//this//SK_" sample.")
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
                call test%disp%show("var_ref")
                call test%disp%show( var_ref )
                call test%disp%show("var")
                call test%disp%show( var )
                call test%disp%show("vdiff")
                call test%disp%show( vdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        PURE function getVarD1(mean, sample, weight) result(var)
            real(TKG), intent(in), optional :: weight(:)
            TYPE_OF_SAMPLE, intent(in) :: sample(:), mean
            TYPE_OF_SAMPLE :: sampleShifted(size(sample, 1, IK))
            real(TKG) :: var
            sampleShifted = getShifted(sample, -mean)
            if (present(weight)) then
                var = real(dot_product(sampleShifted, sampleShifted * weight), TKG) / sum(weight)
            else
                var = real(dot_product(sampleShifted, sampleShifted), TKG) / size(sample, 1, IK)
            end if
        end function

        PURE function getVarD2(mean, sample, dim, weight) result(var)
            real(TKG), intent(in), optional :: weight(:)
            TYPE_OF_SAMPLE, intent(in) :: sample(:,:), mean(:)
            integer(IK), intent(in), optional :: dim
            real(TKG), allocatable :: var(:)
            integer(IK) :: idim
            if (present(dim)) then
                allocate(var(size(sample, 3 - dim, IK)))
                do idim = 1, size(sample, 3 - dim, IK)
                    if (dim == 2) then
                        var(idim) = getVarD1(mean(idim), sample(idim,:), weight)
                    else
                        var(idim) = getVarD1(mean(idim), sample(:,idim), weight)
                    end if
                end do
            else
                block
                    TYPE_OF_SAMPLE :: sampleShifted(size(sample, kind = IK))
                    sampleShifted = reshape(sample, shape(sampleShifted)) - mean(1)
                    if (present(weight)) then
                        var = [real(dot_product(sampleShifted, sampleShifted * weight), TKG) / sum(weight)]
                    else
                        var = [real(dot_product(sampleShifted, sampleShifted), TKG) / size(sampleShifted, 1, IK)]
                    end if
                end block
            end if
        end function

        !%%%%%%%%%%%%%%%%%
#elif   setVarMean_ENABLED
        !%%%%%%%%%%%%%%%%%

        real(TKG) :: rweisum, rweisum_ref
        integer(IK) :: iweisum, iweisum_ref
        real(TKG), allocatable :: rweight(:)
        integer(IK), allocatable :: iweight(:)
        integer(IK) :: itry, nsam, ndim, dim
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:), mean_ref(:), meang(:), mdiff(:)
        real(TKG), allocatable :: var(:), var_ref(:), vdiff(:)
        assertion = .true._LK

        do itry = 1, 50

            dim = 0
            nsam = getUnifRand(2_IK, 7_IK)

            ! test sample D2 DIM interface.

            block

                dim = getChoice([1, 2])
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(mean_ref, ndim)
                call setResized(var_ref, ndim)
                call setResized(mean, ndim)
                call setResized(var, ndim)
                if (dim == 2) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                    meang = sample(:,1)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                    meang = sample(1,:)
                end if
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                var_ref = getVar(sample, dim, iweight)
                mean_ref = getMean(sample, dim, iweight)
                call setVarMean(var, mean, sample, dim, iweight, iweisum, meang)
                call setAssertedSum(__LINE__, SK_"integer-weighted", real(iweisum, TKG), real(iweisum_ref, TKG))
                call setAssertedAvg(__LINE__, SK_"integer-weighted")
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                ! real weighted

                var_ref = getVar(sample, dim, rweight)
                mean_ref = getMean(sample, dim, rweight)
                call setVarMean(var, mean, sample, dim, rweight, rweisum, meang)
                call setAssertedSum(__LINE__, SK_"real-weighted", rweisum, rweisum_ref)
                call setAssertedAvg(__LINE__, SK_"real-weighted")
                call setAssertedVar(__LINE__, SK_"real-weighted")

                ! unweighted

                var_ref = getVar(sample, dim)
                mean_ref = getMean(sample, dim)
                call setVarMean(var, mean, sample, dim, meang)
                call setAssertedAvg(__LINE__, SK_"unweighted")
                call setAssertedVar(__LINE__, SK_"unweighted")

            end block

            ! test sample D2 ALL interface.

            block

                dim = 0
                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(mean_ref, 1_IK)
                call setResized(var_ref, 1_IK)
                call setResized(mean, 1_IK)
                call setResized(var, 1_IK)
                if (getUnifRand()) then
                    sample = getUnifRand(ONE, TWO, ndim, nsam)
                    meang = sample(:,1)
                else
                    sample = getUnifRand(ONE, TWO, nsam, ndim)
                    meang = sample(1,:)
                end if
                rweight = getUnifRand(1._TKG, 9._TKG, nsam * ndim)
                iweight = getUnifRand(1_IK, 9_IK, nsam * ndim)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                var_ref = getVar(sample, iweight)
                mean_ref = getMean(sample, iweight)
                call setVarMean(var(1), mean(1), sample, iweight, iweisum, meang(1))
                call setAssertedSum(__LINE__, SK_"integer-weighted", real(iweisum, TKG), real(iweisum_ref, TKG))
                call setAssertedAvg(__LINE__, SK_"integer-weighted")
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                ! real weighted

                var_ref = getVar(sample, rweight)
                mean_ref = getMean(sample, rweight)
                call setVarMean(var(1), mean(1), sample, rweight, rweisum, meang(1))
                call setAssertedSum(__LINE__, SK_"real-weighted", rweisum, rweisum_ref)
                call setAssertedAvg(__LINE__, SK_"real-weighted")
                call setAssertedVar(__LINE__, SK_"real-weighted")

                ! unweighted

                var_ref = getVar(sample)
                mean_ref = getMean(sample)
                call setVarMean(var(1), mean(1), sample, meang(1))
                call setAssertedAvg(__LINE__, SK_"unweighted")
                call setAssertedVar(__LINE__, SK_"unweighted")

            end block

            ! test sample D1 DIM interface.

            block

                dim = 1_IK
                ndim = 1_IK
                call setResized(var, 1_IK)
                call setResized(mean, 1_IK)
                call setResized(var_ref, 1_IK)
                call setResized(mean_ref, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                meang = sample(1,:)
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                var_ref = getVar(sample(:,1), dim, iweight)
                mean_ref = getMean(sample(:,1), dim, iweight)
                call setVarMean(var(1), mean(1), sample(:,1), dim, iweight, iweisum, meang(1))
                call setAssertedSum(__LINE__, SK_"integer-weighted", real(iweisum, TKG), real(iweisum_ref, TKG))
                call setAssertedAvg(__LINE__, SK_"integer-weighted")
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                ! real weighted

                var_ref = getVar(sample(:,1), dim, rweight)
                mean_ref = getMean(sample(:,1), dim, rweight)
                call setVarMean(var(1), mean(1), sample(:,1), dim, rweight, rweisum, meang(1))
                call setAssertedSum(__LINE__, SK_"real-weighted", rweisum, rweisum_ref)
                call setAssertedAvg(__LINE__, SK_"real-weighted")
                call setAssertedVar(__LINE__, SK_"real-weighted")

                ! unweighted

                var_ref = getVar(sample(:,1), dim)
                mean_ref = getMean(sample(:,1), dim)
                call setVarMean(var(1), mean(1), sample(:,1), dim, meang(1))
                call setAssertedAvg(__LINE__, SK_"unweighted")
                call setAssertedVar(__LINE__, SK_"unweighted")

            end block

            ! test sample D1 ALL interface.

            block

                dim = 0_IK
                ndim = 1_IK
                call setResized(var, 1_IK)
                call setResized(mean, 1_IK)
                call setResized(var_ref, 1_IK)
                call setResized(mean_ref, 1_IK)
                sample = getUnifRand(ONE, TWO, nsam, ndim)
                meang = sample(1,:)
                rweight = getUnifRand(1._TKG, 9._TKG, nsam)
                iweight = getUnifRand(1_IK, 9_IK, nsam)
                iweisum_ref = sum(iweight)
                rweisum_ref = sum(rweight)
                iweisum = 0
                rweisum = 0

                ! integer weighted

                var_ref = getVar(sample(:,1), iweight)
                mean_ref = getMean(sample(:,1), iweight)
                call setVarMean(var(1), mean(1), sample(:,1), iweight, iweisum, meang(1))
                call setAssertedSum(__LINE__, SK_"integer-weighted", real(iweisum, TKG), real(iweisum_ref, TKG))
                call setAssertedAvg(__LINE__, SK_"integer-weighted")
                call setAssertedVar(__LINE__, SK_"integer-weighted")

                ! real weighted

                var_ref = getVar(sample(:,1), rweight)
                mean_ref = getMean(sample(:,1), rweight)
                call setVarMean(var(1), mean(1), sample(:,1), rweight, rweisum, meang(1))
                call setAssertedSum(__LINE__, SK_"real-weighted", rweisum, rweisum_ref)
                call setAssertedAvg(__LINE__, SK_"real-weighted")
                call setAssertedVar(__LINE__, SK_"real-weighted")

                ! unweighted

                var_ref = getVar(sample(:,1))
                mean_ref = getMean(sample(:,1))
                call setVarMean(var(1), mean(1), sample(:,1), meang(1))
                call setAssertedAvg(__LINE__, SK_"unweighted")
                call setAssertedVar(__LINE__, SK_"unweighted")

            end block

        end do

    contains

        subroutine setAssertedSum(line, this, weisum, weisum_ref)
            real(TKG), intent(in) :: weisum, weisum_ref
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

        subroutine setAssertedVar(line, this)
            character(*, SK), intent(in) :: this
            integer, intent(in) :: line
            vdiff = abs(var - var_ref)
            assertion = assertion .and. all(vdiff < rtol)
            call report(line, SK_"The `var` must be computed correctly for "//this//SK_" sample.")
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
                call test%disp%show("meang")
                call test%disp%show( meang )
                call test%disp%show("iweight")
                call test%disp%show( iweight )
                call test%disp%show("rweight")
                call test%disp%show( rweight )
                call test%disp%show("iweisum")
                call test%disp%show( iweisum )
                call test%disp%show("iweisum_ref")
                call test%disp%show( iweisum_ref )
                call test%disp%show("rweisum")
                call test%disp%show( rweisum )
                call test%disp%show("rweisum_ref")
                call test%disp%show( rweisum_ref )
                call test%disp%show("var_ref")
                call test%disp%show( var_ref )
                call test%disp%show("var")
                call test%disp%show( var )
                call test%disp%show("vdiff")
                call test%disp%show( vdiff )
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
#elif   getVarMerged_ENABLED || setVarMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: fracA
        integer(IK) :: itry, ndim, dim
        integer(IK) :: nsam, nsamA, nsamB
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), meanDiff(:)
        real(TKG), allocatable :: varA(:), varB(:), var(:), var_ref(:), vdiff(:)
        assertion = .true._LK
        dim = 2_IK

        do itry = 1, 50

            nsamA = getUnifRand(2_IK, 5_IK)
            nsamB = getUnifRand(2_IK, 5_IK)
            nsam = nsamA + nsamB
            fracA = real(nsamA, TKG) / real(nsam, TKG)

            ! test D2 interface.

            block

                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(var, ndim)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                var_ref = getVar(sample, dim)
                varA = getVar(sample(:,1:nsamA), dim)
                varB = getVar(sample(:,nsamA+1:), dim)
                meanDiff = getMean(sample(:,1:nsamA), dim) - getMean(sample(:,nsamA+1:), dim)

#if             getVarMerged_ENABLED
                var = getVarMerged(varB, varA, meanDiff, fracA)
                vdiff = abs(var - var_ref)
                assertion = assertion .and. all(vdiff < rtol)
                call report(__LINE__, SK_"The new `varMerged` must be computed correctly.")
#elif           setVarMerged_ENABLED
                ! new
                call setVarMerged(var, varB, varA, meanDiff, fracA)
                vdiff = abs(var - var_ref)
                assertion = assertion .and. all(vdiff < rtol)
                call report(__LINE__, SK_"The new `varMerged` must be computed correctly.")
                ! old
                var = varB
                call setVarMerged(var, varA, meanDiff, fracA)
                vdiff = abs(var - var_ref)
                assertion = assertion .and. all(vdiff < rtol)
                call report(__LINE__, SK_"The in-place `varMerged` must be computed correctly.")
#endif

            end block

            ! test D1 interface.

            block

                ndim = 1_IK
                call setResized(var, ndim)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                var_ref = getVar(sample, dim)
                varA = getVar(sample(:,1:nsamA), dim)
                varB = getVar(sample(:,nsamA+1:), dim)
                meanDiff = getMean(sample(:,1:nsamA), dim) - getMean(sample(:,nsamA+1:), dim)

#if             getVarMerged_ENABLED
                var(1) = getVarMerged(varB(1), varA(1), meanDiff(1), fracA)
                vdiff = abs(var - var_ref)
                assertion = assertion .and. all(vdiff < rtol)
                call report(__LINE__, SK_"The new `varMerged` must be computed correctly.")
#elif           setVarMerged_ENABLED
                ! new
                call setVarMerged(var(1), varB(1), varA(1), meanDiff(1), fracA)
                vdiff = abs(var - var_ref)
                assertion = assertion .and. all(vdiff < rtol)
                call report(__LINE__, SK_"The new `varMerged` must be computed correctly.")
                ! old
                var = varB
                call setVarMerged(var(1), varA(1), meanDiff(1), fracA)
                vdiff = abs(var - var_ref)
                assertion = assertion .and. all(vdiff < rtol)
                call report(__LINE__, SK_"The in-place `varMerged` must be computed correctly.")
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
                call test%disp%show("meanDiff")
                call test%disp%show( meanDiff )
                call test%disp%show("varA")
                call test%disp%show( varA )
                call test%disp%show("varB")
                call test%disp%show( varB )
                call test%disp%show("var_ref")
                call test%disp%show( var_ref )
                call test%disp%show("var")
                call test%disp%show( var )
                call test%disp%show("vdiff")
                call test%disp%show( vdiff )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, msg, int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%
#elif   setVarMeanMerged_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        real(TKG) :: fracA
        integer(IK) :: itry, ndim, dim
        integer(IK) :: nsam, nsamA, nsamB
        TYPE_OF_SAMPLE, allocatable :: sample(:,:), mean(:), mean_ref(:), meanA(:), meanB(:), mdiff(:)
        real(TKG), allocatable :: varA(:), varB(:), var(:), var_ref(:), vdiff(:)
        assertion = .true._LK
        dim = 2_IK

        do itry = 1, 50

            nsamA = getUnifRand(2_IK, 5_IK)
            nsamB = getUnifRand(2_IK, 5_IK)
            nsam = nsamA + nsamB
            fracA = real(nsamA, TKG) / real(nsam, TKG)

            ! test D2 interface.

            block

                ndim = getUnifRand(1_IK, 5_IK)
                call setResized(var, ndim)
                call setResized(mean, ndim)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                var_ref = getVar(sample, dim)
                mean_ref = getMean(sample, dim)
                varA = getVar(sample(:,1:nsamA), dim)
                varB = getVar(sample(:,nsamA+1:), dim)
                meanA = getMean(sample(:,1:nsamA), dim)
                meanB = getMean(sample(:,nsamA+1:), dim)

                ! new
                call setVarMeanMerged(var, mean, varB, meanB, varA, meanA, fracA)
                call setAssertedAvg(__LINE__, SK_"new")
                call setAssertedVar(__LINE__, SK_"new")

                ! old
                var = varB
                mean = meanB
                call setVarMeanMerged(var, mean, varA, meanA, fracA)
                call setAssertedAvg(__LINE__, SK_"in-place")
                call setAssertedVar(__LINE__, SK_"in-place")

            end block

            ! test D1 interface.

            block

                ndim = 1_IK
                call setResized(var, ndim)
                call setResized(mean, ndim)
                sample = getUnifRand(ONE, TWO, ndim, nsam)

                var_ref = getVar(sample, dim)
                mean_ref = getMean(sample, dim)
                varA = getVar(sample(:,1:nsamA), dim)
                varB = getVar(sample(:,nsamA+1:), dim)
                meanA = getMean(sample(:,1:nsamA), dim)
                meanB = getMean(sample(:,nsamA+1:), dim)

                ! new
                call setVarMeanMerged(var(1), mean(1), varB(1), meanB(1), varA(1), meanA(1), fracA)
                call setAssertedAvg(__LINE__, SK_"new")
                call setAssertedVar(__LINE__, SK_"new")

                ! old
                var = varB
                mean = meanB
                call setVarMeanMerged(var(1), mean(1), varA(1), meanA(1), fracA)
                call setAssertedAvg(__LINE__, SK_"in-place")
                call setAssertedVar(__LINE__, SK_"in-place")

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

        subroutine setAssertedVar(line, specific)
            character(*, SK), intent(in) :: specific
            integer, intent(in) :: line
            vdiff = abs(var - var_ref)
            assertion = assertion .and. all(vdiff < rtol)
            call report(line, SK_"The "//specific//SK_" `varMerged` must be computed correctly.")
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
                call test%disp%show("varA")
                call test%disp%show( varA )
                call test%disp%show("varB")
                call test%disp%show( varB )
                call test%disp%show("var_ref")
                call test%disp%show( var_ref )
                call test%disp%show("var")
                call test%disp%show( var )
                call test%disp%show("vdiff")
                call test%disp%show( vdiff )
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