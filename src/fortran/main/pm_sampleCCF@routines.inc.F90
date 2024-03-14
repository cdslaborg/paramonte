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
!>  This file contains the implementation details of the 2D routines under the generic interface [pm_sampleCCF](@ref pm_sampleCCF).
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Saturday 4:40 PM, August 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Set the type and kind and the conjugation rule.
#if     CK_ENABLED
#define GET_ABSQ(X)(real(X)**2 + aimag(X)**2)
#define TYPE_OF_SEQ complex(TKC)
#define GET_RE(X)X%re
#elif   RK_ENABLED
#define TYPE_OF_SEQ real(TKC)
#define GET_ABSQ(X)X**2
#define GET_RE(X)X
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getACF_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(TKC) :: normfac
        TYPE_OF_SEQ :: meanf
        class(*), allocatable :: norm_def
        TYPE_OF_SEQ, parameter :: ZERO = 0._TKC
        TYPE_OF_SEQ, allocatable :: ff(:), coef(:)
        integer(IK), allocatable :: factor(:)
        logical(LK) :: unscaled, inf
        integer(IK) :: lagmax
        integer(IK) :: lagmin
        integer(IK) :: lenacf

        lagmax = size(f, 1, IK) - 1_IK
        lagmin = 1_IK - size(f, 1, IK)
        if (present(lag)) then
            if (size(lag, 1, IK) == 0_IK) then
                allocate(acf(0))
                return
            end if
            CHECK_ASSERTION(__LINE__, isAscending(lag), SK_"@setACF(): The condition `isAscending(lag)` must hold. lag = "//getStr(lag))
            lagmax = max(lagmax, lag(size(lag, 1, IK)))
            lagmin = min(lagmin, lag(1))
        end if
        lenacf = lagmax - lagmin + 1_IK

        ! Allocate and pad the arrays.

        allocate(ff(lenacf), coef(lenacf), acf(lenacf))

        ! Shift the padded arrays if needed.

        if (present(norm)) then
            unscaled = logical(same_type_as(norm, zscore_type()) .or. same_type_as(norm, stdscale_type()), LK)
            norm_def = norm
        else
            unscaled = .true._LK
            norm_def = zscore_type() ! Intel cannot digest `same_type_as(norm_def, zscore_type())` if norm_def is set to constant.
        end if

        if (same_type_as(norm_def, meanshift_type()) .or. same_type_as(norm_def, zscore_type())) then
            meanf = getMean(f)
            ff(1 : size(f, 1, IK)) = f - meanf
        elseif (same_type_as(norm_def, stdscale_type()) .or. same_type_as(norm_def, nothing_type())) then
            ff(1 : size(f, 1, IK)) = f
        else
            error stop MODULE_NAME//SK_"@getACF(): Unrecognized `norm`." ! LCOV_EXCL_LINE
        end if
        ff(size(f, 1, IK) + 1 : lenacf) = ZERO

        ! Compute the ACF.

        factor = getFactorFFT(ff, coef)
        call setACF(factor, coef, ff, acf, inf)

        ! Normalize and shift.

        normfac = 1._TKC / lenacf ! default scale.
        if (inf) then
            if (unscaled) normfac = 1._TKC / GET_RE(ff(1))
            if (present(lag)) then
                acf = cshift(ff * normfac, lagmin)
            else
                acf = ff(1 : 1 + lagmax) * normfac
            end if
        else
            if (unscaled) normfac = 1._TKC / GET_RE(acf(1))
            if (present(lag)) then
                acf = cshift(acf * normfac, lagmin)
            else
                acf = acf(1 : 1 + lagmax) * normfac
            end if
        end if
        if (present(lag)) acf = acf(lag - lagmin + 1)
        deallocate(ff, coef) ! essential for gfortran on heap.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setACF_ENABLED && FP_ENABLED && D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam, nsam
        CHECK_ASSERTION(__LINE__, size(f, 1, IK) >= 2_IK, SK_"@setACF(): The condition `1 < size(f)` must hold. size(f) = "//getStr(size(f, 1, IK)))
        CHECK_ASSERTION(__LINE__, size(f, 1, IK) == size(coef, 1, IK), SK_"@setACF(): The condition `size(f) == size(coef)` must hold. size(f), size(coef) = "//getStr([size(f, 1, IK), size(coef, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(f, 1, IK) == size(work, 1, IK), SK_"@setACF(): The condition `size(f) == size(work)` must hold. size(f), size(work) = "//getStr([size(f, 1, IK), size(work, 1, IK)]))
        call setFFTF(factor, coef, f, work, inf)
        nsam = size(f, 1, IK)
#if     CK_ENABLED
        if (inf) then
            ! fft of `f` is in `work`.
            do concurrent(isam = 1 : nsam)
                f(isam)%re = work(isam)%re**2 + work(isam)%im**2
                f(isam)%im = 0._TKC
            end do
        else
            ! fft of `f` is in `f`.
            do concurrent(isam = 1 : nsam)
                f(isam)%re = f(isam)%re**2 + f(isam)%im**2
                f(isam)%im = 0._TKC
            end do
        end if
#elif   RK_ENABLED
        if (mod(nsam, 2_IK) == 0_IK) nsam = nsam - 1_IK
        if (inf) then
            ! fft of `f` is in `work`.
            f(1) = work(1)**2
            do concurrent(isam = 2 : nsam : 2)
                f(isam) = work(isam)**2 + work(isam + 1)**2
                f(isam + 1) = 0._TKC
            end do
            if (nsam < size(f, 1, IK)) f(nsam + 1) = work(nsam + 1)**2
        else
            ! fft of `f` is in `f`.
            f(1) = f(1)**2
            do concurrent(isam = 2 : nsam : 2)
                f(isam) = f(isam)**2 + f(isam + 1)**2
                f(isam + 1) = 0._TKC
            end do
            if (nsam < size(f, 1, IK)) f(nsam + 1) = f(nsam + 1)**2
        end if
#else
#error  "Unrecognized interface."
#endif
        ! result is now in `f`.
        call setFFTR(factor, coef, f, work, inf)
        inf = .not. inf
#undef  SET_CCF

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getCCF_ENABLED && FG_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TYPE_OF_SEQ :: meanf, meang
        TYPE_OF_SEQ, parameter :: ZERO = 0._TKC
        TYPE_OF_SEQ, allocatable :: ff(:), gg(:), coef(:)
        real(TKC) :: normfac, sumsqf, sumsqg
        integer(IK), allocatable :: factor(:)
        class(*), allocatable :: norm_def
        integer(IK) :: lagmax
        integer(IK) :: lagmin
        integer(IK) :: lenccf
        integer(IK) :: iell
        logical(LK) :: inf

        lagmax = size(g, 1, IK) - 1_IK
        lagmin = 1_IK - size(f, 1, IK)
        if (present(lag)) then
            if (size(lag, 1, IK) == 0_IK) then
                allocate(ccf(0))
                return
            end if
            CHECK_ASSERTION(__LINE__, isAscending(lag), SK_"@setCCF(): The condition `isAscending(lag)` must hold. lag = "//getStr(lag))
            lagmax = max(lagmax, lag(size(lag, 1, IK)))
            lagmin = min(lagmin, lag(1))
        end if
        lenccf = lagmax - lagmin + 1_IK

        ! Allocate and pad the arrays.

        allocate(ff(lenccf), gg(lenccf), coef(lenccf), ccf(lenccf))

        ! Shift the padded arrays if needed.

        if (present(norm)) then
            norm_def = norm
        else
            norm_def = zscore_type()
        end if

        normfac = 1._TKC / lenccf
        if (same_type_as(norm_def, meanshift_type()) .or. same_type_as(norm_def, zscore_type())) then
            meanf = getMean(f)
            meang = getMean(g)
            ff(1 : size(f, 1, IK)) = f - meanf
            gg(1 : size(g, 1, IK)) = g - meang
        else
            ff(1 : size(f, 1, IK)) = f
            gg(1 : size(g, 1, IK)) = g
        end if
        ff(size(f, 1, IK) + 1 : lenccf) = ZERO
        gg(size(g, 1, IK) + 1 : lenccf) = ZERO

        if (same_type_as(norm_def, stdscale_type()) .or. same_type_as(norm_def, zscore_type())) then
            sumsqf = 0._TKC
            do iell = 1, size(f, 1, IK)
                sumsqf = sumsqf + GET_ABSQ(ff(iell))
            end do
            sumsqg = 0._TKC
            do iell = 1, size(g, 1, IK)
                sumsqg = sumsqg + GET_ABSQ(gg(iell))
            end do
            normfac = normfac / sqrt(sumsqf * sumsqg)
        elseif (.not. (same_type_as(norm_def, meanshift_type()) .or. same_type_as(norm_def, nothing_type()))) then
            error stop MODULE_NAME//SK_"@getCCF(): Unrecognized `norm`." ! LCOV_EXCL_LINE
        end if

        factor = getFactorFFT(ff, coef)
        call setCCF(factor, coef, ff, gg, ccf, inf)

        if (inf) then
            ccf = cshift(ff * normfac, lagmin)
        else
            ccf = cshift(gg * normfac, lagmin)
        end if
        if (present(lag)) ccf = ccf(lag - lagmin + 1)
        deallocate(ff, gg, coef) ! essential for gfortran on heap.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   setCCF_ENABLED && FP_ENABLED && FG_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: isam, nsam
        CHECK_ASSERTION(__LINE__, size(f, 1, IK) >= 2_IK, SK_"@setCCF(): The condition `1 < size(f)` must hold. size(f) = "//getStr(size(f, 1, IK)))
        CHECK_ASSERTION(__LINE__, size(f, 1, IK) == size(g, 1, IK), SK_"@setCCF(): The condition `size(f) == size(g)` must hold. size(f), size(g) = "//getStr([size(f, 1, IK), size(g, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(f, 1, IK) == size(coef, 1, IK), SK_"@setCCF(): The condition `size(f) == size(coef)` must hold. size(f), size(coef) = "//getStr([size(f, 1, IK), size(coef, 1, IK)]))
        CHECK_ASSERTION(__LINE__, size(f, 1, IK) == size(work, 1, IK), SK_"@setCCF(): The condition `size(f) == size(work)` must hold. size(f), size(work) = "//getStr([size(f, 1, IK), size(work, 1, IK)]))
        call setFFTF(factor, coef, f, work, inf)
        nsam = size(f, 1, IK)
#if     CK_ENABLED
#define SET_CCF(RES,TS1,TS2)\
do concurrent(isam = 1 : nsam); \
RES(isam) = TS1(isam) * conjg(TS2(isam)); \
end do;
        if (inf) then
            ! fft of `f` is in `work`.
            call setFFTF(factor, coef, g, f, inf)
            if (inf) then
                ! fft of `g` is in `f`.
                SET_CCF(f,work,f)
            else
                ! fft of `g` is in `g`.
                SET_CCF(f,work,g)
            end if
        else
            ! fft of `f` is in `f`.
            call setFFTF(factor, coef, g, work, inf)
            if (inf) then
                ! fft of `g` is in `work`.
                SET_CCF(f,f,work)
            else
                ! fft of `g` is in `g`.
                SET_CCF(f,f,g)
            end if
        end if
#elif   RK_ENABLED
#define SET_CCF(RES,TS1,TS2)\
RES(1) = TS2(1) * TS1(1); \
do concurrent(isam = 2 : nsam : 2); \
t1 = TS1(isam) * TS2(isam) + TS1(isam + 1) * TS2(isam + 1); \
t2 = TS1(isam + 1) * TS2(isam) - TS1(isam) * TS2(isam + 1); \
RES(isam + 1) = t2; \
RES(isam) = t1; \
end do; \
if (nsam < size(f, 1, IK)) f(nsam + 1) = TS2(nsam + 1) * TS1(nsam + 1);
        block
            TYPE_OF_SEQ :: t1, t2
            if (mod(nsam, 2_IK) == 0_IK) nsam = nsam - 1_IK
            if (inf) then
                ! fft of `f` is in `work`.
                call setFFTF(factor, coef, g, f, inf)
                if (inf) then
                    ! fft of `g` is in `f`.
                    SET_CCF(f,work,f)
                else
                    ! fft of `g` is in `g`.
                    SET_CCF(f,work,g)
                end if
            else
                ! fft of `f` is in `f`.
                call setFFTF(factor, coef, g, work, inf)
                if (inf) then
                    ! fft of `g` is in `work`.
                    SET_CCF(f,f,work)
                else
                    ! fft of `g` is in `g`.
                    SET_CCF(f,f,g)
                end if
            end if
        end block
#else
#error  "Unrecognized interface."
#endif
        ! result is now in `f`.
        call setFFTR(factor, coef, f, g, inf)
        inf = .not. inf
        !if (inwork) work = g
        !call setFFTR(factor, coef, work, f, inwork)
        !normfac = 1._TKC / real(size(work, 1, IK), TKC)
        !if (inwork) then
        !    do concurrent(isam = 1 : size(work, 1, IK))
        !        work(isam) = f(isam) * normfac
        !    end do
        !else
        !    do concurrent(isam = 1 : size(work, 1, IK))
        !        work(isam) = work(isam) * normfac
        !    end do
        !end if
#undef  SET_CCF

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_SEQ
#undef  GET_ABSQ
#undef  GET_RE