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
!>  This file contains procedure implementations of tests of [test_pm_sampleCCF](@ref test_pm_sampleCCF).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 5:03 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(:, SK), allocatable :: format
        ! Define the comparison precision and tolerance.
        real(TKC), parameter :: rtol = epsilon(1._TKC) * 5000
#if     CK_ENABLED
#define TYPE_OF_SEQ complex(TKC)
#define GET_CONJG(X)conjg(X)
#define COERCE(X)X
        complex(TKC), parameter :: ZERO = (0._TKC, 0._TKC), ONES = (1._TKC, 1._TKC), TWOS = (2._TKC, 2._TKC)
        complex(TKC), parameter :: ctol = (rtol, rtol)
#elif   RK_ENABLED
#define TYPE_OF_SEQ real(TKC)
#define COERCE(X)real(X)
#define GET_CONJG(X)X
        real(TKC), parameter :: ZERO = 0._TKC, ONES = 1._TKC, TWOS = 2._TKC
        real(TKC), parameter :: ctol = rtol
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%
#if     getCCF_ENABLED
        !%%%%%%%%%%%%%

        TYPE_OF_SEQ, allocatable :: f(:), g(:), ccf(:), ccf_ref(:), diff(:)
        integer(IK), allocatable :: lag(:)
        integer(IK) :: nf, ng, itry
#if     CK_ENABLED
        format = getFormat(deliml = '', subsep = SK_'', delimr = 'i', subcount = 2_IK, signed = .true._LK)
#elif   RK_ENABLED
        format = getFormat(mold = [0._TKC])
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK
        do itry = 1, 100

            ! test FG pair interface.

            block

                nf = getUnifRand(2, 40)
                ng = getUnifRand(2, 40)
                f = getUnifRand(ONES, TWOS * 3 / 2, nf)
                g = getUnifRand(ONES, TWOS * 3 / 2, ng)

                ccf = getCCF(f, g)
                ccf_ref = getCCF_ref(f, g)
                call report(__LINE__)

                ccf = getCCF(f, g, norm = zscore)
                ccf_ref = getCCF_ref(f, g, norm = zscore)
                call report(__LINE__, norm = zscore)

                ccf = getCCF(f, g, norm = nothing)
                ccf_ref = getCCF_ref(f, g, norm = nothing)
                call report(__LINE__, norm = nothing)

                ccf = getCCF(f, g, norm = stdscale)
                ccf_ref = getCCF_ref(f, g, norm = stdscale)
                call report(__LINE__, norm = stdscale)

                ccf = getCCF(f, g, norm = meanshift)
                ccf_ref = getCCF_ref(f, g, norm = meanshift)
                call report(__LINE__, norm = meanshift)

            end block

            ! test FG pair interface with arbitrary `lag`.

            block

                nf = getUnifRand(2, 40)
                ng = getUnifRand(2, 40)
                f = getUnifRand(ONES, TWOS * 3 / 2, nf)
                g = getUnifRand(ONES, TWOS * 3 / 2, ng)

                lag = getRange(1_IK, 50_IK)

                ccf = getCCF(f, g, lag)
                ccf_ref = getCCF_ref(f, g, lag)
                call report(__LINE__)

                ccf = getCCF(f, g, lag, norm = zscore)
                ccf_ref = getCCF_ref(f, g, lag, norm = zscore)
                call report(__LINE__, norm = zscore)

                ccf = getCCF(f, g, lag, norm = nothing)
                ccf_ref = getCCF_ref(f, g, lag, norm = nothing)
                call report(__LINE__, norm = nothing)

                ccf = getCCF(f, g, lag, norm = stdscale)
                ccf_ref = getCCF_ref(f, g, lag, norm = stdscale)
                call report(__LINE__, norm = stdscale)

                ccf = getCCF(f, g, lag, norm = meanshift)
                ccf_ref = getCCF_ref(f, g, lag, norm = meanshift)
                call report(__LINE__, norm = meanshift)

                lag = [integer(IK) ::]
                ccf = getCCF(f, g, lag)
                ccf_ref = [TYPE_OF_SEQ ::]
                call report(__LINE__, lag)

            end block

        end do

    contains

        function getCCF_ref(f, g, lag, norm) result(ccf)
            TYPE_OF_SEQ, intent(in), contiguous :: f(:), g(:)
            integer(IK), intent(in), contiguous, optional :: lag(:)
            TYPE_OF_SEQ, allocatable :: ff(:), gg(:), coef(:), ccf(:)
            class(*), intent(in), optional :: norm
            integer(IK), allocatable :: lag_def(:)
            integer(IK), allocatable :: factor(:)
            class(*), allocatable :: norm_def
            real(TKC) :: dumf, dumg, normfac
            TYPE_OF_SEQ :: meanf, meang
            integer(IK) :: lagmin, lagmax
            integer(IK) :: lenccf, lenf, leng
            logical(LK) :: inf
            lenf = size(f, 1, IK)
            leng = size(g, 1, IK)
            lagmax = leng - 1_IK
            lagmin = 1_IK - lenf
            if (present(lag)) then
                if (size(lag, 1, IK) == 0_IK) then
                    call setResized(ccf, 0_IK)
                    return
                end if
                lagmax = max(lagmax, lag(size(lag, 1, IK)))
                lagmin = min(lagmin, lag(1))
                lag_def = lag
            else
                lag_def = getRange(lagmin, lagmax)
            end if
            lenccf = lagmax - lagmin + 1_IK
            call setResized(ff, lenccf)
            call setResized(gg, lenccf)
            call setResized(ccf, lenccf)
            call setResized(coef, lenccf)
            if (present(norm)) then
                norm_def = norm
            else
                norm_def = zscore
            end if
            if (same_type_as(norm_def, meanshift) .or. same_type_as(norm_def, zscore)) then
                meanf = getMean(f)
                meang = getMean(g)
                ff(1 : lenf) = f - meanf
                gg(1 : leng) = g - meang
            else
                ff(1 : lenf) = f
                gg(1 : leng) = g
            end if
            ff(lenf + 1 :) = ZERO
            gg(leng + 1 :) = ZERO
            normfac = 1._TKC / lenccf
            if (same_type_as(norm_def, stdscale) .or. same_type_as(norm_def, zscore)) then
                dumf = real(dot_product(ff(1 : lenf), ff(1 : lenf)), TKC)
                dumg = real(dot_product(gg(1 : leng), gg(1 : leng)), TKC)
                normfac = normfac / sqrt(dumf * dumg)
            end if
            factor = getFactorFFT(ff, coef)
            call setCCF(factor, coef, ff, gg, ccf, inf)
            if (inf) then
                ccf = ff * normfac
            else
                ccf = gg * normfac
            end if
            ccf = cshift(ccf, lagmin)
            ccf = ccf(lag_def - lagmin + 1)
        end function

        subroutine report(line, lag, norm)
            integer, intent(in) :: line
            class(*), intent(in), optional :: norm
            integer(IK), intent(in), optional :: lag(:)
            if (0_IK < size(ccf, 1, IK)) then
                diff = abs(ccf - ccf_ref)
                assertion = assertion .and. all(diff < ctol)
            else
                assertion = assertion .and. size(ccf_ref, 1, IK) == 0_IK
            end if
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[nf, ng, size(ccf, 1, IK)]")
                call test%disp%show( [nf, ng, size(ccf, 1, IK)] )
                call test%disp%show("f")
                call test%disp%show( f , format = format )
                call test%disp%show("g")
                call test%disp%show( g , format = format )
                call test%disp%show("ccf")
                call test%disp%show( ccf , format = format )
                call test%disp%show("ccf_ref")
                call test%disp%show( ccf_ref , format = format )
                call test%disp%show("diff")
                call test%disp%show( diff , format = format )
                call test%disp%show("present(lag)")
                call test%disp%show( present(lag) )
                if (present(lag)) then
                    call test%disp%show("lag")
                    call test%disp%show( lag )
                end if
                call test%disp%show("present(norm)")
                call test%disp%show( present(norm) )
                if (present(norm)) then
                    call test%disp%show("[same_type_as(norm, meanshift), same_type_as(norm, stdscale), same_type_as(norm, zscore), same_type_as(norm, nothing)]")
                    call test%disp%show( [same_type_as(norm, meanshift), same_type_as(norm, stdscale), same_type_as(norm, zscore), same_type_as(norm, nothing)] )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `ccf` must be correctly computed for the specified sequences.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%
#elif   setCCF_ENABLED
        !%%%%%%%%%%%%%

        TYPE_OF_SEQ, allocatable :: f(:), g(:), ff(:), gg(:), coef(:), ccf(:), ccf_ref(:), diff(:)
        integer(IK), allocatable :: factor(:)
        integer(IK) :: nf, ng, nc, itry
        logical(LK) :: inf
#if     CK_ENABLED
        format = getFormat(deliml = '', subsep = SK_'', delimr = 'i', subcount = 2_IK, signed = .true._LK)
#elif   RK_ENABLED
        format = getFormat(mold = [0._TKC])
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK
        do itry = 2, 42, 2

            ! test FG pair interface.

            block
                nf = itry
                ng = itry + 1
                nc = nf + ng - 1_IK

                call setResized(f, nc)
                call setResized(g, nc)
                f(1:nf) = getUnifRand(ONES, TWOS * 3 / 2, nf)
                g(1:ng) = getUnifRand(ONES, TWOS * 3 / 2, ng)
                f(nf + 1:) = ZERO
                g(ng + 1:) = ZERO
                ccf_ref = getCCF_ref(f, g)

                ff = f
                gg = g
                call setResized(ccf, nc)
                call setResized(coef, nc)
                factor = getFactorFFT(f, coef)
                call setCCF(factor, coef, ff, gg, ccf, inf)
                ccf = merge(ff / nc, gg / nc, inf)
                call report(__LINE__)
            end block

        end do

    contains

        function getCCF_ref(f, g) result(ccf)
            TYPE_OF_SEQ, intent(in) :: f(:), g(:)
            complex(TKC), allocatable :: ff(:), gg(:)
            TYPE_OF_SEQ, allocatable :: ccf(:)
            ff = f; gg = g; ccf = COERCE(getFFTI(getFFTF(ff) * conjg(getFFTF(gg))))
        end function

        subroutine report(line)
            integer, intent(in) :: line
            diff = abs(ccf - ccf_ref)
            assertion = assertion .and. all(diff < ctol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[nf, ng, nc]")
                call test%disp%show( [nf, ng, nc] )
                call test%disp%show("f")
                call test%disp%show( f , format = format )
                call test%disp%show("g")
                call test%disp%show( g , format = format )
                call test%disp%show("ccf")
                call test%disp%show( ccf , format = format )
                call test%disp%show("ccf_ref")
                call test%disp%show( ccf_ref , format = format )
                call test%disp%show("diff")
                call test%disp%show( diff , format = format )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `ccf` must be correctly computed for the specified sequences.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%
#elif   getACF_ENABLED
        !%%%%%%%%%%%%%

        TYPE_OF_SEQ, allocatable :: f(:), g(:), acf(:), acf_ref(:), diff(:)
        integer(IK), allocatable :: lag(:)
        integer(IK) :: nf, itry
#if     CK_ENABLED
        format = getFormat(deliml = '', subsep = SK_'', delimr = 'i', subcount = 2_IK, signed = .true._LK)
#elif   RK_ENABLED
        format = getFormat(mold = [0._TKC])
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK

        do itry = 1, 100

            ! test FG pair interface.

            nf = getUnifRand(2, 40)
            f = getUnifRand(ONES, TWOS * 3 / 2, nf)

            call runTestWith(__LINE__)
            call runTestWith(__LINE__, norm = zscore)
            call runTestWith(__LINE__, norm = stdscale)
            call runTestWith(__LINE__, norm = meanshift)

            lag = getRange(1_IK, 50_IK)
            call runTestWith(__LINE__, lag)
            call runTestWith(__LINE__, lag, norm = zscore)
            call runTestWith(__LINE__, lag, norm = stdscale)
            call runTestWith(__LINE__, lag, norm = meanshift)

            lag = [integer(IK) ::]
            call runTestWith(__LINE__, lag)
            call runTestWith(__LINE__, lag, norm = zscore)
            call runTestWith(__LINE__, lag, norm = stdscale)
            call runTestWith(__LINE__, lag, norm = meanshift)

        end do

    contains

        subroutine runTestWith(line, lag, norm)
            integer, intent(in) :: line
            class(*), intent(in), optional :: norm
            integer(IK), intent(in), contiguous, optional :: lag(:)
            integer(IK), allocatable :: lag_def(:)
            acf = getACF(f, lag, norm)
            g = f
            if (present(lag)) then
                acf_ref = getCCF(f, g, lag, norm)
            else
                lag_def = getRange(0_IK, size(f, 1, IK) - 1_IK)
                acf_ref = getCCF(f, g, lag_def, norm)
            end if
            call report(line, lag, norm)
        end subroutine

        subroutine report(line, lag, norm)
            integer, intent(in) :: line
            class(*), intent(in), optional :: norm
            integer(IK), intent(in), optional :: lag(:)
            if (0_IK < size(acf, 1, IK)) then
                diff = abs(acf - acf_ref)
                assertion = assertion .and. all(diff < ctol)
            else
                assertion = assertion .and. size(acf_ref, 1, IK) == 0_IK
            end if
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[nf, size(acf, 1, IK)]")
                call test%disp%show( [nf, size(acf, 1, IK)] )
                call test%disp%show("f")
                call test%disp%show( f , format = format )
                call test%disp%show("acf")
                call test%disp%show( acf , format = format )
                call test%disp%show("acf_ref")
                call test%disp%show( acf_ref , format = format )
                call test%disp%show("diff")
                call test%disp%show( diff , format = format )
                call test%disp%show("present(lag)")
                call test%disp%show( present(lag) )
                if (present(lag)) then
                    call test%disp%show("lag")
                    call test%disp%show( lag )
                end if
                call test%disp%show("present(norm)")
                call test%disp%show( present(norm) )
                if (present(norm)) then
                    call test%disp%show("[same_type_as(norm, meanshift), same_type_as(norm, stdscale), same_type_as(norm, zscore), same_type_as(norm, nothing)]")
                    call test%disp%show( [same_type_as(norm, meanshift), same_type_as(norm, stdscale), same_type_as(norm, zscore), same_type_as(norm, nothing)] )
                end if
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `acf` must be correctly computed for the specified sequences.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%
#elif   setACF_ENABLED
        !%%%%%%%%%%%%%

        TYPE_OF_SEQ, allocatable :: seq(:), f(:), g(:), coef(:), acf(:), acf_ref(:), diff(:)
        integer(IK), allocatable :: factor(:)
        integer(IK) :: nf, nc, itry
        logical(LK) :: inf
#if     CK_ENABLED
        format = getFormat(deliml = '', subsep = SK_'', delimr = 'i', subcount = 2_IK, signed = .true._LK)
#elif   RK_ENABLED
        format = getFormat(mold = [0._TKC])
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK
        do itry = 2, 42, 2

            ! test FG pair interface.

            block
                nf = itry
                nc = nf * 2_IK - 1_IK
                call setResized(seq, nc)
                seq(1:nf) = getUnifRand(ONES, TWOS * 3 / 2, nf)
                seq(nf + 1:) = ZERO
                f = seq
                g = seq

                call setResized(acf, nc)
                call setResized(coef, nc)
                call setResized(acf_ref, nc)
                factor = getFactorFFT(f, coef)
                call setCCF(factor, coef, f, g, acf_ref, inf)
                acf_ref = merge(f, g, inf)
                f = seq
                call setACF(factor, coef, f, acf, inf)
                acf = merge(f, acf, inf)
                call report(__LINE__)
            end block

        end do

    contains

        subroutine report(line)
            integer, intent(in) :: line
            diff = abs(acf - acf_ref)
            assertion = assertion .and. all(diff < ctol)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("[nf, nc]")
                call test%disp%show( [nf, nc] )
                call test%disp%show("seq")
                call test%disp%show( seq , format = format )
                call test%disp%show("acf")
                call test%disp%show( acf , format = format )
                call test%disp%show("acf_ref")
                call test%disp%show( acf_ref , format = format )
                call test%disp%show("diff")
                call test%disp%show( diff , format = format )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The `acf` must be correctly computed for the specified sequence.", int(line, IK))
        end subroutine

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif
#undef  TYPE_OF_SEQ
#undef  GET_CONJG
#undef  COERCE