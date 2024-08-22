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
!>  This include file contains procedure implementations of the tests of [pm_distGenExpGamma](@ref pm_distGenExpGamma).
!>
!>  \author
!>  \AmirShahmoradi, Tuesday 2:06 AM, September 21, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_arraySpace, only: getLinSpace

        real(RKG)   , parameter     :: TOL = epsilon(0._RKG) * 100._RKG
        integer(IK)                 :: i

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGenExpGammaLogPDFNF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKG)   , allocatable   :: kappa(:), invOmega(:), logTarget(:), logTarget_ref(:), diff(:)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
        kappa = getLinSpace(0.1_RKG, 10._RKG, count = 5_IK)
        logTarget_ref = -log_gamma(kappa)
        invOmega = [1._RKG]
        allocate(logTarget, diff, mold = logTarget_ref)

        do i = 1_IK, size(kappa, kind = IK)
            logTarget(i) = getGenExpGammaLogPDFNF(kappa(i))
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The logPDFNF must be computed correctly for the given input scalar value of `shape` and default `invScale`.")
        end do

        do i = 1_IK, size(kappa, kind = IK)
            logTarget(i) = getGenExpGammaLogPDFNF(kappa(i), invOmega(1))
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The logPDFNF must be computed correctly for the given input scalar value of `shape` and `invScale = 1.`.")
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
        kappa = getLinSpace(0.1_RKG, 10._RKG, count = 5_IK)
        invOmega = [2._RKG]
        logTarget_ref = -log_gamma(kappa) + log(invOmega(1))
        allocate(logTarget, diff, mold = logTarget_ref)

        do i = 1_IK, size(kappa, kind = IK)
            logTarget(i) = getGenExpGammaLogPDFNF(kappa(i), invOmega(1))
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            call report()
            call test%assert(assertion, SK_"The logPDFNF must be computed correctly for the given input scalar value of `shape` and `invScale = 2.`.")
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
        kappa = getLinSpace(0.1_RKG, 10._RKG, count = 5_IK)
        invOmega = getLinSpace(10._RKG, 0.1_RKG, count = 5_IK)
        logTarget_ref = -log_gamma(kappa(1)) + log(invOmega)
        logTarget = getGenExpGammaLogPDFNF(kappa(1), invOmega)
        diff = abs(logTarget - logTarget_ref)

        do i = 1_IK, size(kappa, kind = IK)
            call report()
            call test%assert(assertion, SK_"The logPDFNF must be computed correctly for the given input scalar value of `shape` and vector `invScale`.")
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
        kappa = getLinSpace(0.1_RKG, 10._RKG, count = 5_IK)
        invOmega = getLinSpace(10._RKG, 0.1_RKG, count = 5_IK)
        logTarget_ref = -log_gamma(kappa) + log(invOmega)
        logTarget = getGenExpGammaLogPDFNF(kappa, invOmega)
        diff = abs(logTarget - logTarget_ref)

        do i = 1_IK, size(kappa, kind = IK)
            call report()
            call test%assert(assertion, SK_"The logPDFNF must be computed correctly for the given input vector value of `shape` and `invScale`.")
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getGenExpGammaLogPDF_ENABLED || setGenExpGammaLogPDF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     setGenExpGammaLogPDF_ENABLED
        real(RKG)                   :: kappa_current, invOmega_current, logSigma_current
#endif
        real(RKG)   , allocatable   :: logPDFNF(:), kappa(:), invOmega(:), LogSigma(:), Point(:), logTarget(:), logTarget_ref(:), diff(:)

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     getGenExpGammaLogPDF_ENABLED || setGenExpGammaLogPDF_ENABLED
        Point = real([-10,-5,-1,0,2,4,8], kind = RKG)
        logTarget_ref = [ -10.0000453999297624848515355915155600_RKG & ! LCOV_EXCL_LINE
                        , -5.00673794699908546709663604842314809_RKG & ! LCOV_EXCL_LINE
                        , -1.36787944117144232159552377016146087_RKG & ! LCOV_EXCL_LINE
                        , -1.00000000000000000000000000000000000_RKG & ! LCOV_EXCL_LINE
                        , -5.38905609893065022723042746057500802_RKG & ! LCOV_EXCL_LINE
                        , -50.5981500331442390781102612028608809_RKG & ! LCOV_EXCL_LINE
                        , -2972.95798704172827474359209945288863_RKG & ! LCOV_EXCL_LINE
                        ]
#endif
        allocate(logTarget, diff, mold = logTarget_ref)

        do i = 1_IK, size(Point, kind = IK)
#if         getGenExpGammaLogPDF_ENABLED
            logTarget(i) = getGenExpGammaLogPDF(Point(i))
#elif       setGenExpGammaLogPDF_ENABLED
            call setGenExpGammaLogPDF(logTarget(i), Point(i))
#endif
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input scalar `x` with the default parameters.", int(__LINE__, IK))
        end do

#if     getGenExpGammaLogPDF_ENABLED
        logTarget = getGenExpGammaLogPDF(Point)
#elif   setGenExpGammaLogPDF_ENABLED
        call setGenExpGammaLogPDF(logTarget, Point)
#endif
        do i = 1_IK, size(Point, kind = IK)
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input vector `x` with the default parameters.", int(__LINE__, IK))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     getGenExpGammaLogPDF_ENABLED || setGenExpGammaLogPDF_ENABLED
        kappa = [5._RKG]
        logPDFNF = getGenExpGammaLogPDFNF(kappa)
        Point = real([-10,-5,-1,0,2,4,8], kind = RKG)
        logTarget_ref = [ -53.17809923027770810449847719281262070_RKG & ! LCOV_EXCL_LINE
                        , -28.18479177734703108674357764972020190_RKG & ! LCOV_EXCL_LINE
                        , -8.545933271519387941242465371458516020_RKG & ! LCOV_EXCL_LINE
                        , -4.178053830347945619646941601297055340_RKG & ! LCOV_EXCL_LINE
                        , -0.567109929278595846877369061872063355_RKG & ! LCOV_EXCL_LINE
                        , -37.77620386349218469775720280415793240_RKG & ! LCOV_EXCL_LINE
                        , -2944.136040872076220363239041054185630_RKG & ! LCOV_EXCL_LINE
                        ]
#endif
        allocate(logTarget, diff, mold = logTarget_ref)

        do i = 1_IK, size(Point, kind = IK)
#if         getGenExpGammaLogPDF_ENABLED
            logTarget(i:i) = getGenExpGammaLogPDF(Point(i), kappa)
#elif       setGenExpGammaLogPDF_ENABLED
            call setGenExpGammaLogPDF(logTarget(i:i), Point(i), logPDFNF, kappa)
#endif
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input scalar `x` with `shape = [5]`.", int(__LINE__, IK))
        end do

#if     getGenExpGammaLogPDF_ENABLED
        logTarget = getGenExpGammaLogPDF(Point, kappa(1))
#elif   setGenExpGammaLogPDF_ENABLED
        call setGenExpGammaLogPDF(logTarget, Point, logPDFNF(1), kappa(1))
#endif
        do i = 1_IK, size(Point, kind = IK)
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input vector `x` with `shape = [5]`.", int(__LINE__, IK))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     getGenExpGammaLogPDF_ENABLED || setGenExpGammaLogPDF_ENABLED
        kappa = [5._RKG]
        invOmega = [1._RKG]
        LogSigma = [0._RKG]
        logPDFNF = getGenExpGammaLogPDFNF(kappa)
        Point = real([-10,-5,-1,0,2,4,8], kind = RKG)
        logTarget_ref = [ -53.17809923027770810449847719281262070_RKG & ! LCOV_EXCL_LINE
                        , -28.18479177734703108674357764972020190_RKG & ! LCOV_EXCL_LINE
                        , -8.545933271519387941242465371458516020_RKG & ! LCOV_EXCL_LINE
                        , -4.178053830347945619646941601297055340_RKG & ! LCOV_EXCL_LINE
                        , -0.567109929278595846877369061872063355_RKG & ! LCOV_EXCL_LINE
                        , -37.77620386349218469775720280415793240_RKG & ! LCOV_EXCL_LINE
                        , -2944.136040872076220363239041054185630_RKG & ! LCOV_EXCL_LINE
                        ]
#endif
        allocate(logTarget, diff, mold = logTarget_ref)

        do i = 1_IK, size(Point, kind = IK)
#if         getGenExpGammaLogPDF_ENABLED
            logTarget(i:i) = getGenExpGammaLogPDF(Point(i), kappa, invOmega(1), LogSigma)
#elif       setGenExpGammaLogPDF_ENABLED
            call setGenExpGammaLogPDF(logTarget(i:i), Point(i), logPDFNF, kappa, invOmega(1), LogSigma)
#endif
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input scalar `x` with `shape, invScale, loc = [5], [1], [0]`.", int(__LINE__, IK))
        end do

#if     getGenExpGammaLogPDF_ENABLED
        logTarget = getGenExpGammaLogPDF(Point, kappa(1), invOmega(1), LogSigma(1))
#elif   setGenExpGammaLogPDF_ENABLED
        call setGenExpGammaLogPDF(logTarget, Point, logPDFNF(1), kappa(1), invOmega(1), LogSigma(1))
#endif
        do i = 1_IK, size(Point, kind = IK)
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input vector `x` with `shape, invScale, loc = [5], [1], [0]`.", int(__LINE__, IK))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     getGenExpGammaLogPDF_ENABLED || setGenExpGammaLogPDF_ENABLED
        kappa = [5._RKG]
        invOmega = [1._RKG/7._RKG]
        LogSigma = [0._RKG]
        logPDFNF = getGenExpGammaLogPDFNF(kappa, invOmega)
        Point = real([-10,-5,-1,0,2,4,8], kind = RKG)
        logTarget_ref = [ -12.5064721587021775768864997445425634_RKG & ! LCOV_EXCL_LINE
                        , -9.18493421038878347841123612554047209_RKG & ! LCOV_EXCL_LINE
                        , -6.70512759343915483796951253560755200_RKG & ! LCOV_EXCL_LINE
                        , -6.12396397940325892475229434474023509_RKG & ! LCOV_EXCL_LINE
                        , -5.02610474827918033062690825274710576_RKG & ! LCOV_EXCL_LINE
                        , -4.03761607469555673334531727846312397_RKG & ! LCOV_EXCL_LINE
                        , -2.54539302868736732595888578628241702_RKG & ! LCOV_EXCL_LINE
                        ]
#endif
        allocate(logTarget, diff, mold = logTarget_ref)

        do i = 1_IK, size(Point, kind = IK)
#if         getGenExpGammaLogPDF_ENABLED
            logTarget(i:i) = getGenExpGammaLogPDF(Point(i), kappa, invOmega(1))
#elif       setGenExpGammaLogPDF_ENABLED
            call setGenExpGammaLogPDF(logTarget(i:i), Point(i), logPDFNF, kappa, invOmega(1))
#endif
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input scalar `x` with `shape, invScale = [5], [7]`.", int(__LINE__, IK))
        end do

#if     getGenExpGammaLogPDF_ENABLED
        logTarget = getGenExpGammaLogPDF(Point, kappa(1), invOmega(1), LogSigma(1))
#elif   setGenExpGammaLogPDF_ENABLED
        call setGenExpGammaLogPDF(logTarget, Point, logPDFNF(1), kappa(1), invOmega(1), LogSigma(1))
#endif
        do i = 1_IK, size(Point, kind = IK)
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input vector `x` with `shape, invScale = [5], [1/7]`.", int(__LINE__, IK))
        end do

#if     getGenExpGammaLogPDF_ENABLED
        logTarget = getGenExpGammaLogPDF(Point, kappa(1), invOmega(1), LogSigma(1))
#elif   setGenExpGammaLogPDF_ENABLED
        call setGenExpGammaLogPDF(logTarget, Point, logPDFNF(1), kappa(1), invOmega(1), LogSigma(1))
#endif
        do i = 1_IK, size(Point, kind = IK)
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input vector `x` with `shape, invScale, loc = [5], [1/7], [0]`.", int(__LINE__, IK))
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

#if     getGenExpGammaLogPDF_ENABLED || setGenExpGammaLogPDF_ENABLED
        kappa = [0.5_RKG]
        invOmega = [8._RKG]
        LogSigma = [3._RKG]
        logPDFNF = getGenExpGammaLogPDFNF(kappa, invOmega)
        Point = real([-10,-5,-1,0,2,4,8], kind = RKG)
        logTarget_ref = getGenExpGammaLogPDF(Point - LogSigma(1), kappa(1), invOmega(1))
#endif
        allocate(logTarget, diff, mold = logTarget_ref)

        do i = 1_IK, size(Point, kind = IK)
#if         getGenExpGammaLogPDF_ENABLED
            logTarget(i:i) = getGenExpGammaLogPDF(Point(i), kappa, invOmega(1), LogSigma(1))
#elif       setGenExpGammaLogPDF_ENABLED
            call setGenExpGammaLogPDF(logTarget(i:i), Point(i), logPDFNF, kappa, invOmega(1), LogSigma(1))
#endif
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input scalar `x` with `shape, invScale = [5], [7]`.", int(__LINE__, IK))
        end do

#if     getGenExpGammaLogPDF_ENABLED
        logTarget = getGenExpGammaLogPDF(Point, kappa(1), invOmega(1), LogSigma(1))
#elif   setGenExpGammaLogPDF_ENABLED
        call setGenExpGammaLogPDF(logTarget, Point, logPDFNF(1), kappa(1), invOmega(1), LogSigma(1))
#endif
        do i = 1_IK, size(Point, kind = IK)
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input vector `x` with `shape, invScale = [5], [1/7]`.", int(__LINE__, IK))
        end do

#if     getGenExpGammaLogPDF_ENABLED
        logTarget = getGenExpGammaLogPDF(Point, kappa(1), invOmega(1), LogSigma(1))
#elif   setGenExpGammaLogPDF_ENABLED
        call setGenExpGammaLogPDF(logTarget, Point, logPDFNF(1), kappa(1), invOmega(1), LogSigma(1))
#endif
        do i = 1_IK, size(Point, kind = IK)
            diff(i) = abs(logTarget(i) - logTarget_ref(i))
            assertion = assertion .and. diff(i) <= TOL
            call report()
            call test%assert(assertion, SK_"The `logPDF` must be computed correctly for a input vector `x` with `shape, invScale, loc = [5], [1/7], [0]`.", int(__LINE__, IK))
        end do

#if     setGenExpGammaLogPDF_ENABLED
        block

            use pm_val2str, only: getStr
            use pm_quadPack, only: isFailedQuad
            use pm_distUnif, only: setUnifRand
            use pm_distGenExpGamma, only: getGenExpGammaCDF

            real(RKG) :: lb, ub, integral_def, integral, abserr
            lb = -10._RKG
            ub = +10._RKG

            integral_def = getGenExpGammaCDF(ub) - getGenExpGammaCDF(lb)
            assertion = .not. isFailedQuad(getFunc1, lb, ub, integral, abserr)
            call test%assert(assertion, SK_"The integral of the PDF over its support must not fail.", int(__LINE__, IK))
            assertion = assertion .and. abs(integral - integral_def) <= abserr
            call test%assert(assertion, SK_"The integral of the PDF must equal unity. integral, abserr = "//getStr([integral, abserr]), int(__LINE__, IK))

            call setUnifRand(kappa_current, epsilon(0._RKG), 5._RKG)
            integral_def = getGenExpGammaCDF(ub, kappa = kappa_current) - getGenExpGammaCDF(lb, kappa = kappa_current)
            assertion = .not. isFailedQuad(getFunc2, lb, ub, integral, abserr)
            call test%assert(assertion, SK_"The integral of the PDF over its support must not fail.", int(__LINE__, IK))
            assertion = assertion .and. abs(integral - integral_def) <= abserr
            call test%assert(assertion, SK_"The integral of the PDF must equal unity. integral, abserr, kappa_current = "//getStr([integral, abserr, kappa_current]), int(__LINE__, IK))

            call setUnifRand(kappa_current, epsilon(0._RKG), 5._RKG)
            call setUnifRand(invOmega_current, epsilon(0._RKG), 5._RKG)
            integral_def = getGenExpGammaCDF(ub, kappa = kappa_current, invOmega = invOmega_current) - getGenExpGammaCDF(lb, kappa = kappa_current, invOmega = invOmega_current)
            assertion = .not. isFailedQuad(getFunc3, lb, ub, integral, abserr)
            call test%assert(assertion, SK_"The integral of the PDF over its support must not fail.", int(__LINE__, IK))
            assertion = assertion .and. abs(integral - integral_def) <= abserr
            call test%assert(assertion, SK_"The integral of the PDF must equal unity. integral, abserr, kappa_current, invOmega_current = "//getStr([integral, abserr, kappa_current, invOmega_current]), int(__LINE__, IK))

            call setUnifRand(kappa_current, epsilon(0._RKG), 5._RKG)
            call setUnifRand(invOmega_current, epsilon(0._RKG), 5._RKG)
            call setUnifRand(logSigma_current, -10._RKG, 10._RKG)
            integral_def = getGenExpGammaCDF(ub, kappa_current, invOmega_current, logSigma_current) - getGenExpGammaCDF(lb, kappa_current, invOmega_current, logSigma_current)
            assertion = .not. isFailedQuad(getFunc4, lb, ub, integral, abserr)
            call test%assert(assertion, SK_"The integral of the PDF over its support must not fail.", int(__LINE__, IK))
            assertion = assertion .and. abs(integral - integral_def) <= abserr
            call test%assert(assertion, SK_"The integral of the PDF must equal unity. integral, abserr, kappa_current, invOmega_current, logSigma_current = "//getStr([integral, abserr, kappa_current, invOmega_current, logSigma_current]), int(__LINE__, IK))

        end block
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

#if     setGenExpGammaLogPDF_ENABLED

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getFunc1(x) result(func)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: func
            call setGenExpGammaLogPDF(func, x)
            func = exp(func)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getFunc2(x) result(func)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: func
            call setGenExpGammaLogPDF(func, x, getGenExpGammaLogPDFNF(kappa_current), kappa_current)
            func = exp(func)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getFunc3(x) result(func)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: func
            call setGenExpGammaLogPDF(func, x, getGenExpGammaLogPDFNF(kappa_current, invOmega_current), kappa_current, invOmega_current)
            func = exp(func)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function getFunc4(x) result(func)
            real(RKG)   , intent(in)    :: x
            real(RKG)                   :: func
            call setGenExpGammaLogPDF(func, x, getGenExpGammaLogPDFNF(kappa_current, invOmega_current), kappa_current, invOmega_current, logSigma_current)
            func = exp(func)
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif
        subroutine reset()
            if (allocated(diff)) deallocate(diff)
            if (allocated(logTarget)) deallocate(logTarget)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                assertion = assertion .and. diff(i) <= TOL
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "kappa          ", kappa
                write(test%disp%unit,"(*(g0,:,', '))") "invOmega       ", invOmega
                write(test%disp%unit,"(*(g0,:,', '))") "logTarget_ref  ", logTarget_ref
                write(test%disp%unit,"(*(g0,:,', '))") "logTarget      ", logTarget
                write(test%disp%unit,"(*(g0,:,', '))") "diff           ", diff
                write(test%disp%unit,"(*(g0,:,', '))") "TOL            ", TOL
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
